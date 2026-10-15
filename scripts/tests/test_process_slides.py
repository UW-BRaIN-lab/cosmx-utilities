#!/usr/bin/env python3
"""Tests for scripts/process-slides.py discovery and preflight logic (no AWS).

Runnable either under pytest or directly:  uv run python scripts/tests/test_process_slides.py
"""
import importlib.util
import json
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / "process-slides.py"
_spec = importlib.util.spec_from_file_location("process_slides", _SCRIPT)
ps = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ps)


IMAGE = "ghcr.io/uw-brain-lab/cosmx-utilities:headless-latest"


def _task_def(image=IMAGE, task_role="arn:aws:iam::1:role/R",
              exec_role="arn:aws:iam::1:role/E", ephemeral=200):
    return {
        "ephemeralStorage": {"sizeInGiB": ephemeral},
        "taskRoleArn": task_role,
        "executionRoleArn": exec_role,
        "containerDefinitions": [{"name": "process-slide", "image": image}],
    }


def test_pinned_fields_reads_container_and_top_level():
    """Pinned fields live at both levels of a task definition."""
    got = ps._pinned_fields(_task_def())
    assert got["image"] == IMAGE
    assert got["taskRoleArn"] == "arn:aws:iam::1:role/R"
    assert got["ephemeralStorage"] == {"sizeInGiB": 200}


def test_no_drift_when_registration_matches(monkeypatch=None):
    """An in-sync registration must not block a launch."""
    ps._template_task_definition = lambda: _task_def()
    ps._get_ecs = lambda: type("E", (), {
        "describe_task_definition": staticmethod(
            lambda **kw: {"taskDefinition": _task_def()})})()
    assert ps.task_definition_drift() == []


def test_detects_stale_image():
    """The real failure: registration left pointing at a pre-rename image, so
    every task dies on an image pull that says nothing about drift."""
    stale = "ghcr.io/keene-lab/cosmx-utilities:headless-latest"
    ps._template_task_definition = lambda: _task_def()
    ps._get_ecs = lambda: type("E", (), {
        "describe_task_definition": staticmethod(
            lambda **kw: {"taskDefinition": _task_def(image=stale)})})()
    drift = ps.task_definition_drift()
    assert [d[0] for d in drift] == ["image"], drift
    assert drift[0][1] == stale and drift[0][2] == IMAGE


def test_detects_role_and_storage_drift():
    """Roles and ephemeral storage are pinned too — a run silently under-sized
    on disk fails late and expensively."""
    ps._template_task_definition = lambda: _task_def()
    ps._get_ecs = lambda: type("E", (), {
        "describe_task_definition": staticmethod(
            lambda **kw: {"taskDefinition": _task_def(
                task_role="arn:aws:iam::1:role/WRONG", ephemeral=21)})})()
    fields = {d[0] for d in ps.task_definition_drift()}
    assert fields == {"taskRoleArn", "ephemeralStorage"}, fields


def test_unreadable_registration_does_not_block():
    """If ECS cannot be queried, warn but let the launch proceed — the check is
    a guard rail, not a new way for the pipeline to refuse to run."""
    ps._template_task_definition = lambda: _task_def()
    def boom():
        raise RuntimeError("no credentials")
    ps._get_ecs = lambda: type("E", (), {
        "describe_task_definition": staticmethod(lambda **kw: boom())})()
    assert ps.task_definition_drift() == []


def test_missing_template_does_not_block():
    """A checkout without the template file still launches."""
    ps._template_task_definition = lambda: None
    assert ps.task_definition_drift() == []


# ---- Trustworthy re-runs: versioned marker + replacing upload ---------------

class _FakeS3Marker:
    """S3 stand-in serving one slide's _SUCCESS body (or absence)."""

    def __init__(self, body):
        self.body = body

    def get_object(self, Bucket, Key):
        if self.body is None:
            raise RuntimeError("NoSuchKey")
        import io
        return {"Body": io.BytesIO(self.body)}


def _slide():
    return ps.Slide(bucket="b", base_path="Study/Run/DecodedFiles/Slide/Scan")


def _with_marker(body):
    ps._get_s3 = lambda: _FakeS3Marker(body)
    return _slide()


CURRENT = "sha256:aaaa"
OLDER = "sha256:bbbb"


def _marker_body(version):
    return json.dumps({"code_version": version, "completed_at": "t"}).encode()


def test_skip_when_same_code_produced_the_output():
    """The normal case: output made by the code about to run is genuinely done."""
    slide = _with_marker(_marker_body(CURRENT))
    assert ps.is_already_processed(slide, CURRENT) is True


def test_do_not_skip_when_different_code_produced_the_output():
    """The failure this exists for: a run finished and wrote _SUCCESS, but the
    data was wrong. After fixing the code, the stale marker must not protect the
    bad output from being replaced."""
    slide = _with_marker(_marker_body(OLDER))
    assert ps.is_already_processed(slide, CURRENT) is False


def test_missing_marker_is_not_processed():
    slide = _with_marker(None)
    assert ps.is_already_processed(slide, CURRENT) is False


def test_legacy_zero_byte_marker_still_skips():
    """Markers written before this change carry no version. Falling back to
    presence keeps existing outputs from being re-stitched wholesale."""
    slide = _with_marker(b"")
    assert ps.is_already_processed(slide, CURRENT) is True


def test_unresolvable_expected_version_falls_back_to_presence():
    """If the running version cannot be determined, behave as before rather than
    re-running everything."""
    slide = _with_marker(_marker_body(OLDER))
    assert ps.is_already_processed(slide, "") is True


def test_unparseable_marker_treated_as_versionless():
    slide = _with_marker(b"not json at all")
    assert ps.is_already_processed(slide, CURRENT) is True


def test_registry_digest_ignores_non_ghcr_and_malformed():
    """Resolution must degrade to "unknown" rather than raising mid-launch."""
    assert ps._registry_digest("docker.io/library/x:1") == ""
    assert ps._registry_digest("bare-name") == ""


# ---- Auto-sizing by transcript-file size -------------------------------------

class _FakeS3Head:
    """S3 stand-in answering head_object with a fixed size, or failing."""

    def __init__(self, size=None, error=None):
        self.size = size
        self.error = error
        self.requested = []

    def head_object(self, Bucket, Key):
        self.requested.append((Bucket, Key))
        if self.error:
            raise self.error
        return {"ContentLength": self.size}


def _sized(size=None, error=None):
    fake = _FakeS3Head(size=size, error=error)
    ps._get_s3 = lambda: fake
    return fake


def _maddie_slide():
    return ps.Slide(
        bucket="keene-cosmx-data",
        base_path=("CosMx-Maddie/20261002_SORL1_3D_40donors/DecodedFiles/"
                   "7678A77484A67316A66802A6/20260903_213537_S1"))


def test_tx_file_key_is_the_runs_flatfile_for_this_slide():
    key = _maddie_slide().tx_file_key
    print(f"  tx_file_key = {key}")
    assert key == ("CosMx-Maddie/20261002_SORL1_3D_40donors/flatFiles/"
                   "7678A77484A67316A66802A6/7678A77484A67316A66802A6_tx_file.csv.gz")


def test_tx_file_key_handles_an_experiment_directory_level():
    """Slides nested under an experiment directory keep their run's flatFiles."""
    slide = ps.Slide(bucket="b", base_path="Study/Exp/Run/DecodedFiles/S/Scan")
    assert slide.tx_file_key == "Study/Exp/Run/flatFiles/S/S_tx_file.csv.gz"


def test_the_two_slides_that_ran_out_of_memory_get_the_large_profile():
    """The real failures: 3.98 GB and 3.12 GB transcript files, killed at 4 vCPU."""
    for size in (3_978_455_480, 3_121_813_666):
        fake = _sized(size)
        got = ps.auto_fargate_size(_maddie_slide())
        print(f"  tx_file={size / 1e9:.2f} GB -> {got}")
        assert got == (ps.LARGE_SLIDE_CPU, ps.LARGE_SLIDE_MEMORY)
        assert fake.requested == [("keene-cosmx-data", _maddie_slide().tx_file_key)]


def test_small_slides_that_passed_keep_the_default():
    for size in (2_198_341_860, 2_053_511_023, 967_873_281):
        _sized(size)
        got = ps.auto_fargate_size(_maddie_slide())
        print(f"  tx_file={size / 1e9:.2f} GB -> {got}")
        assert got == (None, None)


def test_largest_slide_that_passed_is_deliberately_over_provisioned():
    """2.75 GB finished at 4 vCPU / 30 GB, but it is near the cliff, so it is
    given the large profile: an unneeded bigger task costs less than a failed one."""
    _sized(2_747_701_342)
    assert ps.auto_fargate_size(_maddie_slide()) == (
        ps.LARGE_SLIDE_CPU, ps.LARGE_SLIDE_MEMORY)


def test_threshold_sits_below_the_largest_passing_slide():
    """Near the boundary the cost of over-provisioning is small and a failure is
    not, so the threshold must not exceed the largest slide known to pass."""
    assert ps.LARGE_SLIDE_TX_FILE_BYTES <= 2_747_701_342
    assert ps.LARGE_SLIDE_TX_FILE_BYTES < 3_121_813_666


def test_boundary_is_inclusive():
    _sized(ps.LARGE_SLIDE_TX_FILE_BYTES)
    assert ps.auto_fargate_size(_maddie_slide())[0] == ps.LARGE_SLIDE_CPU
    _sized(ps.LARGE_SLIDE_TX_FILE_BYTES - 1)
    assert ps.auto_fargate_size(_maddie_slide()) == (None, None)


def test_unreadable_transcript_file_keeps_the_default():
    """A missing or denied tx file must not block the launch or guess a size."""
    _sized(error=RuntimeError("404"))
    assert ps.tx_file_bytes(_maddie_slide()) is None
    assert ps.auto_fargate_size(_maddie_slide()) == (None, None)


def test_large_profile_is_the_largest_benchmark_profile():
    assert (ps.LARGE_SLIDE_CPU, ps.LARGE_SLIDE_MEMORY) == ("16384", "122880")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
