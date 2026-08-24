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


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
