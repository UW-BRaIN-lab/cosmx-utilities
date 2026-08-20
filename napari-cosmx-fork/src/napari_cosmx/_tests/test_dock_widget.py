import os
import pytest
from qtpy.QtWidgets import QGroupBox
from pathlib import Path

from napari_cosmx._dock_widget import GeminiQWidget, LayerGroupManagerWidget, make_stitching_widget


def _group_titles(widget):
    return {group_box.title() for group_box in widget.findChildren(QGroupBox)}


def test_make_stitching_widget_standalone(make_napari_viewer):
    viewer = make_napari_viewer()
    widget = make_stitching_widget(viewer)

    assert isinstance(widget, GeminiQWidget)
    assert widget.gem is None
    assert widget.show_stitching_widget is True
    assert widget.show_data_widgets is False

    titles = _group_titles(widget)
    assert 'Stitch Images' in titles
    assert 'Dataset' not in titles
    assert 'Morphology Images' not in titles
    assert 'RNA Transcripts' not in titles
    assert 'Protein Expression' not in titles
    assert 'Color Cells' not in titles


def test_make_stitching_widget_uses_current_viewer_when_not_injected(make_napari_viewer, monkeypatch):
    viewer = make_napari_viewer()
    monkeypatch.setattr('napari_cosmx._dock_widget.napari.current_viewer', lambda: viewer)

    widget = make_stitching_widget()

    assert isinstance(widget, GeminiQWidget)
    assert widget.viewer is viewer
    assert widget.gem is None
    assert widget.show_stitching_widget is True
    assert widget.show_data_widgets is False


def test_widget_requires_gem_when_data_widgets_enabled(make_napari_viewer):
    viewer = make_napari_viewer()
    with pytest.raises(AssertionError, match='gem is required'):
        GeminiQWidget(viewer, gem=None, show_data_widgets=True, show_stitching_widget=False)


def test_layer_group_manager_init_without_renamed_event(qtbot):
    class _Signal:
        def connect(self, _callback):
            return None

    class _Events:
        def __init__(self):
            self.removed = _Signal()
            # intentionally no `renamed` attribute to mimic older napari

    class _Layers:
        def __init__(self):
            self.events = _Events()

        def __iter__(self):
            return iter(())

    class _Viewer:
        def __init__(self):
            self.layers = _Layers()

    widget = LayerGroupManagerWidget(_Viewer())
    qtbot.addWidget(widget)
    assert widget is not None


def test_candidate_primary_folders_from_rna_path(make_napari_viewer, tmp_path):
    viewer = make_napari_viewer()
    widget = GeminiQWidget(viewer, gem=None, show_data_widgets=False, show_stitching_widget=True)
    # Build a real OS-native path so os.path.dirname splits correctly on both platforms.
    rna_path = tmp_path / "slide" / "AnalysisResults" / "RunA"
    rna_path.mkdir(parents=True)
    candidates = widget._candidate_primary_folders(str(rna_path))
    candidate_names = [os.path.basename(c) for c in candidates]
    assert candidate_names[0] == "RunA"
    assert "AnalysisResults" in candidate_names
    assert "slide" in candidate_names


def test_apply_primary_folder_state_enables_output_for_valid_slide(make_napari_viewer, qtbot, tmp_path):
    viewer = make_napari_viewer()
    widget = GeminiQWidget(viewer, gem=None, show_data_widgets=False, show_stitching_widget=True)
    qtbot.addWidget(widget)

    slide_root = Path(tmp_path)
    (slide_root / 'CellStatsDir').mkdir()
    (slide_root / 'RunSummary').mkdir()

    assert widget.browser_output_button.isEnabled() is False
    widget._apply_primary_folder_state(str(slide_root))
    assert widget.browser_output_button.isEnabled() is True


def test_resolve_analysis_dir_accepts_analysisresults_with_fov_dirs(make_napari_viewer, tmp_path):
    viewer = make_napari_viewer()
    widget = GeminiQWidget(viewer, gem=None, show_data_widgets=False, show_stitching_widget=True)

    analysis_results = tmp_path / 'AnalysisResults'
    (analysis_results / 'FOV00001').mkdir(parents=True)
    (analysis_results / 'FOV00002').mkdir(parents=True)

    resolved, reason = widget._resolve_analysis_dir(str(analysis_results))

    assert resolved == str(analysis_results)
    assert 'contains FOV*' in reason


def test_resolve_analysis_dir_accepts_nested_analysisresults_with_fov_dirs(make_napari_viewer, tmp_path):
    viewer = make_napari_viewer()
    widget = GeminiQWidget(viewer, gem=None, show_data_widgets=False, show_stitching_widget=True)

    root = tmp_path / 'RNA3D'
    analysis_results = root / 'AnalysisResults'
    (analysis_results / 'FOV00001').mkdir(parents=True)
    (analysis_results / 'FOV00002').mkdir(parents=True)

    resolved, reason = widget._resolve_analysis_dir(str(root))

    assert resolved == str(analysis_results)
    assert 'contains FOV*' in reason
