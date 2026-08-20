"""
Definition of custom QWidget to provide interface for Gemini module.
"""
import napari
from qtpy.QtCore import Qt, QSize, QItemSelectionModel
from qtpy.QtWidgets import (
    QCheckBox,
    QWidget,
    QApplication,
    QGridLayout,
    QPushButton,
    QLabel,
    QComboBox,
    QGroupBox,
    QVBoxLayout,
    QListWidget,
    QListWidgetItem,
    QAbstractItemView,
    QLineEdit,
    QFileDialog,
    QHBoxLayout,
    QLineEdit
)

from qtpy.QtGui import (
    QIcon,
    QPixmap,
    QColor,
    QImage
)

from napari.utils.notifications import (
    notification_manager,
    show_info,
)
from napari.experimental import link_layers
from napari.layers import Labels, Image
from napari.utils.colormaps import AVAILABLE_COLORMAPS, label_colormap
from napari_cosmx._colors import categorical_color_map, colorable_columns
from napari.utils.colormaps.vendored import colors as color_utils
import numpy as np
import pandas as pd
from os import path, listdir, system, sep, name
import os
import glob
import re
from pathlib import Path
from sysconfig import get_path
from ntpath import sep
import subprocess

class GeminiQWidget(QWidget):
    def __init__(self, napari_viewer, gem=None, show_stitching_widget=True, show_data_widgets=True):
        super().__init__()
        self.viewer = napari_viewer
        self.gem = gem
        self.show_stitching_widget = bool(show_stitching_widget)
        self.show_data_widgets = bool(show_data_widgets)

        vbox_layout = QVBoxLayout(self)
        vbox_layout.setContentsMargins(9, 9, 9, 9)

        self.setLayout(vbox_layout)

        if self.show_data_widgets:
            assert self.gem is not None, "gem is required when show_data_widgets=True"
            self.createDatasetSummaryWidget()
            self.createMorphologyImageWidget()
            if self.gem.has_rna:
                self.createTranscriptsWidget()
            if self.gem.has_protein:
                self.createProteinExpressionWidget()
            self.createMetadataWidget()

        if self.show_stitching_widget:
            self.createStitchingWidget()

    def createDatasetSummaryWidget(self):
        groupBox = QGroupBox(self, title="Dataset")
        vbox = QVBoxLayout(groupBox)

        modalities = getattr(self.gem, 'modalities', ())
        modalities_text = ", ".join(modalities) if len(modalities) > 0 else "none detected"
        summary = QLabel(
            f"name: {self.gem.name}\n"
            f"folder: {self.gem.folder}\n"
            f"modalities: {modalities_text}\n"
            f"dimensions: {'3D' if self.gem.is_3d else '2D'}"
        )
        summary.setWordWrap(True)
        vbox.addWidget(summary)
        group_mgr_btn = QPushButton("Layer Groups\u2026")
        group_mgr_btn.setToolTip("Open the Layer Group Manager to organize and toggle layer visibility by group")
        group_mgr_btn.clicked.connect(self._open_group_manager)
        vbox.addWidget(group_mgr_btn)
        groupBox.setLayout(vbox)
        self.layout().addWidget(groupBox)

    def _open_group_manager(self):
        """Open the Layer Group Manager as a separate dock widget."""
        widgets = [w for w in self.viewer.window._dock_widgets.values()
            if isinstance(w.widget(), LayerGroupManagerWidget)]
        if widgets:
            for w in widgets:
                w.show()
        else:
            self.viewer.window.add_dock_widget(
                LayerGroupManagerWidget(self.viewer),
                area='right',
                name='Layer Groups'
            )

    def _z_mapping_summary(self):
        if not getattr(self.gem, 'is_3d', False):
            return ""
        mapping = self.gem.z_plane_mapping() if hasattr(self.gem, 'z_plane_mapping') else {i: z for i, z in enumerate(self.gem.z_slices)}
        pairs = [f"{idx}->{z}" for idx, z in mapping.items()]
        if len(pairs) > 8:
            pairs = pairs[:6] + ['...'] + pairs[-2:]
        return "z index->value: " + ", ".join(pairs)

    def _parse_z_plane_text(self, text, allow_multi=False):
        value = (text or "").strip()
        if value == "":
            return None
        lower = value.lower()
        if lower in {"auto", "current", "any", "all"}:
            if lower in {"any", "all"} and not allow_multi:
                raise ValueError("'all' is only supported for transcript/selection actions.")
            return 'all' if lower == 'any' else lower
        if '-' in value:
            if not allow_multi:
                raise ValueError("z ranges are only supported for transcript/selection actions.")
            a, b = [part.strip() for part in value.split('-', 1)]
            return (int(a), int(b))
        if ',' in value:
            if not allow_multi:
                raise ValueError("multiple z values are only supported for transcript/selection actions.")
            return [int(part.strip()) for part in value.split(',') if part.strip() != ""]
        return int(value)

    def _on_morph_click(self):
        try:
            z_plane = self._parse_z_plane_text(self.morphZPlaneLineEdit.text(), allow_multi=False) if hasattr(self, 'morphZPlaneLineEdit') else None
            self.gem.add_channel(self.channelsComboBox.currentText(),
                self.channelsColormapComboBox.currentText(), z_plane=z_plane)
        except Exception as e:
            show_info(f"Could not add morphology layer: {e}")

    def _on_expr_click(self):
        try:
            z_plane = self._parse_z_plane_text(self.proteinZPlaneLineEdit.text(), allow_multi=False) if hasattr(self, 'proteinZPlaneLineEdit') else None
            self.gem.add_protein(self.proteinComboBox.currentText(), self.colormapComboBox.currentText(), z_plane=z_plane)
        except Exception as e:
            show_info(f"Could not add protein layer: {e}")

    def _on_rna_click(self):
        try:
            z_plane = self._parse_z_plane_text(self.rnaZPlaneLineEdit.text(), allow_multi=True) if hasattr(self, 'rnaZPlaneLineEdit') else None
            point_size = 5
            if hasattr(self, 'rnaPointSizeLineEdit'):
                size_text = self.rnaPointSizeLineEdit.text().strip()
                if size_text != "":
                    point_size = float(size_text)
                if point_size <= 0:
                    raise ValueError("point size must be > 0")
            self.gem.plot_transcripts(gene=self.targetsComboBox.currentText(),
                color=self.colorsComboBox.currentText(), point_size=point_size, z_plane=z_plane)
        except Exception as e:
            show_info(f"Could not add transcript layer: {e}")

    def _get_label_colors(self):
        labels_layer = self._selected_labels_layer()
        if labels_layer is None or labels_layer.features is None:
            return {}
        if self.metaComboBox.currentText() not in labels_layer.features.columns:
            return {}
        if not self.showSelectedCheckbox.isChecked():
            colors = [i.icon().pixmap(QSize(1, 1)).toImage().pixelColor(0,0).name() for i in self.labelListWidget.items()]
            items = [i.text() for i in self.labelListWidget.items()]
            items = pd.Series(items).astype(labels_layer.features[self.metaComboBox.currentText()].dtype)
            return dict(zip(items, colors))
        else:
            colors = [i.icon().pixmap(QSize(1, 1)).toImage().pixelColor(0,0).name() for i in self.labelListWidget.selectedItems()]
            items = [i.text() for i in self.labelListWidget.selectedItems()]
            items = pd.Series(items).astype(labels_layer.features[self.metaComboBox.currentText()].dtype)
            return dict(zip(items, colors))

    def _label_layers(self):
        return [layer for layer in self.viewer.layers if isinstance(layer, Labels)]

    def _refresh_labels_layer_choices(self):
        if not hasattr(self, 'labelsLayerComboBox'):
            return
        current = self.labelsLayerComboBox.currentText()
        self.labelsLayerComboBox.blockSignals(True)
        self.labelsLayerComboBox.clear()
        layer_names = [layer.name for layer in self._label_layers()]
        self.labelsLayerComboBox.addItems(layer_names)
        if current in layer_names:
            self.labelsLayerComboBox.setCurrentText(current)
        elif self.gem.cells_layer is not None and self.gem.cells_layer.name in layer_names:
            self.labelsLayerComboBox.setCurrentText(self.gem.cells_layer.name)
        self.labelsLayerComboBox.blockSignals(False)

    def _selected_labels_layer(self):
        if not hasattr(self, 'labelsLayerComboBox'):
            return self.gem.cells_layer
        layer_name = self.labelsLayerComboBox.currentText().strip()
        if layer_name == "":
            return self.gem.cells_layer
        if layer_name not in self.viewer.layers:
            return None
        layer = self.viewer.layers[layer_name]
        return layer if isinstance(layer, Labels) else None

    def _labels_layer_for_z_value(self, z_value):
        z_value = int(z_value)
        pattern = re.compile(rf"\[z={z_value}\]$")
        for layer in self._label_layers():
            if pattern.search(layer.name):
                return layer
        return None

    def _layer_z_value(self, layer):
        if layer is None:
            return None
        z_value = layer.metadata.get('cosmx_z_plane') if hasattr(layer, 'metadata') else None
        if z_value is not None:
            return int(z_value)
        match = re.search(r"\[z=\s*(-?\d+)\s*\]$", getattr(layer, 'name', ''))
        if match is None:
            return None
        return int(match.group(1))

    def _resolve_coloring_labels_layer(self):
        labels_layer = self._selected_labels_layer()
        if labels_layer is None:
            show_info("Selected labels layer is not available in viewer.")
            return None
        if not self.gem.is_3d or not hasattr(self, 'colorZPlaneLineEdit'):
            return labels_layer

        z_text = self.colorZPlaneLineEdit.text().strip()
        if z_text == "":
            return labels_layer

        try:
            parsed = self._parse_z_plane_text(z_text, allow_multi=False)
        except Exception as e:
            show_info(f"Invalid color z-plane input: {e}")
            return None

        if parsed in {'auto', 'any'}:
            show_info("Color Cells z-plane accepts a single plane or 'current'.")
            return None
        if parsed == 'current':
            current_step = getattr(self.viewer.dims, 'current_step', ())
            z_index = int(current_step[0]) if len(current_step) > 0 else 0
            z_index = max(0, min(z_index, len(self.gem.z_slices) - 1))
            z_value = int(self.gem.z_slices[z_index])
        else:
            z_value = int(parsed)
            if z_value in self.gem.z_slice_to_index:
                z_value = int(z_value)
            elif 0 <= z_value < len(self.gem.z_slices):
                z_value = int(self.gem.z_slices[z_value])
            else:
                show_info(f"Requested z={z_value} is not valid for this dataset.")
                return None

        z_layer = self._labels_layer_for_z_value(z_value)
        if z_layer is None:
            try:
                z_layer = self.gem.add_cell_labels(z_plane=z_value)
                self._refresh_labels_layer_choices()
                if hasattr(self, 'labelsLayerComboBox') and z_layer is not None:
                    self.labelsLayerComboBox.setCurrentText(z_layer.name)
            except Exception as e:
                show_info(
                    f"Could not load labels layer for z={z_value}: {e}. "
                    "Try adding the labels layer first (add_cell_labels)."
                )
                return None
        return z_layer

    def _labels_layer_changed(self, _text):
        labels_layer = self._selected_labels_layer()
        z_value = self._layer_z_value(labels_layer)
        if hasattr(self, 'colorZPlaneLineEdit'):
            self.colorZPlaneLineEdit.blockSignals(True)
            self.colorZPlaneLineEdit.setText("" if z_value is None else str(int(z_value)))
            self.colorZPlaneLineEdit.blockSignals(False)
        self._meta_changed(self.metaComboBox.currentText())

    def _labels_selected(self):
        if self.showSelectedCheckbox.isChecked():
            labels_layer = self._resolve_coloring_labels_layer()
            if labels_layer is None or labels_layer.features is None:
                return
            meta_col = self.metaComboBox.currentText()
            colors = self._get_label_colors()
            if len(colors) == 0:
                return
            selected_items = [i.text() for i in self.labelListWidget.selectedItems()]
            selected_items = pd.Series(selected_items).astype(labels_layer.features[meta_col].dtype)
            cells = labels_layer.features[labels_layer.features[meta_col].isin(selected_items)]['index']
            keep_prev = self.keepColoredLayersCheckbox.isChecked() if hasattr(self, 'keepColoredLayersCheckbox') else False
            self.gem.color_cells(meta_col, colors, subset=cells, labels_layer=labels_layer, keep_previous=keep_prev)
    
    def _show_selected_changed(self, state):
        if self.showSelectedCheckbox.isChecked():
            self._labels_selected()
        else:
            self._meta_changed(self.metaComboBox.currentText())

    def _meta_changed(self, text):
        self.showSelectedCheckbox.setChecked(False)
        labels_layer = self._resolve_coloring_labels_layer()
        if labels_layer is None or labels_layer.features is None:
            self.labelListWidget.setHidden(True)
            self.showSelectedCheckbox.setHidden(True)
            return
        if text == "" or text is None or text not in labels_layer.features.columns:
            self.labelListWidget.setHidden(True)
            self.showSelectedCheckbox.setHidden(True)
            return
        if not self.gem.is_categorical_metadata(text):
            self.labelListWidget.setHidden(True)
            self.showSelectedCheckbox.setHidden(True)
        else:
            self.updateLabelsWidget(text)
            self.labelListWidget.setHidden(False)
            self.showSelectedCheckbox.setHidden(False)
        if text != "" and text is not None:
            keep_prev = self.keepColoredLayersCheckbox.isChecked() if hasattr(self, 'keepColoredLayersCheckbox') else False
            self.gem.color_cells(text, labels_layer=labels_layer, keep_previous=keep_prev)

    def _channel_changed(self, text):
        cmap = self.gem.omero(text, auto=False)['color']
        all_items = [self.channelsColormapComboBox.itemText(i)
            for i in range(self.channelsColormapComboBox.count())]
        if cmap in all_items:
            self.channelsColormapComboBox.setCurrentIndex(all_items.index(cmap))

    def _protein_changed(self, text):
        cmap = self.gem.omero(text, protein=True, auto=False)['color']
        all_items = [self.colormapComboBox.itemText(i)
            for i in range(self.colormapComboBox.count())]
        if cmap in all_items:
            self.colormapComboBox.setCurrentIndex(all_items.index(cmap))

    def _set_workflow_checkboxes(self, read_targets=None, stitch_protein=None, respect_manual=True, reset_manual=False):
        """Programmatically update workflow toggles with optional manual-override protection."""
        if not hasattr(self, '_workflow_manual_override'):
            self._workflow_manual_override = {'rna': False, 'protein': False}
        if reset_manual:
            self._workflow_manual_override = {'rna': False, 'protein': False}

        self._suppress_workflow_checkbox_callbacks = True
        try:
            if read_targets is not None and hasattr(self, 'runReadTargetsCheckbox'):
                if (not respect_manual) or (not self._workflow_manual_override.get('rna', False)):
                    self.runReadTargetsCheckbox.setChecked(bool(read_targets))
            if stitch_protein is not None and hasattr(self, 'stitchProteinCheckbox'):
                if (not respect_manual) or (not self._workflow_manual_override.get('protein', False)):
                    self.stitchProteinCheckbox.setChecked(bool(stitch_protein))
        finally:
            self._suppress_workflow_checkbox_callbacks = False

    def _on_read_targets_checkbox_changed(self, _state):
        if not getattr(self, '_suppress_workflow_checkbox_callbacks', False):
            self._workflow_manual_override['rna'] = True
        self._update_stitch_plan_preview()

    def _on_stitch_protein_checkbox_changed(self, _state):
        if not getattr(self, '_suppress_workflow_checkbox_callbacks', False):
            self._workflow_manual_override['protein'] = True
        self._update_stitch_plan_preview()

    def _check_folder_validity(self) -> bool:
        """Check if the selected primary stitching folder is valid.

        A valid folder must contain CellStatsDir and RunSummary. RNA and protein
        inputs are detected separately to support RNA-only, protein-only, and
        multi-omics workflows.

        Returns:
            bool: True if valid, False if not valid
        """
        isValid = True
        if not path.isdir(os.path.join(self.stitching_folder, 'CellStatsDir')):
            print("No valid CellStatsDir")
            isValid = False
        if not path.isdir(os.path.join(self.stitching_folder, 'RunSummary')):
            print("No valid RunSummary")
            isValid = False 
        return isValid

    def _normalize_selected_dir_path(self, folder_path):
        """Normalize a selected directory path while preserving UNC semantics on Windows."""
        if folder_path is None:
            return ""
        folder_path = str(folder_path).strip()
        if folder_path == "":
            return ""
        return os.path.normpath(folder_path)

    def _apply_primary_folder_state(self, folder_path):
        """Apply primary folder selection and refresh dependent UI state."""
        self.stitching_folder = str(folder_path)
        self.folder_path.setText(folder_path)
        is_valid_folder = self._check_folder_validity()
        if is_valid_folder:
            self._autodetect_stitch_dims()
            has_rna_inputs, rna_reason, _analysis_dir = self._detect_rna_inputs()
            has_protein_inputs, protein_reason = self._detect_protein_inputs()
            self._set_workflow_checkboxes(
                read_targets=has_rna_inputs,
                stitch_protein=has_protein_inputs,
                respect_manual=False,
                reset_manual=True,
            )
            if hasattr(self, 'rnaDetectLabel'):
                self.rnaDetectLabel.setText(f"RNA input detection: {rna_reason}")
            if hasattr(self, 'proteinDetectLabel'):
                self.proteinDetectLabel.setText(f"Protein input detection: {protein_reason}")
            self.selected_folder_label.setText(f"Selected folder: {folder_path}")
            self._update_stitch_plan_preview()
            self.browser_output_button.setEnabled(True)
        else:
            self.selected_folder_label.setText(
                f"Error: the selected folder, {folder_path}, is not a valid slide. Please select another folder."
            )
            self.browser_output_button.setEnabled(False)

    def _candidate_primary_folders(self, folder_path):
        """Return candidate slide roots to probe for CellStatsDir/RunSummary."""
        candidates = []
        current = os.path.normpath(folder_path)
        candidates.append(current)
        parent = os.path.dirname(current)
        if parent and parent != current:
            candidates.append(parent)
        grandparent = os.path.dirname(parent) if parent else None
        if grandparent and grandparent not in candidates and grandparent != parent:
            candidates.append(grandparent)
        return candidates
    
    def _browse_folder(self):
        """ opens finder folder to select, 
            checks validity of folder, 
            and enables stitch button if valid.
        """
        folder_path = QFileDialog.getExistingDirectory(self, "Select Folder")
        if not folder_path:
            return
        folder_path = self._normalize_selected_dir_path(folder_path)
        self.folder_path.setText(folder_path)
        msg = "The selected folder is" + str(folder_path)
        print(msg)
        self._apply_primary_folder_state(str(folder_path))

    def _browse_rna_input_folder(self):
        rna_input_path = QFileDialog.getExistingDirectory(self, "Select RNA Analysis Folder")
        if rna_input_path:
            rna_input_path = self._normalize_selected_dir_path(rna_input_path)
            self.rnaInputPathLineEdit.setText(rna_input_path)
            has_rna_inputs, rna_reason, _analysis_dir = self._detect_rna_inputs()
            self.rnaDetectLabel.setText(f"RNA input detection: {rna_reason}")
            self._set_workflow_checkboxes(read_targets=has_rna_inputs, respect_manual=True)
            # If the user picks a slide root or nearby analysis path first,
            # opportunistically detect and set the primary stitching folder.
            if (not hasattr(self, 'stitching_folder')) or self.stitching_folder in ["", "No folder selected."]:
                for candidate in self._candidate_primary_folders(rna_input_path):
                    prev = getattr(self, 'stitching_folder', "No folder selected.")
                    self.stitching_folder = candidate
                    if self._check_folder_validity():
                        self._apply_primary_folder_state(candidate)
                        break
                    self.stitching_folder = prev
            self._update_stitch_plan_preview()

    def _browse_output_folder(self):
        output_folder_path = QFileDialog.getExistingDirectory(self, "Select Folder")
        print(output_folder_path)
        if output_folder_path:
            output_folder_path = self._normalize_selected_dir_path(output_folder_path)
            self.folder_output_path.setText(output_folder_path)
        print(self.folder_output_path.text())
        self.selected_output_folder_label.setText(f"Result will be sent to: {output_folder_path}")
        self._update_stitch_plan_preview()
        self.stitch_button.setEnabled(True) # enable the stitch button

    def _detect_input_ndim_from_labels(self):
        if not hasattr(self, 'stitching_folder'):
            return None, "No folder selected."
        cell_stats_dir = os.path.join(self.stitching_folder, "CellStatsDir")
        if not path.isdir(cell_stats_dir):
            return None, "CellStatsDir not found."

        pattern_2d = re.compile(r"CELLLABELS_F[0-9]+\.TIF")
        pattern_3d = re.compile(r"CELLLABELS_F[0-9]+_Z[0-9]+\.TIF")
        count_2d = 0

        for root, dirs, files in os.walk(cell_stats_dir):
            for filename in files:
                upper = filename.upper()
                if pattern_3d.match(upper):
                    # For preview/detection, a single 3D-formatted labels file is enough.
                    return 3, f"Detected 3D CellLabels pattern (example: {filename})"
                elif pattern_2d.match(upper):
                    count_2d += 1

        if count_2d > 0:
            return 2, f"Detected {count_2d} CellLabels files with 2D pattern"
        return None, "No CellLabels_F*.tif files found for auto-detection"

    def _autodetect_stitch_dims(self):
        detected_ndim, reason = self._detect_input_ndim_from_labels()
        self.detected_input_ndim = detected_ndim
        if detected_ndim is not None:
            self.inputDimComboBox.setCurrentText(str(detected_ndim))
            self.outputDimComboBox.setCurrentText(str(detected_ndim))
            self.detected_mode_label.setText(
                f"Detected input ndim: {detected_ndim} ({reason}). Default plan set to {detected_ndim}D -> {detected_ndim}D. You can override before stitching."
            )
        else:
            self.detected_mode_label.setText(
                f"Could not auto-detect input ndim ({reason}). Defaulting to current selector values; you can set them manually."
            )
        self._update_stitch_plan_preview()

    def _detect_protein_inputs(self):
        if not hasattr(self, 'stitching_folder'):
            return False, "No folder selected."

        analysis_parent = os.path.join(self.stitching_folder, 'AnalysisResults')
        if path.isdir(analysis_parent):
            for subdir in [i for i in listdir(analysis_parent) if not i.startswith('.')]:
                analysis_dir = os.path.join(analysis_parent, subdir)
                if not path.isdir(analysis_dir):
                    continue
                fov_dirs = [x for x in Path(analysis_dir).glob("FOV*") if x.is_dir()]
                if len(fov_dirs) == 0:
                    continue
                checks = [path.isdir(os.path.join(str(fov_dir), "ProteinImages")) for fov_dir in fov_dirs]
                if all(checks):
                    return True, f"Found ProteinImages for all FOVs under AnalysisResults/{subdir}"
                if any(checks):
                    return False, f"Partial ProteinImages under AnalysisResults/{subdir}; some FOV directories are missing ProteinImages"

        protein_dir = os.path.join(self.stitching_folder, "ProteinDir")
        if path.isdir(protein_dir):
            fov_dirs = [x for x in Path(protein_dir).glob("FOV*") if x.is_dir()]
            if len(fov_dirs) == 0:
                return False, "ProteinDir exists but no FOV* directories were found"
            checks = [path.isdir(os.path.join(str(fov_dir), "ProteinImages")) for fov_dir in fov_dirs]
            if all(checks):
                return True, "Found ProteinImages for all FOVs under ProteinDir"
            if any(checks):
                return False, "Partial ProteinImages under ProteinDir; some FOV directories are missing ProteinImages"
            return False, "ProteinDir found but no FOV*/ProteinImages directories detected"

        return False, "No ProteinDir or valid AnalysisResults/*/FOV*/ProteinImages structure found"

    def _resolve_analysis_dir(self, analysis_root_or_dir):
        if analysis_root_or_dir is None or analysis_root_or_dir == "":
            return None, "No RNA analysis path provided."
        if not path.isdir(analysis_root_or_dir):
            return None, f"Path does not exist: {analysis_root_or_dir}"

        def _has_fov_dirs(folder):
            entries = [
                i for i in listdir(folder)
                if not i.startswith('.') and path.isdir(os.path.join(folder, i))
            ]
            return any(i.upper().startswith('FOV') for i in entries)

        base_name = os.path.basename(os.path.normpath(analysis_root_or_dir))
        if base_name.lower() == 'analysisresults':
            if _has_fov_dirs(analysis_root_or_dir):
                return analysis_root_or_dir, "Using AnalysisResults directly (contains FOV* folders)"
            subdirs = [i for i in listdir(analysis_root_or_dir) if not i.startswith('.') and path.isdir(os.path.join(analysis_root_or_dir, i))]
            if len(subdirs) == 1:
                resolved = os.path.join(analysis_root_or_dir, subdirs[0])
                return resolved, f"Using AnalysisResults child: {subdirs[0]}"
            return None, "AnalysisResults must contain exactly one subdirectory, or provide AnalysisDir directly."

        nested_analysis_results = os.path.join(analysis_root_or_dir, 'AnalysisResults')
        if path.isdir(nested_analysis_results):
            if _has_fov_dirs(nested_analysis_results):
                return nested_analysis_results, "Using nested AnalysisResults directly (contains FOV* folders)"
            subdirs = [i for i in listdir(nested_analysis_results) if not i.startswith('.') and path.isdir(os.path.join(nested_analysis_results, i))]
            if len(subdirs) == 1:
                resolved = os.path.join(nested_analysis_results, subdirs[0])
                return resolved, f"Using nested AnalysisResults child: {subdirs[0]}"
            return None, "Nested AnalysisResults must contain exactly one subdirectory, or provide AnalysisDir directly."

        if _has_fov_dirs(analysis_root_or_dir):
            return analysis_root_or_dir, "Using provided AnalysisDir (contains FOV* folders)"

        # Assume user provided AnalysisDir directly.
        return analysis_root_or_dir, "Using provided AnalysisDir"

    def _detect_rna_inputs(self):
        if not hasattr(self, 'stitching_folder'):
            return False, "No folder selected.", None

        override_path = self.rnaInputPathLineEdit.text().strip() if hasattr(self, 'rnaInputPathLineEdit') else ""
        rna_root = override_path if override_path != "" else self.stitching_folder
        analysis_dir, reason = self._resolve_analysis_dir(rna_root)
        if analysis_dir is None:
            return False, reason, None

        # Confirm RNA data by looking for target call CSV files (same patterns as read-targets CLI).
        # This prevents protein-only folders (which also have AnalysisResults/FOV*) from being
        # incorrectly detected as RNA inputs.
        target_call_pattern = "[a-zA-Z0-9_]*__complete_code_cell_target_call_coord.csv"
        hit = (
            glob.glob(os.path.join(analysis_dir, "FOV[0-9]*", "FOV_Analysis_Summary", target_call_pattern))
            or glob.glob(os.path.join(analysis_dir, "FOV[0-9]*", target_call_pattern))
            or glob.glob(os.path.join(analysis_dir, target_call_pattern))
        )
        if not hit:
            return False, f"AnalysisResults found but no target call CSV files detected (protein-only folder?)", None
        return True, reason, analysis_dir

    def _update_stitch_plan_preview(self):
        input_ndim = int(self.inputDimComboBox.currentText())
        output_ndim = int(self.outputDimComboBox.currentText())
        zslice_text = self.zsliceLineEdit.text().strip()
        z_step_text = self.zStepLineEdit.text().strip()

        if input_ndim == 3 and output_ndim == 2:
            if zslice_text == "":
                zslice_msg = "middle z-plane will be used (if not overridden by --zslice)"
            else:
                zslice_msg = f"zslice={zslice_text}"
        else:
            zslice_msg = "zslice not specified"

        z_step_msg = f"z-step={z_step_text} um" if z_step_text != "" else "z-step auto from image metadata (fallback to default)"
        rna_msg = "read-targets: on" if getattr(self, 'runReadTargetsCheckbox', None) is not None and self.runReadTargetsCheckbox.isChecked() else "read-targets: off"
        protein_msg = "protein stitch: on" if getattr(self, 'stitchProteinCheckbox', None) is not None and self.stitchProteinCheckbox.isChecked() else "protein stitch: off"
        self.planned_mode_label.setText(
            f"Planned stitch mode: {input_ndim}D input -> {output_ndim}D output | {zslice_msg} | {z_step_msg} | {rna_msg} | {protein_msg}"
        )

    def _stitch_images_in_widget(self):
        """ Lower-level function that does stitching.
        """
        def _set_status(text):
            self.stitch_message_label.setText(text)
            # Flush paint/events so users see stage updates while commands run.
            QApplication.processEvents()

        self.stitch_button.setText("Running in the background...")
        msg = "Stitching started.\n"
        _set_status(msg)
        if os.name == 'nt':
            msg += 'Windows detected. Normalizing file paths.\n'
            _set_status(msg)
            self.stitching_folder = self._normalize_selected_dir_path(self.stitching_folder)
        stitch_script_path = os.path.join(get_path("scripts"), "stitch-images")
        if os.name == 'nt':
            stitch_script_path += '.exe'
            stitch_script_path = self._normalize_selected_dir_path(stitch_script_path)
        print(stitch_script_path)
        read_targets_script_path = os.path.join(get_path("scripts"), "read-targets")
        if os.name == 'nt':
            read_targets_script_path += '.exe'
            read_targets_script_path = self._normalize_selected_dir_path(read_targets_script_path)
        stitch_expression_script_path = os.path.join(get_path("scripts"), "stitch-expression")
        if os.name == 'nt':
            stitch_expression_script_path += '.exe'
            stitch_expression_script_path = self._normalize_selected_dir_path(stitch_expression_script_path)
        if not path.isfile(stitch_script_path):
            _set_status(msg + "Could not find stitch-images in path " + stitch_script_path)
            print("could not fine stitch script in path: " + stitch_script_path)
            return None
        else:
            msg += "Found console script stitch-images.\n"
            _set_status(msg)

        has_rna_inputs, rna_reason, analysis_dir = self._detect_rna_inputs()
        has_protein_inputs, protein_reason = self._detect_protein_inputs()
        if not has_rna_inputs and not has_protein_inputs:
            msg += "No RNA or protein inputs detected. Need AnalysisResults (or external RNA path) and/or ProteinDir / AnalysisResults/*/ProteinImages.\n"
            _set_status(msg)
            self.stitch_button.setText("Stitch")
            return None
        run_read_targets = getattr(self, 'runReadTargetsCheckbox', None) is not None and self.runReadTargetsCheckbox.isChecked()
        if run_read_targets and not has_rna_inputs:
            msg += f"RNA inputs not detected, skipping read-targets: {rna_reason}.\n"
            run_read_targets = False
            _set_status(msg)
        if run_read_targets and not path.isfile(read_targets_script_path):
            msg += "Could not find read-targets script; RNA target stitching will be skipped.\n"
            run_read_targets = False
            _set_status(msg)
        elif run_read_targets:
            msg += "Found console script read-targets.\n"
            _set_status(msg)
        run_stitch_expression = getattr(self, 'stitchProteinCheckbox', None) is not None and self.stitchProteinCheckbox.isChecked()
        if run_stitch_expression and not has_protein_inputs:
            msg += f"Protein inputs not detected, skipping stitch-expression: {protein_reason}.\n"
            run_stitch_expression = False
            _set_status(msg)
        if run_stitch_expression and not path.isfile(stitch_expression_script_path):
            msg += "Could not find stitch-expression script; protein stitching will be skipped.\n"
            run_stitch_expression = False
            _set_status(msg)
        elif run_stitch_expression:
            msg += "Found console script stitch-expression.\n"
            _set_status(msg)
        CellStatsDir = os.path.join(self.stitching_folder, "CellStatsDir")
        msg += "CellStatsDir is " + CellStatsDir + ".\n"
        _set_status(msg)
        RunSummaryDir = os.path.join(self.stitching_folder, 'RunSummary')
        msg += "RunSummaryDir is " + RunSummaryDir + ".\n"
        _set_status(msg)
        if run_read_targets and analysis_dir is not None:
            msg += "AnalysisDir is " + analysis_dir + ".\n"
        msg += "Running stitch-images...\n"
        _set_status(msg)
        cmd = [
            stitch_script_path, "-i", CellStatsDir, "-f", RunSummaryDir, "-o", self.folder_output_path.text()
        ]
        input_ndim = int(self.inputDimComboBox.currentText())
        output_ndim = int(self.outputDimComboBox.currentText())
        if input_ndim == 2 and output_ndim == 3:
            msg += "\nInvalid mode selection: cannot stitch 2D input into 3D output."
            _set_status(msg)
            self.stitch_button.setText("Stitch")
            return None
        self._update_stitch_plan_preview()
        cmd.extend(["--input-ndim", str(input_ndim), "--output-ndim", str(output_ndim)])

        zslice_text = self.zsliceLineEdit.text().strip()
        if zslice_text != "":
            cmd.extend(["--zslice", zslice_text])

        z_step_text = self.zStepLineEdit.text().strip()
        if z_step_text != "":
            cmd.extend(["--z-step-um", z_step_text])

        print(cmd)
        self._run_command(cmd)
        msg += "\nFinished stitching images.\n"
        _set_status(msg)

        if run_read_targets and analysis_dir is not None:
            msg += "Running read-targets...\n"
            _set_status(msg)
            cmd2 = [
                read_targets_script_path, analysis_dir, "-o", self.folder_output_path.text()
            ]
            print(cmd2)
            self._run_command(cmd2)
            msg += "Finished read-targets.\n"
            _set_status(msg)

        if run_stitch_expression:
            has_protein_inputs, protein_reason = self._detect_protein_inputs()
            if has_protein_inputs:
                msg += "\nRunning stitch-expression...\n"
                _set_status(msg)
                cmd3 = [
                    stitch_expression_script_path,
                    "-i", self.stitching_folder,
                    "-o", self.folder_output_path.text(),
                    "--ndim", str(output_ndim),
                ]
                if z_step_text != "":
                    cmd3.extend(["--z-step-um", z_step_text])
                print(cmd3)
                self._run_command(cmd3)
                msg += "Finished stitch-expression.\n"
                _set_status(msg)
            else:
                msg += f"\nSkipping stitch-expression: {protein_reason}.\n"
                _set_status(msg)

        msg += "\nStitching workflow finished.\nSee output folder for results.\n"
        _set_status(msg)
        self.stitch_button.setText("Reselect input/output to stitch again")
        self.stitch_button.setEnabled(False)

    def _run_command(self, command):
        subprocess.run(command)

    def update_metadata(self, path):
        # TODO: would be nice to merge metadata, but replace for now
        self.gem.read_metadata(path)
        if self.gem.metadata is not None:
            self.metaComboBox.clear()
            self.metaComboBox.addItems(colorable_columns(self.gem.metadata.columns))

    def createMorphologyImageWidget(self):
        groupBox = QGroupBox(self, title="Morphology Images")

        btn = QPushButton("Add layer")
        btn.clicked.connect(self._on_morph_click)

        boxp = QComboBox(groupBox)
        boxp.addItems(self.gem.channels)

        colormaps = [i for i in list(AVAILABLE_COLORMAPS.keys()) if i not in ['label_colormap', 'custom']]
        boxc = QComboBox(groupBox)
        boxc.addItems(colormaps)

        self.channelsComboBox = boxp
        self.channelsColormapComboBox = boxc

        vbox = QVBoxLayout(groupBox)
        grid = QGridLayout()

        grid.addWidget(QLabel('channel:'), 1, 0)
        grid.addWidget(self.channelsComboBox, 1, 1)
        grid.addWidget(QLabel('colormap:'), 2, 0)
        grid.addWidget(self.channelsColormapComboBox, 2, 1)
        if self.gem.is_3d:
            self.morphZPlaneLineEdit = QLineEdit(groupBox)
            self.morphZPlaneLineEdit.setPlaceholderText("Optional: current or z value (e.g. 3)")
            grid.addWidget(QLabel('z plane:'), 3, 0)
            grid.addWidget(self.morphZPlaneLineEdit, 3, 1)
            mapping_label = QLabel(self._z_mapping_summary())
            mapping_label.setWordWrap(True)
            vbox.addWidget(mapping_label)
        vbox.addLayout(grid)
        vbox.addWidget(btn)
        groupBox.setLayout(vbox)

        self._channel_changed(self.channelsComboBox.currentText())
        self.channelsComboBox.currentTextChanged.connect(self._channel_changed)

        self.layout().addWidget(groupBox)

    def createProteinExpressionWidget(self):
        groupBox = QGroupBox(self, title="Protein Expression")

        btn = QPushButton("Add layer")
        btn.clicked.connect(self._on_expr_click)

        boxp = QComboBox(groupBox)
        boxp.addItems(self.gem.proteins)

        colormaps = [i for i in list(AVAILABLE_COLORMAPS.keys()) if i not in ['label_colormap', 'custom']]
        boxc = QComboBox(groupBox)
        boxc.addItems(colormaps)

        self.proteinComboBox = boxp
        self.colormapComboBox = boxc

        vbox = QVBoxLayout(groupBox)
        grid = QGridLayout()

        grid.addWidget(QLabel('protein:'), 1, 0)
        grid.addWidget(self.proteinComboBox, 1, 1)
        grid.addWidget(QLabel('colormap:'), 2, 0)
        grid.addWidget(self.colormapComboBox, 2, 1)
        if self.gem.is_3d:
            self.proteinZPlaneLineEdit = QLineEdit(groupBox)
            self.proteinZPlaneLineEdit.setPlaceholderText("Optional: current or z value (e.g. 3)")
            grid.addWidget(QLabel('z plane:'), 3, 0)
            grid.addWidget(self.proteinZPlaneLineEdit, 3, 1)
        vbox.addLayout(grid)
        vbox.addWidget(btn)
        groupBox.setLayout(vbox)

        self._protein_changed(self.proteinComboBox.currentText())
        self.proteinComboBox.currentTextChanged.connect(self._protein_changed)

        self.layout().addWidget(groupBox)

    def createTranscriptsWidget(self):
        groupBox = QGroupBox(self, title="RNA Transcripts")

        btn = QPushButton("Add layer")
        btn.clicked.connect(self._on_rna_click)

        boxp = QComboBox(groupBox)
        boxp.addItems(self.gem.genes)

        colors = ['white', 'red', 'green', 'blue',
            'magenta', 'yellow', 'cyan']
        boxc = QComboBox(groupBox)
        boxc.addItems(colors)

        self.targetsComboBox = boxp
        self.colorsComboBox = boxc

        vbox = QVBoxLayout(groupBox)
        grid = QGridLayout()

        grid.addWidget(QLabel('target:'), 1, 0)
        grid.addWidget(self.targetsComboBox, 1, 1)
        grid.addWidget(QLabel('color:'), 2, 0)
        grid.addWidget(self.colorsComboBox, 2, 1)
        self.rnaPointSizeLineEdit = QLineEdit(groupBox)
        self.rnaPointSizeLineEdit.setText("5")
        self.rnaPointSizeLineEdit.setPlaceholderText("Point size (e.g. 5)")
        grid.addWidget(QLabel('point size:'), 3, 0)
        grid.addWidget(self.rnaPointSizeLineEdit, 3, 1)
        if self.gem.is_3d:
            self.rnaZPlaneLineEdit = QLineEdit(groupBox)
            self.rnaZPlaneLineEdit.setPlaceholderText("Optional: auto/current/all/3/2-4/2,4")
            grid.addWidget(QLabel('z plane(s):'), 4, 0)
            grid.addWidget(self.rnaZPlaneLineEdit, 4, 1)
        vbox.addLayout(grid)
        vbox.addWidget(btn)
        groupBox.setLayout(vbox)

        self.layout().addWidget(groupBox)

    def createMetadataWidget(self):
        groupBox = QGroupBox(self, title="Color Cells")

        boxc = QComboBox(groupBox)
        boxc.toolTip = "Open a _metadata.csv file to populate dropdown"
        if self.gem.metadata is not None:
            boxc.addItems(colorable_columns(self.gem.metadata.columns))
        else:
            boxc.addItems([])

        self.metaComboBox = boxc

        vbox = QVBoxLayout(groupBox)
        grid = QGridLayout()

        self.labelsLayerComboBox = QComboBox(groupBox)
        self._refresh_labels_layer_choices()
        self.labelsLayerComboBox.currentTextChanged.connect(self._labels_layer_changed)

        self.colorZPlaneLineEdit = QLineEdit(groupBox)
        self.colorZPlaneLineEdit.setPlaceholderText("Optional (3D): current or z value/index")
        self.colorZPlaneLineEdit.textChanged.connect(lambda _t: self._meta_changed(self.metaComboBox.currentText()))

        grid.addWidget(QLabel('labels layer:'), 0, 0)
        grid.addWidget(self.labelsLayerComboBox, 0, 1)
        grid.addWidget(QLabel('color z plane:'), 1, 0)
        grid.addWidget(self.colorZPlaneLineEdit, 1, 1)
        grid.addWidget(QLabel('column:'), 2, 0)
        grid.addWidget(self.metaComboBox, 2, 1)
        vbox.addLayout(grid)
        if self.gem.is_3d:
            vbox.addWidget(QLabel(self._z_mapping_summary()))
        vbox.addWidget(QLabel("Color Cells uses the selected labels layer as source."))

        listl = QListWidget(groupBox)
        listl.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.labelListWidget = listl
        self.showSelectedCheckbox = QCheckBox('Only show selected labels')
        self.showSelectedCheckbox.stateChanged.connect(self._show_selected_changed)
        
        self.keepColoredLayersCheckbox = QCheckBox('Keep previous colored layers')
        self.keepColoredLayersCheckbox.setToolTip('When checked, new colored metadata layers are added instead of replacing the previous one')

        self.labelListWidget.setHidden(True)
        self.showSelectedCheckbox.setHidden(True)
        self._meta_changed(self.metaComboBox.currentText())
        if self.gem.cells_layer is not None:
            self.gem.cells_layer.visible = False
        self.viewer.layers.events.inserted.connect(lambda _e: self._refresh_labels_layer_choices())
        self.viewer.layers.events.removed.connect(lambda _e: self._refresh_labels_layer_choices())
        self.metaComboBox.currentTextChanged.connect(self._meta_changed)

        listl.itemSelectionChanged.connect(self._labels_selected)
        vbox.addWidget(self.labelListWidget)
        vbox.addWidget(self.showSelectedCheckbox)
        vbox.addWidget(self.keepColoredLayersCheckbox)
        groupBox.setLayout(vbox)
        self.layout().addWidget(groupBox)

    def createStitchingWidget(self):
        # set initial value
        self.stitching_folder = "No folder selected." 
        self.folder_output_path = "No output folder selected."
        groupBox = QGroupBox(self, title="Stitch Images")

        # Widget description text
        folder_info_label = QLabel(
            "Step 1: Select PRIMARY stitching folder (must contain CellStatsDir and RunSummary)."
        )
        folder_info_label.setWordWrap(True)
        primary_folder_hint = QLabel(
            "The widget auto-detects RNA/protein inputs from the primary folder and auto-selects "
            "read-targets / stitch-expression checkboxes. For split multi-omics inputs, choose the protein run folder as primary."
        )
        primary_folder_hint.setWordWrap(True)
        additional_rna_hint = QLabel(
            "Step 2 (optional, multi-omics): Select ADDITIONAL RNA folder for RNA input path."
        )
        additional_rna_hint.setWordWrap(True)
        self.detected_input_ndim = None

        # Ingestion mode controls
        self.inputDimComboBox = QComboBox(groupBox)
        self.inputDimComboBox.addItems(["2", "3"])
        self.inputDimComboBox.setCurrentText("2")

        self.outputDimComboBox = QComboBox(groupBox)
        self.outputDimComboBox.addItems(["2", "3"])
        self.outputDimComboBox.setCurrentText("2")

        self.zsliceLineEdit = QLineEdit(groupBox)
        self.zsliceLineEdit.setPlaceholderText("Optional (3D input -> 2D output)")

        self.zStepLineEdit = QLineEdit(groupBox)
        self.zStepLineEdit.setPlaceholderText("Optional, e.g. 0.8")

        self.detected_mode_label = QLabel("Auto-detection runs after you select an input folder.")
        self.detected_mode_label.setWordWrap(True)
        self.planned_mode_label = QLabel("Planned stitch mode: 2D input -> 2D output")
        self.planned_mode_label.setWordWrap(True)
        self._workflow_manual_override = {'rna': False, 'protein': False}
        self._suppress_workflow_checkbox_callbacks = False
        self.stitchProteinCheckbox = QCheckBox("Run stitch-expression for protein (auto-detected)")
        self.stitchProteinCheckbox.setChecked(False)
        self.stitchProteinCheckbox.stateChanged.connect(self._on_stitch_protein_checkbox_changed)
        self.runReadTargetsCheckbox = QCheckBox("Run read-targets for RNA (auto-detected)")
        self.runReadTargetsCheckbox.setChecked(False)
        self.runReadTargetsCheckbox.stateChanged.connect(self._on_read_targets_checkbox_changed)
        self.rnaInputPathLineEdit = QLineEdit(groupBox)
        self.rnaInputPathLineEdit.setPlaceholderText("Optional override: RNA AnalysisDir / AnalysisResults / parent folder")
        self.rnaInputPathLineEdit.textChanged.connect(lambda _t: self._update_stitch_plan_preview())
        self.rnaBrowseButton = QPushButton("Browse Additional RNA Folder (Optional)")
        self.rnaBrowseButton.clicked.connect(self._browse_rna_input_folder)
        self.rnaDetectLabel = QLabel("RNA input detection runs after selecting the primary folder.")
        self.rnaDetectLabel.setWordWrap(True)
        self.proteinDetectLabel = QLabel("Protein input detection runs after selecting the primary folder.")
        self.proteinDetectLabel.setWordWrap(True)

        self.inputDimComboBox.currentTextChanged.connect(lambda _t: self._update_stitch_plan_preview())
        self.outputDimComboBox.currentTextChanged.connect(lambda _t: self._update_stitch_plan_preview())
        self.zsliceLineEdit.textChanged.connect(lambda _t: self._update_stitch_plan_preview())
        self.zStepLineEdit.textChanged.connect(lambda _t: self._update_stitch_plan_preview())

        # Button for browser
        self.folder_path = QLineEdit()
        self.folder_path.setEnabled(False)

        self.browse_button = QPushButton("Browse Primary Folder")
        self.browse_button.clicked.connect(self._browse_folder)

        # Visual feedback for user (selected input path)
        self.selected_folder_label = QLabel(self.stitching_folder)
        self.selected_folder_label.setWordWrap(True)

        # Button for output folder
        self.folder_output_path = QLineEdit()
        self.folder_output_path.setEnabled(False)

        self.browser_output_button = QPushButton("Choose Output Folder")
        self.browser_output_button.clicked.connect(self._browse_output_folder)
        self.browser_output_button.setEnabled(False)

        # Visual feedback for user (selected output path)
        self.selected_output_folder_label = QLabel(self.folder_output_path)
        self.selected_output_folder_label.setWordWrap(True)

        # Button for calling stitching function
        self.stitch_button = QPushButton("Stitch")
        self.stitch_button.clicked.connect(self._stitch_images_in_widget)
        self.stitch_button.setEnabled(False)

        # Visual feedback for user (stitch path)
        self.stitch_message_label = QLabel("")
        self.stitch_message_label.setWordWrap(True)

        # Configure panels
        vbox = QVBoxLayout(groupBox)
        grid = QGridLayout()
        grid.addWidget(folder_info_label, 1, 0)
        grid.addWidget(QLabel("Input ndim:"), 2, 0)
        grid.addWidget(self.inputDimComboBox, 2, 1)
        grid.addWidget(QLabel("Output ndim:"), 3, 0)
        grid.addWidget(self.outputDimComboBox, 3, 1)
        grid.addWidget(QLabel("Z slice:"), 4, 0)
        grid.addWidget(self.zsliceLineEdit, 4, 1)
        grid.addWidget(QLabel("Z step (um):"), 5, 0)
        grid.addWidget(self.zStepLineEdit, 5, 1)
        vbox.addLayout(grid)
        vbox.addWidget(primary_folder_hint)
        vbox.addWidget(self.browse_button)
        vbox.addWidget(self.selected_folder_label)
        vbox.addWidget(self.detected_mode_label)
        vbox.addWidget(self.planned_mode_label)
        vbox.addWidget(self.runReadTargetsCheckbox)
        vbox.addWidget(self.stitchProteinCheckbox)
        vbox.addWidget(self.rnaDetectLabel)
        vbox.addWidget(self.proteinDetectLabel)
        vbox.addWidget(additional_rna_hint)
        vbox.addWidget(self.rnaInputPathLineEdit)
        vbox.addWidget(self.rnaBrowseButton)
        vbox.addWidget(self.browser_output_button)
        vbox.addWidget(self.selected_output_folder_label)
        vbox.addWidget(self.stitch_button)
        vbox.addWidget(self.stitch_message_label)
        groupBox.setLayout(vbox)    

        self._update_stitch_plan_preview()

        # render
        self.layout().addWidget(groupBox)

    def updateLabelsWidget(self, meta_col):
        self.labelListWidget.clear()
        vals = sorted(np.unique(self.gem.metadata[meta_col]))
        if self.gem.adata is not None and meta_col + "_colors" in self.gem.adata.uns:
            # get colors from AnnData object
            color = dict(zip(self.gem.adata.obs[meta_col].cat.categories, self.gem.adata.uns[meta_col + "_colors"]))
            cols = np.vstack((np.zeros(4, dtype='float64'),
                color_utils.to_rgba_array([color[i] for i in vals])))
        else:
            # Prefer a per-column '<col>_color' then legacy 'hex_color' for
            # consistent legend colors across slides; otherwise auto-generate.
            hex_map = categorical_color_map(self.gem.metadata, meta_col)
            if hex_map is not None:
                cols = np.vstack((np.zeros(4, dtype='float64'),
                    color_utils.to_rgba_array([hex_map.get(v, '#808080') for v in vals])))
            else:
                cols = label_colormap(len(vals)+1).colors
        for i,n in enumerate(vals):
            pmap = QPixmap(24, 24)
            rgba = cols[i+1]
            color = np.round(255 * rgba).astype(int)
            pmap.fill(QColor(*list(color)))
            icon = QIcon(pmap)
            qitem = QListWidgetItem(icon, str(n), self.labelListWidget)


def make_stitching_widget(viewer=None):
    """Factory for launching Stitch Images as a standalone dock widget.

    Some napari versions/plugin entry paths may not inject a viewer argument
    into widget factory functions. In that case, fall back to the currently
    active viewer.
    """
    if viewer is None:
        viewer = napari.current_viewer()
        if viewer is None:
            raise RuntimeError("No active napari viewer found.")
    return GeminiQWidget(
        viewer,
        gem=None,
        show_stitching_widget=True,
        show_data_widgets=False,
    )


class LayerGroupManagerWidget(QWidget):
    """Manage named visibility groups of napari layers.

    Create groups, assign any combination of loaded layers to them, then
    show or hide the entire group with one click.  Groups are stored in
    memory for the session; they are not persisted to disk.
    """

    def __init__(self, napari_viewer):
        super().__init__()
        self.viewer = napari_viewer
        self._groups = {}  # {group_name: [layer_name, ...]}
        self._build_ui()
        events = self.viewer.layers.events
        events.removed.connect(self._on_layer_removed)
        # Compatibility across napari versions:
        # - newer versions expose `layers.events.renamed`
        # - older versions may only expose `layers.events.changed`
        renamed_event = getattr(events, 'renamed', None)
        if renamed_event is not None:
            renamed_event.connect(self._on_layer_renamed)
        else:
            changed_event = getattr(events, 'changed', None)
            if changed_event is not None:
                changed_event.connect(self._on_layer_renamed)

    def _build_ui(self):
        vbox = QVBoxLayout(self)
        vbox.setContentsMargins(9, 9, 9, 9)

        # New group row
        new_row = QHBoxLayout()
        self.newGroupLineEdit = QLineEdit()
        self.newGroupLineEdit.setPlaceholderText("New group name\u2026")
        new_btn = QPushButton("New Group")
        new_btn.clicked.connect(self._create_group)
        new_row.addWidget(self.newGroupLineEdit)
        new_row.addWidget(new_btn)
        vbox.addLayout(new_row)

        # Group list
        vbox.addWidget(QLabel("Groups:"))
        self.groupListWidget = QListWidget()
        self.groupListWidget.setMaximumHeight(130)
        self.groupListWidget.currentRowChanged.connect(self._group_selection_changed)
        vbox.addWidget(self.groupListWidget)

        # Per-group visibility / delete controls
        ctrl_row = QHBoxLayout()
        show_btn = QPushButton("Show")
        show_btn.setToolTip("Make all layers in this group visible")
        show_btn.clicked.connect(lambda: self._set_group_visibility(True))
        hide_btn = QPushButton("Hide")
        hide_btn.setToolTip("Make all layers in this group invisible")
        hide_btn.clicked.connect(lambda: self._set_group_visibility(False))
        del_btn = QPushButton("Delete Group")
        del_btn.clicked.connect(self._delete_group)
        ctrl_row.addWidget(show_btn)
        ctrl_row.addWidget(hide_btn)
        ctrl_row.addWidget(del_btn)
        vbox.addLayout(ctrl_row)

        # Layers in the currently selected group
        vbox.addWidget(QLabel("Layers in group:"))
        self.groupLayersWidget = QListWidget()
        self.groupLayersWidget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        vbox.addWidget(self.groupLayersWidget)

        # Assignment controls
        assign_row = QHBoxLayout()
        add_btn = QPushButton("Add selected viewer layers")
        add_btn.setToolTip("Assign the layers currently selected in the napari viewer to this group")
        add_btn.clicked.connect(self._add_selected_layers)
        remove_btn = QPushButton("Remove")
        remove_btn.setToolTip("Remove highlighted items from this group")
        remove_btn.clicked.connect(self._remove_selected_from_group)
        assign_row.addWidget(add_btn)
        assign_row.addWidget(remove_btn)
        vbox.addLayout(assign_row)

        self.setLayout(vbox)

    def _current_group_name(self):
        item = self.groupListWidget.currentItem()
        return item.text() if item is not None else None

    def _create_group(self):
        name = self.newGroupLineEdit.text().strip()
        if name == "":
            show_info("Enter a group name first.")
            return
        if name in self._groups:
            show_info(f"Group '{name}' already exists.")
            return
        self._groups[name] = []
        self.groupListWidget.addItem(name)
        self.newGroupLineEdit.clear()
        self.groupListWidget.setCurrentRow(self.groupListWidget.count() - 1)

    def _delete_group(self):
        name = self._current_group_name()
        if name is None:
            return
        del self._groups[name]
        self.groupListWidget.takeItem(self.groupListWidget.currentRow())
        self.groupLayersWidget.clear()

    def _group_selection_changed(self, _row):
        self._refresh_group_layers()

    def _refresh_group_layers(self):
        self.groupLayersWidget.clear()
        name = self._current_group_name()
        if name is None:
            return
        for layer_name in self._groups[name]:
            self.groupLayersWidget.addItem(layer_name)

    def _add_selected_layers(self):
        group_name = self._current_group_name()
        if group_name is None:
            show_info("Create or select a group first.")
            return
        selected = list(self.viewer.layers.selection)
        if not selected:
            show_info("Select one or more layers in the napari viewer first.")
            return
        existing = set(self._groups[group_name])
        added = 0
        for layer in selected:
            if layer.name not in existing:
                self._groups[group_name].append(layer.name)
                existing.add(layer.name)
                added += 1
        self._refresh_group_layers()
        if added == 0:
            show_info("All selected layers are already in this group.")

    def _remove_selected_from_group(self):
        group_name = self._current_group_name()
        if group_name is None:
            return
        remove_names = {item.text() for item in self.groupLayersWidget.selectedItems()}
        self._groups[group_name] = [n for n in self._groups[group_name] if n not in remove_names]
        self._refresh_group_layers()

    def _set_group_visibility(self, visible):
        group_name = self._current_group_name()
        if group_name is None:
            show_info("Select a group first.")
            return
        layer_names = set(self._groups[group_name])
        for layer in self.viewer.layers:
            if layer.name in layer_names:
                layer.visible = visible

    def _on_layer_removed(self, event):
        removed_name = event.value.name
        for group in self._groups.values():
            if removed_name in group:
                group.remove(removed_name)
        self._refresh_group_layers()

    def _on_layer_renamed(self, event):
        # Keep groups consistent after a rename: prune entries that no longer
        # match any viewer layer (old name disappears, new name will be added
        # manually by the user if desired).
        viewer_names = {layer.name for layer in self.viewer.layers}
        for group_name in self._groups:
            self._groups[group_name] = [n for n in self._groups[group_name] if n in viewer_names]
        self._refresh_group_layers()