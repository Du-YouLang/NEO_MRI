"""Run inside 3D Slicer with --python-script; opens an editable local session."""
import datetime
import hashlib
import json
from pathlib import Path
import sys
import uuid

import numpy as np
import qt
import slicer


def start(session_path):
    previous = getattr(slicer.modules, 'neoLesionSession', None)
    if previous:
        slicer.util.mainWindow().removeDockWidget(previous['panel'])
        previous['panel'].close()
        previous['panel'].deleteLater()
    session = json.loads(Path(session_path).read_text(encoding='utf-8'))
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != session['editor_engine_sha256']:
        raise RuntimeError('Editor bridge changed; run lesion-prepare for the updated code')
    editor_dir = Path(session['editor_dir'])
    editor_dir.mkdir(parents=True, exist_ok=True)
    reference = slicer.util.loadVolume(session['reference'])
    backgrounds = [reference]
    for row in session['volumes'][1:]:
        backgrounds.append(slicer.util.loadVolume(row['path']))
    segmentation = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLSegmentationNode', 'NEO_' + session['lesion_kind'])
    segmentation.CreateDefaultDisplayNodes()
    segmentation.SetReferenceImageGeometryParameterFromVolumeNode(reference)
    segment_id = segmentation.GetSegmentation().AddEmptySegment('lesion', session['lesion_kind'], [1.0, .3, .1])
    initial = session['initial_mask']
    last = editor_dir / 'last_draft.json'
    if last.is_file():
        draft = json.loads(last.read_text(encoding='utf-8'))
        p = Path(draft['mask']).resolve()
        if (draft['prepared_key'] != session['prepared_key'] or not p.is_relative_to(editor_dir.resolve())
                or hashlib.sha256(p.read_bytes()).hexdigest() != draft['mask_sha256']):
            raise RuntimeError('Saved draft identity changed; import or inspect it explicitly')
        initial = str(p)
    mask_node = slicer.util.loadVolume(initial)
    slicer.util.updateSegmentBinaryLabelmapFromArray(slicer.util.arrayFromVolume(mask_node).astype(np.uint8), segmentation, segment_id, reference)
    slicer.mrmlScene.RemoveNode(mask_node)
    slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)
    slicer.util.setSliceViewerLayers(background=reference, foreground=None)
    slicer.util.resetThreeDViews()
    slicer.util.resetSliceViews()
    slicer.util.selectModule('SegmentEditor')
    editor = slicer.modules.segmenteditor.widgetRepresentation().self().editor
    editor.setSegmentationNode(segmentation)
    editor.setSourceVolumeNode(reference)
    editor.setCurrentSegmentID(segment_id)
    panel = qt.QDockWidget('NEO MRI — manual lesion / 手工病灶编辑', slicer.util.mainWindow())
    panel.setMinimumWidth(280)
    panel.setMaximumWidth(400)
    content = qt.QWidget()
    layout = qt.QVBoxLayout(content)
    info = qt.QLabel('Edit the lesion in Segment Editor. Save a draft below.\n保存只生成草稿；返回 neo-mri 导入并审核后才能导出。')
    info.setWordWrap(True)
    layout.addWidget(info)
    choice = qt.QComboBox()
    for row in session['volumes']:
        choice.addItem(row['role'])
    def change_background(index):
        node = backgrounds[index]
        slicer.util.setSliceViewerLayers(background=node, foreground=None)
        editor.setSourceVolumeNode(node)
    choice.connect('currentIndexChanged(int)', change_background)
    layout.addWidget(choice)
    save_button = qt.QPushButton('Save draft / 保存草稿（非审核）')
    layout.addWidget(save_button)
    feedback = qt.QLabel('Drafts are versioned; reopening resumes the last saved draft.')
    feedback.setWordWrap(True)
    layout.addWidget(feedback)
    layout.addStretch(1)
    def save_draft():
        try:
            if segmentation.GetSegmentation().GetNumberOfSegments() != 1:
                raise RuntimeError('This session supports one lesion label; keep stroke and WMH in separate sessions')
            data = slicer.util.arrayFromSegmentBinaryLabelmap(segmentation, segment_id, reference)
            snapshot = editor_dir / (datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S') + '_' + uuid.uuid4().hex[:8])
            snapshot.mkdir()
            label = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLLabelMapVolumeNode')
            try:
                label.CopyOrientation(reference)
                slicer.util.updateVolumeFromArray(label, (data > 0).astype(np.uint8))
                mask_path = snapshot / 'lesion_draft.nii.gz'
                if not slicer.util.saveNode(label, str(mask_path)):
                    raise RuntimeError('Failed to save NIfTI')
                if not slicer.util.saveNode(segmentation, str(snapshot / 'lesion.seg.nrrd')):
                    raise RuntimeError('Failed to save editable segmentation')
            finally:
                slicer.mrmlScene.RemoveNode(label)
            record = dict(prepared_key=session['prepared_key'], mask=str(mask_path),
                          mask_sha256=hashlib.sha256(mask_path.read_bytes()).hexdigest(),
                          editor='3D Slicer', editor_version=slicer.app.applicationVersion,
                          editor_script_sha256=session['editor_engine_sha256'],
                          status='DRAFT_NOT_REVIEWED', saved_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
            temporary = editor_dir / 'last_draft.json.tmp'
            temporary.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding='utf-8')
            temporary.replace(last)
            feedback.setText('Draft saved / 草稿已保存（未审核）\nSnapshot: ' + snapshot.name + '\nReturn to neo-mri lesion-import, then lesion-review.')
        except Exception as exc:
            slicer.util.errorDisplay(str(exc))
    save_button.connect('clicked()', save_draft)
    panel.setWidget(content)
    slicer.util.mainWindow().addDockWidget(qt.Qt.RightDockWidgetArea, panel)
    panel.show()
    # Keep PythonQt wrappers and callbacks alive for the full session.
    slicer.modules.neoLesionSession = dict(panel=panel, save=save_draft, background=change_background,
        reference=reference, segmentation=segmentation, segment_id=segment_id, backgrounds=backgrounds)


if __name__ == '__main__':
    start(sys.argv[-1])
