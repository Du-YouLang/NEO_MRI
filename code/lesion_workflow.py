"""Model-free manual lesion workflow. No registration, inference, or GLM calls."""
from pathlib import Path
import json
import re
import shutil
import html
import gzip

import nibabel as nib
import numpy as np

from lesion_state import StageError, digest, lock, now, read, sha, write

CODE = Path(__file__).resolve().parent
# A wheel installs flat modules into site-packages, not a source code/ folder.
# Protect the module directory there without treating its parent as a checkout.
REPOSITORY = CODE.parent if CODE.name == 'code' and (CODE.parent / 'pyproject.toml').is_file() else CODE

# This exact predecessor differed only in validation relevant to the present
# update. Its byte-verified preparations remain usable AFTER current validation;
# newly created preparations and imports record the current engine hash.
COMPATIBLE_PREPARATION_ENGINES = (
    'e0e0e874eefa02d5220734387876a4604bb0fcade013b8483078eccdbfec75d8',
    '5bee7a85abb6d39899d23531110ca81ac0cc365f5a9574cd4c011449515ac446',
    '5bee7a85abb6d39899d23531110ca81ac0cc365f5a9574cd4c011449515ac446',
)


def same_grid(a, b):
    return (a.header.get_xyzt_units()[0] == b.header.get_xyzt_units()[0]
            and a.shape == b.shape and np.allclose(a.affine, b.affine, rtol=0, atol=1e-5))


def image(path):
    try:
        im = nib.load(path)
    except nib.filebasedimages.ImageFileError as exc:
        raise StageError('Cannot read NIfTI image: ' + str(path)) from exc
    if len(im.shape) != 3 or any(n < 1 for n in im.shape):
        raise StageError('A nonempty 3D image is required: ' + str(path))
    if im.header.get_xyzt_units()[0] != 'mm':
        raise StageError('NIfTI spatial units must explicitly be mm; convert coordinates/geometry correctly before importing: ' + str(path))
    if not np.isfinite(im.affine).all() or abs(np.linalg.det(im.affine[:3, :3])) < 1e-10:
        raise StageError('Invalid image affine: ' + str(path))
    q, qc = im.get_qform(coded=True)
    s, sc = im.get_sform(coded=True)
    if not qc and not sc:
        raise StageError('Image has no declared qform/sform spatial mapping: ' + str(path))
    if qc and sc and not np.allclose(q, s, rtol=0, atol=1e-3):
        raise StageError('Conflicting qform/sform; review image geometry: ' + str(path))
    if not np.isfinite(im.get_fdata(dtype=np.float32)).all():
        raise StageError('Image contains nonfinite values: ' + str(path))
    return im


def load_config(path, root):
    root = Path(root).resolve()
    def resolve(value):
        p = Path(value)
        return p.resolve() if p.is_absolute() else (root / p).resolve()
    path = resolve(path)
    c = read(path)
    if not isinstance(c, dict):
        raise StageError('Lesion configuration must be a JSON object')
    allowed = {'schema', 'template_only', 'case_id', 'reference_t1', 'review_images',
               'initial_mask', 'initial_label', 'workspace', 'editor', 'lesion_kind'}
    if set(c) - allowed or c.get('schema') != 'neo-lesion-1':
        raise StageError('Unknown lesion schema or configuration keys')
    if c.get('template_only'):
        raise StageError('Copy and fill the lesion template; remove template_only')
    if not re.fullmatch(r'[A-Za-z0-9_-]+', c['case_id']) or not re.search(r'[A-Za-z0-9]', c['case_id']):
        raise StageError('Invalid case_id')
    if c.get('lesion_kind') not in ['stroke', 'wmh', 'other']:
        raise StageError('Declare lesion_kind: stroke, wmh, or other')
    c['workspace'] = str(resolve(c['workspace']))
    ws = Path(c['workspace'])
    # Scientific outputs must not be placed in this source checkout.
    if ws == REPOSITORY or ws.is_relative_to(REPOSITORY):
        raise StageError('Lesion workspace must be outside the code repository')
    c['reference_t1'] = str(resolve(c['reference_t1']))
    for row in c.get('review_images', []):
        if set(row) != {'modality', 'image', 'space'} or row['modality'] not in ['T2w', 'FLAIR'] or row['space'] != 'selected-T1w':
            raise StageError('Review images require modality T2w/FLAIR and space selected-T1w')
        row['image'] = str(resolve(row['image']))
    if c.get('initial_mask'):
        c['initial_mask'] = str(resolve(c['initial_mask']))
    editor = c.get('editor', {})
    if not isinstance(editor, dict) or set(editor) - {'type', 'executable'} or editor.get('type', 'slicer') != 'slicer':
        raise StageError('Editor must be a Slicer configuration object')
    if editor.get('executable'):
        editor['executable'] = str(resolve(editor['executable']))
    c['editor'] = editor
    sources = [c['reference_t1']] + [x['image'] for x in c.get('review_images', [])]
    if c.get('initial_mask'):
        sources.append(c['initial_mask'])
    if c.get('initial_label') is not None and not c.get('initial_mask'):
        raise StageError('initial_label requires initial_mask')
    for source in sources:
        p = Path(source)
        if ws == p or ws in p.parents or p in ws.parents:
            raise StageError('Workspace and source inputs must be separate')
    return c


def manifest(folder):
    return [dict(path=p.relative_to(folder).as_posix(), sha256=sha(p))
            for p in sorted(folder.rglob('*')) if p.is_file()]


def verify(folder, rows):
    for row in rows:
        p = (folder / row['path']).resolve()
        if not p.is_relative_to(folder.resolve()) or not p.is_file() or sha(p) != row['sha256']:
            raise StageError('Registered lesion file missing or changed: ' + str(p))


def copy_verified(source, target, expected):
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
    if sha(target) != expected:
        raise StageError('Input changed while copying: ' + str(source))


def save_mask(data, ref, path):
    header = ref.header.copy()
    header.set_data_dtype(np.uint8)
    header.set_slope_inter(1, 0)
    out = nib.Nifti1Image(data.astype(np.uint8), ref.affine, header)
    out.set_qform(ref.get_qform(), int(ref.header['qform_code']))
    out.set_sform(ref.get_sform(), int(ref.header['sform_code']))
    out.to_filename(path)


def binary_mask(path, ref, label=None, allow_reorientation=False):
    im = image(path)
    if not same_grid(im, ref) and allow_reorientation:
        orientation = nib.orientations.ornt_transform(nib.orientations.io_orientation(im.affine),
                                                      nib.orientations.io_orientation(ref.affine))
        im = im.as_reoriented(orientation)
    if not same_grid(im, ref):
        raise StageError('Mask must match selected T1 grid; no implicit resampling or affine replacement')
    a = im.get_fdata()
    if label is None:
        if not np.isin(a, [0, 1]).all():
            raise StageError('Mask is not binary; select an explicit positive integer label')
        return a.astype(bool)
    if type(label) is not int or label <= 0 or not np.equal(a, np.floor(a)).all() or (a < 0).any():
        raise StageError('Explicit label requires a nonnegative integer labelmap and a positive label')
    if not np.any(a == label):
        raise StageError('Selected label is absent')
    return a == label


def preparation_spec(c):
    sources = [dict(role='T1w', path=c['reference_t1'])]
    sources += [dict(role=x['modality'], path=x['image']) for x in c.get('review_images', [])]
    if c.get('initial_mask'):
        sources.append(dict(role='initial_mask', path=c['initial_mask']))
    for row in sources:
        row['sha256'] = sha(Path(row['path']))
    return dict(case_id=c['case_id'], lesion_kind=c['lesion_kind'], sources=sources,
                initial_label=c.get('initial_label'), engine_sha256=sha(Path(__file__)),
                editor_engine_sha256=sha(CODE / 'lesion_slicer.py'))


def preparation_candidates(c):
    spec = preparation_spec(c)
    return [spec] + [dict(spec, engine_sha256=engine) for engine in COMPATIBLE_PREPARATION_ENGINES
                     if engine != spec['engine_sha256']]


def validate_preparation_inputs(c):
    ref = image(c['reference_t1'])
    for row in c.get('review_images', []):
        if not same_grid(image(row['image']), ref):
            raise StageError('Review image must already match T1 grid; register explicitly before prepare')
    initial = binary_mask(c['initial_mask'], ref, c.get('initial_label')) if c.get('initial_mask') else np.zeros(ref.shape, bool)
    return ref, initial


def require_prepared(c):
    for spec in preparation_candidates(c):
        key = digest(spec)
        folder = Path(c['workspace']) / 'prepared' / key
        if not folder.exists():
            continue
        if not (folder / 'receipt.json').is_file():
            raise StageError('Incomplete preparation preserved; inspect it before retrying')
        receipt = read(folder / 'receipt.json')
        if receipt['key'] != key or receipt['spec'] != spec:
            raise StageError('Preparation receipt does not match inputs')
        # Legacy recognition never exempts inputs from the new geometry/unit checks.
        validate_preparation_inputs(c)
        verify(folder, receipt['files'])
        return folder, receipt
    raise StageError('Run lesion-prepare for the current inputs first')


def prepare(c, args):
    spec = preparation_spec(c)
    key = digest(spec)
    folder = Path(c['workspace']) / 'prepared' / key
    if any((Path(c['workspace']) / 'prepared' / digest(candidate)).exists() for candidate in preparation_candidates(c)):
        folder, receipt = require_prepared(c)
        return dict(status='REUSED_PREPARED_INPUTS', key=receipt['key'], output=str(folder))
    ref, initial = validate_preparation_inputs(c)
    folder.mkdir(parents=True)
    volumes = []
    for i, row in enumerate(spec['sources']):
        if row['role'] == 'initial_mask':
            continue
        name = f'{i:02d}_{row["role"]}.nii.gz'
        source = Path(row['path'])
        if source.name.endswith('.nii.gz'):
            copy_verified(source, folder / name, row['sha256'])
        elif source.name.endswith('.nii'):
            # Compress original NIfTI bytes, preserving header, scaling, and voxels.
            with source.open('rb') as original, (folder / name).open('wb') as output:
                with gzip.GzipFile(fileobj=output, mode='wb', mtime=0) as compressed:
                    shutil.copyfileobj(original, compressed)
        else:
            raise StageError('Reference/background input must be .nii or .nii.gz')
        if sha(Path(row['path'])) != row['sha256']:
            raise StageError('Input changed during preparation')
        volumes.append(dict(role=row['role'], path=str(folder / name)))
    if c.get('initial_mask') and sha(Path(c['initial_mask'])) != spec['sources'][-1]['sha256']:
        raise StageError('Initial mask changed during preparation')
    save_mask(initial, ref, folder / 'initial_mask.nii.gz')
    session = dict(case_id=c['case_id'], lesion_kind=c['lesion_kind'], prepared_key=key,
                   editor_engine_sha256=spec['editor_engine_sha256'],
                   reference=volumes[0]['path'], volumes=volumes, initial_mask=str(folder / 'initial_mask.nii.gz'),
                   editor_dir=str(Path(c['workspace']) / 'editor' / key))
    write(folder / 'session.json', session)
    write(folder / 'receipt.json', dict(key=key, spec=spec, files=manifest(folder), created_UTC=now()))
    return dict(status='PREPARED_FOR_MANUAL_EDIT', key=key, output=str(folder))


def edit(c, args):
    folder, pre = require_prepared(c)
    executable = getattr(args, 'editor_exe', None) or c['editor'].get('executable')
    if not executable or not Path(executable).is_file():
        raise StageError('Set editor.executable or --editor-exe to an installed Slicer executable')
    import subprocess
    command = [str(executable), '--no-splash', '--python-script', str(CODE / 'lesion_slicer.py'),
               '--', str(folder / 'session.json')]
    if getattr(args, 'dry_run', False):
        return dict(status='EDITOR_PLAN_ONLY', arguments=command)
    process = subprocess.Popen(command, shell=False)
    return dict(status='EDITOR_OPENED_EDIT_NOT_YET_IMPORTED', pid=process.pid, prepared_key=pre['key'])


def qc_images(ref_path, mask, volumes, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ref = image(ref_path)
    # Reorient both the image and label with nibabel; no interpolation or registration.
    orient = nib.orientations.ornt_transform(nib.orientations.io_orientation(ref.affine), nib.orientations.axcodes2ornt(('R', 'A', 'S')))
    canonical_mask = nib.orientations.apply_orientation(mask, orient)
    if mask.any():
        center = np.argwhere(canonical_mask).mean(0).round().astype(int)
    else:
        center = np.array(canonical_mask.shape) // 2
    for i, row in enumerate(volumes):
        a = nib.orientations.apply_orientation(image(row['path']).get_fdata(dtype=np.float32), orient)
        low, high = np.percentile(a, [2, 98])
        if high <= low:
            high = low + 1
        fig, axes = plt.subplots(1, 3, figsize=(10, 4))
        for axis, ax in enumerate(axes):
            bg = np.rot90(np.take(a, center[axis], axis=axis))
            label = np.rot90(np.take(canonical_mask, center[axis], axis=axis))
            # Radiological flips are not introduced: RAS array axes increase R/A/S.
            canonical_zooms = nib.as_closest_canonical(ref).header.get_zooms()
            remaining = [k for k in range(3) if k != axis]
            ax.imshow(bg, cmap='gray', vmin=low, vmax=high,
                      aspect=canonical_zooms[remaining[1]] / canonical_zooms[remaining[0]])
            ax.imshow(np.ma.masked_where(~label, label), cmap='autumn', alpha=.45, vmin=0, vmax=1,
                      aspect=canonical_zooms[remaining[1]] / canonical_zooms[remaining[0]])
            ax.set_title(['Sagittal: horizontal A, top S', 'Coronal: horizontal R, top S', 'Axial: horizontal R, top A'][axis])
            ax.set_xlabel(f'Canonical voxel index {center[axis]}')
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(row['role'] + ' / selected-T1 grid; selected slices only, full editor review required')
        fig.tight_layout()
        fig.savefig(destination / f'{i:02d}_{row["role"]}.png', dpi=120)
        plt.close(fig)


def import_mask(c, args):
    folder, pre = require_prepared(c)
    session = read(folder / 'session.json')
    source = getattr(args, 'mask', None)
    editor_evidence = None
    if source is None:
        pointer = Path(session['editor_dir']) / 'last_draft.json'
        if not pointer.is_file():
            raise StageError('Save a draft in Slicer first, or specify --mask')
        draft = read(pointer)
        if draft['prepared_key'] != pre['key']:
            raise StageError('Draft belongs to different prepared inputs')
        source = draft['mask']
        editor_evidence = {k: draft.get(k) for k in ['editor', 'editor_version', 'editor_script_sha256', 'saved_UTC']}
        p = Path(source).resolve()
        if not p.is_relative_to(Path(session['editor_dir']).resolve()) or sha(p) != draft['mask_sha256']:
            raise StageError('Saved draft mask missing or changed')
    source = Path(source).resolve()
    source_sha = sha(source)
    ref = image(session['reference'])
    allow_reorientation = editor_evidence is not None or getattr(args, 'reorient_to_reference', False)
    if editor_evidence and editor_evidence['editor_script_sha256'] != pre['spec']['editor_engine_sha256']:
        raise StageError('Draft was produced by a different editor bridge')
    source_image = image(source)
    data = binary_mask(source, ref, getattr(args, 'label', None), allow_reorientation)
    lossless_reorientation = not same_grid(source_image, ref)
    if sha(source) != source_sha:
        raise StageError('Mask changed during import')
    spec = dict(prepared_key=pre['key'], source_mask_sha256=source_sha,
                selected_label=getattr(args, 'label', None), lesion_kind=c['lesion_kind'],
                editor_evidence=editor_evidence, lossless_reorientation_applied=lossless_reorientation,
                source_affine=source_image.affine.tolist(), engine_sha256=sha(Path(__file__)))
    key = digest(spec)
    destination = Path(c['workspace']) / 'versions' / key
    if destination.exists():
        verify_version(c, key, pre['key'])
    else:
        destination.mkdir(parents=True)
        save_mask(data, ref, destination / 'lesion_T1w_mask.nii.gz')
        qc = destination / 'qc'; qc.mkdir()
        qc_images(session['reference'], data, session['volumes'], qc)
        count = int(data.sum())
        write(destination / 'provenance.json', dict(spec=spec, source_mask=str(source),
              status='FILES_VALIDATED_REVIEW_PENDING', candidate_voxels=count,
              candidate_volume_ml=count * abs(float(np.linalg.det(ref.affine[:3, :3]))) / 1000,
              mask_empty=count == 0, automatic_brain_clipping=False, registration_performed=False,
              interpretation='File geometry and labels checked; anatomical validity requires human review', created_UTC=now()))
        escaped = html.escape(c['case_id'])
        images = ''.join(f'<figure><img src="qc/{p.name}" alt="{html.escape(p.stem)}"></figure>' for p in sorted(qc.glob('*.png')))
        (destination / 'report.html').write_text(
            '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>病灶导入核查</title>'
            '<style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 20px}img{max-width:100%}code{overflow-wrap:anywhere}</style>'
            f'<h1>{escaped} · 病灶导入核查</h1><p>此页是导入时的静态快照，文件检查通过；审核状态请运行 lesion-status。</p>'
            f'<p>标签类型：{html.escape(c["lesion_kind"])}；阳性体素：{count}；候选体积：{count * abs(float(np.linalg.det(ref.affine[:3,:3]))) / 1000:.6f} mL。</p>'
            '<p>下图只展示选定切片，不能替代逐层、多序列审核。未执行配准、脑掩膜裁剪或模型推理。</p>'
            + images + '<p><a href="lesion_T1w_mask.nii.gz">下载掩膜</a> · <a href="provenance.json">来源记录</a></p></html>', encoding='utf-8')
        write(destination / 'receipt.json', dict(key=key, spec=spec, files=manifest(destination)))
    write(Path(c['workspace']) / 'current.json', dict(version=key, prepared_key=pre['key']))
    return dict(status='FILES_VALIDATED_REVIEW_PENDING', version=key, output=str(destination))


def verify_version(c, key=None, prepared_key=None):
    if key is None:
        key = read(Path(c['workspace']) / 'current.json')['version']
    if not re.fullmatch('[0-9a-f]{64}', key):
        raise StageError('Invalid lesion version')
    folder = Path(c['workspace']) / 'versions' / key
    record = read(folder / 'receipt.json')
    if record['key'] != key or digest(record['spec']) != key:
        raise StageError('Lesion receipt identity mismatch')
    if prepared_key and record['spec']['prepared_key'] != prepared_key:
        raise StageError('Lesion version belongs to different reference inputs')
    verify(folder, record['files'])
    return folder, record


def review(c, args):
    _, pre = require_prepared(c)
    folder, record = verify_version(c, getattr(args, 'version', None), pre['key'])
    reviewer = getattr(args, 'reviewer', None)
    note = getattr(args, 'review_note', None)
    if not reviewer or not reviewer.strip() or not note or not note.strip():
        raise StageError('Explicit --reviewer and --review-note are required')
    im = image(folder / 'lesion_T1w_mask.nii.gz')
    empty = not np.any(im.get_fdata())
    if empty and not getattr(args, 'confirm_empty', False):
        raise StageError('Empty segmentation requires --confirm-empty after visual review')
    decision = dict(status='REVIEWED_MANUAL_ATTESTATION', version=record['key'], prepared_key=pre['key'],
                    mask_sha256=sha(folder / 'lesion_T1w_mask.nii.gz'), reviewer=reviewer.strip(),
                    note=note.strip(), empty_confirmed=empty, reviewed_UTC=now(),
                    limitation='Recorded user attestation, not authentication or clinical certification')
    review_key = digest(decision)
    target = Path(c['workspace']) / 'reviews' / record['key'] / (review_key + '.json')
    write(target, decision)
    write(Path(c['workspace']) / 'review_pointers' / (record['key'] + '.json'),
          dict(path=str(target), sha256=sha(target)))
    return dict(status=decision['status'], version=record['key'], review=str(target))


def require_review(c, record, folder):
    pointer_path = Path(c['workspace']) / 'review_pointers' / (record['key'] + '.json')
    if not pointer_path.is_file():
        raise StageError('This lesion version is not reviewed; run lesion-review first')
    pointer = read(pointer_path)
    path = Path(pointer['path']).resolve()
    base = (Path(c['workspace']) / 'reviews' / record['key']).resolve()
    if not path.is_relative_to(base) or sha(path) != pointer['sha256']:
        raise StageError('Review record missing or changed')
    decision = read(path)
    if (decision['status'] != 'REVIEWED_MANUAL_ATTESTATION' or decision['version'] != record['key']
            or decision['prepared_key'] != record['spec']['prepared_key']
            or decision['mask_sha256'] != sha(folder / 'lesion_T1w_mask.nii.gz')):
        raise StageError('Review does not match the current mask and reference inputs')
    return decision


def export(c, args):
    prepared, pre = require_prepared(c)
    folder, record = verify_version(c, getattr(args, 'version', None), pre['key'])
    decision = require_review(c, record, folder)
    purpose = getattr(args, 'purpose', None)
    if purpose != 'analysis':
        raise StageError('Select --purpose analysis; CFM is deferred in the structural release')
    mask = folder / 'lesion_T1w_mask.nii.gz'
    a = image(mask).get_fdata()
    if not getattr(args, 'export_dir', None):
        raise StageError('Specify a fresh --export-dir')
    destination = Path(args.export_dir).resolve()
    ws = Path(c['workspace']).resolve()
    if destination.is_relative_to(REPOSITORY) or destination == ws or destination in ws.parents:
        raise StageError('Export must not overwrite the repository or workflow workspace')
    if destination.exists():
        raise StageError('Export directory already exists; select a fresh version')
    destination.mkdir(parents=True)
    session = read(prepared / 'session.json')
    copy_verified(mask, destination / mask.name, decision['mask_sha256'])
    write(destination / 'lesion_review.json', decision)
    write(destination / 'provenance.json', dict(case_id=c['case_id'], lesion_kind=c['lesion_kind'], purpose=purpose,
          prepared_key=pre['key'], version=record['key'], reference_t1_sha256=sha(Path(session['reference'])),
          original_reference_t1_sha256=pre['spec']['sources'][0]['sha256'], exported_UTC=now(),
          actual_CFM_execution_verified=False, preprocessing_runs=0, GLM_fits=0))
    write(destination / 'export_receipt.json', dict(files=manifest(destination)))
    return dict(status='REVIEWED_MASK_EXPORTED', purpose=purpose, output=str(destination),
                note='Reviewed mask and provenance only; no preprocessing or statistical analysis launched')


def status(c, args):
    folder, pre = require_prepared(c)
    result = dict(status='PREPARED_FOR_MANUAL_EDIT', prepared_key=pre['key'], output=str(folder))
    pointer = Path(c['workspace']) / 'current.json'
    if pointer.exists():
        key = read(pointer)['version']
        folder, record = verify_version(c, key, pre['key'])
        result.update(status='FILES_VALIDATED_REVIEW_PENDING', version=key, output=str(folder))
        if (Path(c['workspace']) / 'review_pointers' / (key + '.json')).exists():
            result['review'] = require_review(c, record, folder)
            result['status'] = 'REVIEWED_MANUAL_ATTESTATION'
    return result


def dispatch(c, stage, args):
    with lock(c):
        return {'lesion-prepare': prepare, 'lesion-edit': edit, 'lesion-import': import_mask,
                'lesion-review': review, 'lesion-export': export, 'lesion-status': status}[stage](c, args)
