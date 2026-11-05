"""Offline SynthStroke adapter. Produces candidates, never review decisions."""
from pathlib import Path
import os
import subprocess
import uuid

import numpy as np
import lesion_workflow as lesion
from lesion_state import StageError, read, write, sha, now, lock

WORKER = Path(__file__).with_name('synthstroke_worker.py')


def load_config(path, root):
    root = Path(root).resolve()
    def resolve(p):
        return str((root / Path(p)).resolve())
    c = read(resolve(path))
    allowed = {'schema', 'template_only', 'lesion_config', 'workspace', 'python',
               'python_paths', 'model_source', 'model_dir', 'variant', 'device', 'threads', 'patch_size'}
    if not isinstance(c, dict) or set(c) - allowed or c.get('schema') != 'neo-synthstroke-1':
        raise StageError('Unknown SynthStroke schema or keys')
    if c.get('template_only'):
        raise StageError('Fill the SynthStroke template and remove template_only')
    for key in ['lesion_config', 'workspace', 'python', 'model_source', 'model_dir']:
        c[key] = resolve(c[key])
    if not isinstance(c.get('python_paths', []), list):
        raise StageError('python_paths must be a list of external dependency directories')
    c['python_paths'] = [resolve(p) for p in c.get('python_paths', [])]
    c.setdefault('device', 'cpu'); c.setdefault('threads', 6); c.setdefault('patch_size', 128)
    if c.get('variant') not in ['baseline', 'synth'] or c['device'] not in ['cpu', 'cuda']:
        raise StageError('Select baseline/synth and cpu/cuda explicitly; no fallback')
    if any(type(c[k]) is not int or c[k] <= 0 for k in ['threads', 'patch_size']):
        raise StageError('threads and patch_size must be positive integers')
    lc = lesion.load_config(c['lesion_config'], root)
    if lc['lesion_kind'] != 'stroke' or lc.get('initial_mask'):
        raise StageError('Use a stroke lesion configuration with initial_mask null; preserve existing drafts separately')
    ws = Path(c['workspace'])
    if ws.is_relative_to(lesion.REPOSITORY):
        raise StageError('Inference workspace must be outside the code repository')
    sources = [Path(lc['workspace']), Path(lc['reference_t1']), Path(c['lesion_config']),
               Path(c['python']), Path(c['model_source']), Path(c['model_dir'])]
    sources += [Path(x['image']) for x in lc.get('review_images', [])]
    sources += [Path(p) for p in c['python_paths']]
    for p in sources:
        if ws == p or ws in p.parents or p in ws.parents:
            raise StageError('Inference workspace must be separate from lesion workspace, inputs, and runtime assets')
    required = [Path(c['python']), Path(c['model_source']) / 'synthstroke_model.py',
                Path(c['model_dir']) / 'config.json', Path(c['model_dir']) / 'model.safetensors']
    if any(not p.is_file() for p in required) or any(not Path(p).is_dir() for p in c['python_paths']):
        raise StageError('Local Python, model source, dependency directories, config and safetensors must exist')
    mc = read(required[-2])
    if (mc.get('spatial_dims'), mc.get('in_channels'), mc.get('out_channels'), mc.get('model_type')) != (3, 1, 2 if c['variant'] == 'baseline' else 6, c['variant']):
        raise StageError('Only matching single-channel baseline (2 outputs) or synth (6 outputs) checkpoints are supported')
    ref = lesion.image(lc['reference_t1'])
    for row in lc.get('review_images', []):
        if not lesion.same_grid(lesion.image(row['image']), ref):
            raise StageError('Review images must already match selected T1 grid')
    assets = required + sorted(Path(c['model_source']).rglob('*.py')) + [WORKER, Path(__file__)]
    assets += [Path(lc['reference_t1'])] + [Path(x['image']) for x in lc.get('review_images', [])]
    c['_hashes'] = {str(p): sha(p) for p in assets}
    c['_lesion'] = lc
    return c


def run(c, dry_run=False):
    plan = dict(backend='synthstroke', variant=c['variant'], device=c['device'],
                reference_t1=c['_lesion']['reference_t1'], python=c['python'],
                model_dir=c['model_dir'], workspace=c['workspace'], review_required=True,
                input_sha256=c['_hashes'], downloads=False, other_models=False)
    if dry_run:
        return dict(status='PLANNED_SYNTHSTROKE_NO_EXECUTION', **plan)
    with lock(c):
        folder = Path(c['workspace']) / uuid.uuid4().hex
        folder.mkdir(parents=True)
        job = {k: c[k] for k in ['python_paths', 'model_source', 'model_dir', 'variant', 'device', 'threads', 'patch_size']}
        job.update(reference_t1=plan['reference_t1'], output=str(folder))
        write(folder / 'job.json', job)
        write(folder / 'request.json', dict(**plan, started_UTC=now()))
        env = os.environ.copy()
        env.update(HF_HUB_OFFLINE='1', HF_HUB_DISABLE_TELEMETRY='1', WANDB_DISABLED='true', PYTHONNOUSERSITE='1')
        # Avoid inheriting the workflow's dependency search path into the model runtime.
        env.pop('PYTHONPATH', None); env.pop('PYTHONHOME', None)
        command = [c['python'], '-s', str(WORKER), str(folder / 'job.json')]
        try:
            with (folder / 'inference.log').open('w', encoding='utf-8') as log:
                result = subprocess.run(command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise StageError('SynthStroke worker failed; inspect ' + str(folder / 'inference.log'))
            for path, expected in c['_hashes'].items():
                if sha(Path(path)) != expected:
                    raise StageError('Input or runtime asset changed during inference: ' + path)
            mask = folder / 'candidate_T1w_mask.nii.gz'
            ref = lesion.image(plan['reference_t1'])
            lesion.binary_mask(mask, ref)
            score = lesion.image(folder / 'lesion_softmax_score.nii.gz')
            if not lesion.same_grid(score, ref) or not np.all((score.get_fdata() >= 0) & (score.get_fdata() <= 1)):
                raise StageError('Invalid native T1 softmax score')
            handoff = dict(c['_lesion'], initial_mask=str(mask))
            write(folder / 'lesion.with-candidate.json', handoff)
            write(folder / 'receipt.json', dict(status='CANDIDATE_REVIEW_PENDING', command=command,
                  request=plan, worker=read(folder / 'worker_result.json'),
                  files=lesion.manifest(folder), completed_UTC=now()))
        except Exception as exc:
            write(folder / 'failed.json', dict(status='FAILED_NO_REVIEW', reason=str(exc)))
            raise
        return dict(status='CANDIDATE_REVIEW_PENDING', output=str(folder), mask=str(mask),
                    lesion_config=str(folder / 'lesion.with-candidate.json'),
                    next_stage='lesion-prepare, then lesion-edit with the generated lesion config')
