"""Run only in the explicitly selected external SynthStroke Python environment."""
import copy
import importlib.metadata as md
import json
from pathlib import Path
import sys


def main(job_path):
    job = json.loads(Path(job_path).read_text(encoding='utf-8'))
    sys.path[:0] = [job['model_source']] + job['python_paths']
    import numpy as np
    import nibabel as nib
    import torch
    import monai as mn
    from monai.data import MetaTensor
    from safetensors.torch import load_file
    from synthstroke_model import SynthStrokeModel, ModelConfig

    torch.set_num_threads(job['threads']); torch.manual_seed(0)
    if job['device'] == 'cuda' and not torch.cuda.is_available():
        raise ValueError('CUDA explicitly selected but unavailable')
    cfg = ModelConfig.from_dict(json.loads((Path(job['model_dir']) / 'config.json').read_text(encoding='utf-8')))
    model = SynthStrokeModel(cfg)
    state = load_file(str(Path(job['model_dir']) / 'model.safetensors'))
    if not any(k.startswith('unet.') for k in state):
        state = {'unet.' + k: v for k, v in state.items()}
    model.load_state_dict(state, strict=True)
    model.to(job['device']).eval()
    del state
    # MONAI retains native affine BEFORE RAS and spacing transforms. Invert logits,
    # then apply softmax/argmax in native T1 space (never interpolate a binary mask).
    prep = mn.transforms.Compose([
        mn.transforms.LoadImaged(keys='img', image_only=True),
        mn.transforms.EnsureChannelFirstd(keys='img'),
        mn.transforms.Orientationd(keys='img', axcodes='RAS'),
        mn.transforms.Spacingd(keys='img', pixdim=(1, 1, 1), mode='bilinear'),
        mn.transforms.HistogramNormalized(keys='img'),
        mn.transforms.NormalizeIntensityd(keys='img', nonzero=False, channel_wise=True)])
    ref = nib.load(job['reference_t1'])
    if ref.header.get_xyzt_units()[0] != 'mm':
        raise ValueError('SynthStroke requires NIfTI spatial units explicitly set to mm')
    image = prep({'img': job['reference_t1']})['img'].to(dtype=torch.float32)
    with torch.inference_mode():
        logits = mn.inferers.sliding_window_inference(image[None], roi_size=(job['patch_size'],) * 3,
                 sw_batch_size=1, predictor=model, overlap=.5, mode='gaussian',
                 sw_device=job['device'], device='cpu')
        output = MetaTensor(logits[0].cpu(), affine=image.affine, meta=copy.deepcopy(image.meta),
                            applied_operations=copy.deepcopy(image.applied_operations))
        with mn.transforms.utils.allow_missing_keys_mode(prep):
            native = prep.inverse({'img': output})['img']
        probs = torch.softmax(native, dim=0).cpu().numpy()
    if probs.shape[1:] != ref.shape or not np.allclose(native.affine.numpy(), ref.affine, rtol=0, atol=1e-4):
        raise ValueError('Inverse transform did not recover native T1 geometry')
    if not np.isfinite(probs).all() or not np.allclose(probs.sum(0), 1, atol=1e-5):
        raise ValueError('Nonfinite or invalid softmax output')
    label = 1 if job['variant'] == 'baseline' else 5
    mask = (probs.argmax(0) == label).astype(np.uint8)
    for name, data, dtype in [('candidate_T1w_mask.nii.gz', mask, np.uint8),
                              ('lesion_softmax_score.nii.gz', probs[label], np.float32)]:
        header = ref.header.copy(); header.set_data_dtype(dtype); header.set_slope_inter(1, 0)
        out = nib.Nifti1Image(data.astype(dtype), ref.affine, header)
        out.set_qform(ref.get_qform(), int(ref.header['qform_code']))
        out.set_sform(ref.get_sform(), int(ref.header['sform_code']))
        out.to_filename(Path(job['output']) / name)
    result = dict(python=sys.executable, versions={p: md.version(p) for p in
                  ['torch', 'monai', 'numpy', 'nibabel', 'safetensors', 'huggingface_hub']},
                  variant=job['variant'], stroke_label=label, native_grid_verified=True,
                  candidate_voxels=int(mask.sum()), inference=dict(patch_size=job['patch_size'],
                  overlap=.5, blending='gaussian logits', TTA=False, device=job['device']),
                  preprocessing=['full-head T1', 'RAS', '1mm', 'histogram normalization', 'whole-image z normalization'],
                  selection='native-space multiclass argmax; no brain clipping or component removal',
                  note='Candidate only; softmax score is not calibrated confidence or clinical evidence')
    (Path(job['output']) / 'worker_result.json').write_text(json.dumps(result, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main(sys.argv[1])
