"""Capture actual Nipype input-chain evidence inside the pinned fMRIPrep image."""
from pathlib import Path
import argparse,hashlib,json,shutil
import nibabel as nib
import numpy as np

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def record(p):return dict(path=str(p),sha256=sha(p))
def same_values(a,b):
    a=nib.load(a);b=nib.load(b)
    return a.shape==b.shape and np.allclose(a.affine,b.affine,atol=1e-4,rtol=0) and np.array_equal(a.get_fdata(dtype=np.float32),b.get_fdata(dtype=np.float32))
def capture(work,approved_t1,approved_mask,output,review_sha):
    from nipype.utils.filemanip import loadpkl
    work=Path(work);output=Path(output)
    if output.exists() and any(output.iterdir()):raise FileExistsError('Fresh runtime evidence required')
    output.mkdir(parents=True,exist_ok=True)
    found={}
    for f in work.rglob('result_*.pklz'):
        if f.parent.name not in ['apply_mask','trunc_mov','registration']:continue
        result=loadpkl(str(f));inputs=result.inputs
        if f.parent.name=='apply_mask' and isinstance(inputs,dict) and 'in_mask' in inputs:
            if sha(inputs['in_mask'])==sha(approved_mask):found['apply']=(f,result)
        elif f.parent.name=='trunc_mov' and isinstance(inputs,dict) and 'op1' in inputs:
            if Path(str(inputs['op1'])).is_file() and same_values(inputs['op1'],approved_t1):found['trunc']=(f,result)
        elif f.parent.name=='registration' and isinstance(inputs,dict) and 'moving_mask' in inputs:
            m=inputs['moving_mask']
            if isinstance(m,str) and Path(m).is_file() and sha(m)==sha(approved_mask):found['registration']=(f,result)
    if set(found)!={'apply','trunc','registration'}:raise RuntimeError('Actual mask/normalization input chain missing: '+str(list(found)))
    apply=found['apply'][1];trunc=found['trunc'][1];reg=found['registration'][1]
    if not same_values(apply.inputs['in_file'],approved_t1):raise RuntimeError('Mask application used another T1')
    if sha(trunc.outputs.output_image)!=sha(reg.inputs['moving_image']):raise RuntimeError('Normalization moving input not linked to observed truncation')
    a=nib.load(apply.inputs['in_file']);mask=nib.load(apply.inputs['in_mask']);brain=nib.load(apply.outputs.out_file)
    if a.shape!=brain.shape or not np.allclose(a.affine,brain.affine,atol=1e-4,rtol=0) or not np.allclose(a.get_fdata()*mask.get_fdata(),brain.get_fdata(),atol=1e-5,rtol=1e-6):raise RuntimeError('Observed brain image differs from mask application')
    outputs={}
    for key,source in [('observed_mask',reg.inputs['moving_mask']),('observed_t1',apply.inputs['in_file']),('masked_brain',apply.outputs.out_file)]:
        target=output/(key+'.nii.gz');shutil.copy2(source,target);outputs[key]=record(target)
    sources=[]
    for key,(source,result) in found.items():
        target=output/(key+'_result.pklz');shutil.copy2(source,target);sources.append(record(target))
    receipt=dict(schema='neo-mask-runtime-evidence-1',review_sha256=review_sha,normalization_input_chain_verified=True,
        expected_mask_sha256=sha(approved_mask),node_result_files=sources,**outputs,
        observed_apply_mask=dict(in_mask=apply.inputs['in_mask'],in_file=apply.inputs['in_file'],out_file=apply.outputs.out_file),
        observed_registration=dict(moving_mask=reg.inputs['moving_mask'],moving_image=reg.inputs['moving_image'],truncation_input=trunc.inputs['op1']),
        note='Observed reused inputs and actual mask application, not newly emitted native derivatives; no clinical approval.')
    (output/'runtime_evidence.json').write_text(json.dumps(receipt,indent=2))
    return receipt
def main():
    p=argparse.ArgumentParser();p.add_argument('--work',required=True);p.add_argument('--t1',required=True);p.add_argument('--mask',required=True);p.add_argument('--output',required=True);p.add_argument('--review-sha',required=True);a=p.parse_args()
    capture(a.work,a.t1,a.mask,a.output,a.review_sha)
if __name__=='__main__':main()
