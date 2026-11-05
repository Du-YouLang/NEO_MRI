"""Check built artifacts, including accidentally retained build/lib modules."""
from pathlib import Path
import sys,tomllib,zipfile
root=Path(__file__).resolve().parents[1]
allowed=set(tomllib.loads((root/'pyproject.toml').read_text(encoding='utf-8'))['tool']['setuptools']['py-modules'])
for path in sys.argv[1:]:
    with zipfile.ZipFile(path) as z:
        modules={Path(n).stem for n in z.namelist() if n.endswith('.py')}
        if modules!=allowed:raise ValueError('Wheel module scope mismatch: '+str(modules^allowed))
        if any(n.endswith(('.nii','.nii.gz','.dcm','.h5','.ipynb')) for n in z.namelist()):raise ValueError('Non-code payload in wheel')
        ep=z.read('neo_mri-0.2.0rc1.dist-info/entry_points.txt').decode()
        if 'neo_mri_cli' in ep or 'staged_glm' in ep:raise ValueError('Functional entrypoint in structural wheel')
    print('Verified structural wheel:',path,'modules:',len(modules))
