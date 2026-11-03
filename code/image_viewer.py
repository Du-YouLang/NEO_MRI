"""Generic saved-image HTML viewer. Never registers, fits, or infers space."""
from pathlib import Path
import argparse
import nibabel as nib
import numpy as np

def export(background,output,overlay=None,space='UNDECLARED',threshold=2.3):
    from nilearn import plotting
    bg=nib.load(background)
    if len(bg.shape)!=3:raise ValueError('Background must be 3D')
    if not np.isfinite(bg.affine).all():raise ValueError('Invalid affine')
    if overlay:
        image=nib.load(overlay)
        if image.shape!=bg.shape or not np.allclose(image.affine,bg.affine,rtol=0,atol=1e-5):
            raise ValueError('Overlay must match background grid; resample with verified transforms before viewing')
        if threshold<0:raise ValueError('Nonnegative display threshold required')
        view=plotting.view_img(image,bg_img=bg,threshold=threshold,symmetric_cmap=True,cmap='cold_hot',
                               title=f'{space}: saved statistic (display threshold only)')
    else:
        view=plotting.view_img(bg,bg_img=False,threshold=None,symmetric_cmap=False,cmap='gray',title=f'{space}: anatomy')
    output=Path(output)
    if output.exists():raise FileExistsError('Choose a new viewer filename')
    output.parent.mkdir(parents=True,exist_ok=True);view.save_as_html(str(output));return output

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--background',required=True);parser.add_argument('--output',required=True)
    parser.add_argument('--overlay');parser.add_argument('--space',required=True,help='Explicit native T1/MNI template identity')
    parser.add_argument('--threshold',type=float,default=2.3)
    args=parser.parse_args(argv)
    print(export(args.background,args.output,args.overlay,args.space,args.threshold));return 0
