"""Structural-only path and protocol contract; no functional-analysis dependency."""
from pathlib import Path
from anatomy_review import sha,read
from state_io import write_json as write
ROOT=Path.cwd()
class StageError(RuntimeError):pass
def resolve(path):
    p=Path(path);return p.resolve() if p.is_absolute() else (ROOT/p).resolve()
def source_root(c):return resolve(c['preprocessing']['project'])
def load(path):
    global ROOT
    p=Path(path).resolve();ROOT=p.parent;c=read(p)
    if c.get('template_only'):raise StageError('Fill a private copy and remove template_only before use')
    if c.get('schema')!='neo-t1-anatomy-1':raise StageError('Structural T1 configuration required')
    return c
def arguments(p):
    if p.get('normalization')!='standard' or p.get('sdc')!='off' or p.get('stc')!='off':
        raise StageError('Structural release supports standard T1 normalization only, no CFM/SDC/STC')
    if not p.get('t1w') or p.get('fieldmaps') or p.get('bold_metadata'):raise StageError('T1-only protocol required')
    return ['--ignore','slicetiming','fieldmaps','flair']
def validate_inputs(root,p):
    arguments(p);file=(Path(root)/'bids'/p['t1w']).resolve()
    if not file.is_relative_to((Path(root)/'bids').resolve()) or not file.is_file():raise StageError('Selected raw T1 missing/outside BIDS')
    import nibabel as nib
    if len(nib.load(file).shape)!=3:raise StageError('T1 must be a 3D scalar image')
