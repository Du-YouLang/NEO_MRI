from pathlib import Path
import sys,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from state_io import read_json,write_json

def test_transient_state_read_is_retried_without_running_science(tmp_path,monkeypatch):
    p=tmp_path/'state.json';p.write_text('{"state":"DONE"}')
    original=Path.read_text;calls=[]
    def read(self,**kw):
        calls.append(self)
        if len(calls)==1:raise PermissionError('synthetic sharing violation')
        return original(self,**kw)
    monkeypatch.setattr(Path,'read_text',read)
    assert read_json(p,attempts=2,delay=0)==dict(state='DONE') and len(calls)==2

def test_persistent_state_error_stops_after_bound(tmp_path,monkeypatch):
    calls=[]
    def read(*a,**kw):calls.append(1);raise PermissionError('durable denial, never bypassed')
    monkeypatch.setattr(Path,'read_text',read)
    with pytest.raises(PermissionError):read_json(tmp_path/'state.json',attempts=3,delay=0)
    assert len(calls)==3

def test_replace_sharing_conflict_is_bounded(tmp_path,monkeypatch):
    p=tmp_path/'state.json';p.write_text('{"state":"OLD"}')
    replace=Path.replace;calls=[]
    def changed(self,target):
        calls.append(self)
        if len(calls)==1:raise PermissionError('synthetic transient replace conflict')
        return replace(self,target)
    monkeypatch.setattr(Path,'replace',changed)
    write_json(p,dict(state='NEW'),attempts=2,delay=0)
    assert read_json(p)==dict(state='NEW') and len(calls)==2
