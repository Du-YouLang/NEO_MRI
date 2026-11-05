"""Bounded control-file I/O retries; never retries scientific computation."""
from pathlib import Path
import json,os,time

def read_json(path,attempts=20,delay=.1):
    if not 1<=attempts<=100 or not 0<=delay<=1:raise ValueError('Bounded state-read policy required')
    for attempt in range(attempts):
        try:return json.loads(Path(path).read_text(encoding='utf-8'))
        except (PermissionError,FileNotFoundError,json.JSONDecodeError):
            if attempt+1==attempts:raise
            time.sleep(delay)

def write_json(path,value,attempts=20,delay=.1):
    if not 1<=attempts<=100 or not 0<=delay<=1:raise ValueError('Bounded state-write policy required')
    path=Path(path);tmp=path.with_name(path.name+'.'+str(os.getpid())+'.tmp')
    tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    for attempt in range(attempts):
        try:tmp.replace(path);return
        except PermissionError:
            if attempt+1==attempts:raise
            time.sleep(delay)
