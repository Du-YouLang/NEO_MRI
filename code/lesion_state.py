"""Small state primitives for the independent lesion workflow."""
from contextlib import contextmanager
from pathlib import Path
import datetime
import hashlib
import json
import os


class StageError(RuntimeError):
    pass


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':'),
                                   ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


@contextmanager
def lock(c):
    path = Path(c['workspace']) / '.stage.lock'
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        stream = path.open('x', encoding='utf-8')
    except FileExistsError:
        raise StageError(f'Another stage owns {path}; inspect its PID before removing a stale lock')
    try:
        with stream:
            json.dump(dict(pid=os.getpid(), started_UTC=now()), stream)
        yield
    finally:
        path.unlink(missing_ok=True)
