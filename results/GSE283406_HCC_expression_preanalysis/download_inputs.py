from pathlib import Path
import concurrent.futures
import hashlib
import json
import os
import requests

ROOT = Path(__file__).resolve().parent
WORK = Path(os.environ.get('GSE283406_WORK', ROOT/'work'))
INPUT = WORK/'inputs'
INPUT.mkdir(parents=True, exist_ok=True)

def download(item):
    path = INPUT/item['name']
    if not path.exists():
        response = requests.get(item['url'], stream=True, timeout=(20,120))
        response.raise_for_status()
        temporary = path.with_suffix(path.suffix+'.tmp')
        with temporary.open('wb') as handle:
            for chunk in response.iter_content(1024*1024):
                handle.write(chunk)
        temporary.rename(path)
    assert path.stat().st_size == item['bytes'], path
    assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256'], path
    return path.name

if __name__ == '__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for name in pool.map(download, json.loads((ROOT/'source_manifest.json').read_text())):
            print('Verified:', name, flush=True)
