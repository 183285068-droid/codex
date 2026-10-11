"""Fetch only frozen prerequisite files from the already-published base analysis."""
from pathlib import Path
import requests,json,hashlib,concurrent.futures
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent/'GSE245908_spatial_PGAM5';COMMIT='9cccc00079d26e55eb9d4f1d746cca9394366b76'
def get(record):
 p=BASE/record['path'];p.parent.mkdir(parents=True,exist_ok=True)
 if not p.exists():
  url=f'https://raw.githubusercontent.com/183285068-droid/codex/{COMMIT}/results/GSE245908_spatial_PGAM5/{record["path"]}'
  with requests.get(url,timeout=(20,120),stream=True) as r:
   r.raise_for_status()
   with p.with_suffix(p.suffix+'.tmp').open('wb') as f:
    for chunk in r.iter_content(1024*1024):f.write(chunk)
  p.with_suffix(p.suffix+'.tmp').rename(p)
 assert p.stat().st_size==record['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==record['sha256'],f'Frozen input mismatch: {p}'
 return record['path']
if __name__=='__main__':
 records=json.loads((ROOT/'source_inputs/base_requirements.json').read_text())
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
  for path in ex.map(get,records):print('Verified',path)
