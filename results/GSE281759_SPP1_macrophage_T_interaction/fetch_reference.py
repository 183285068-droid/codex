"""Fetch pinned prior marker AnnData objects; preserve existing files and verify SHA256."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,hashlib,requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
ROOT=Path(__file__).resolve().parent;REF=ROOT.parent/'GSE281759_spatial_HCC_PGAM5';REF.mkdir(exist_ok=True)
def get(item):
    p=REF/item['name']
    if not p.exists():
        with requests.Session() as sess:
            sess.mount('https://',HTTPAdapter(max_retries=Retry(total=4,backoff_factor=1,status_forcelist=[429,500,502,503,504])))
            resp=sess.get(item['url'],timeout=(20,180));resp.raise_for_status();tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_bytes(resp.content)
            assert tmp.stat().st_size==item['bytes'] and hashlib.sha256(tmp.read_bytes()).hexdigest()==item['sha256']
            tmp.replace(p)
    assert p.stat().st_size==item['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'],p
    return p.name
if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=2) as pool:
        for name in pool.map(get,json.loads((ROOT/'reference_manifest.json').read_text())):print('Verified reference',name,flush=True)
