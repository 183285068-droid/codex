"""Fetch official inputs with TLS verification, retries and exact byte/SHA256 checks."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
ROOT=Path(__file__).resolve().parent
WORK=Path(os.environ.get('GSE281759_WORK',ROOT/'work'))
INPUT=WORK/'inputs';INPUT.mkdir(parents=True,exist_ok=True)
def verified(p,item):
    return p.exists() and p.stat().st_size==item['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
def download(item):
    p=INPUT/item['name']
    if not verified(p,item):
        with requests.Session() as session:
            session.mount('https://',HTTPAdapter(max_retries=Retry(total=4,backoff_factor=1,status_forcelist=[429,500,502,503,504])))
            with session.get(item['url'],stream=True,timeout=(20,180)) as response:
                response.raise_for_status();tmp=p.with_suffix(p.suffix+'.tmp')
                with tmp.open('wb') as handle:
                    for chunk in response.iter_content(1024*1024):handle.write(chunk)
            assert verified(tmp,item),item['name']
            tmp.replace(p)
    return p.name
if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=2) as pool:
        for name in pool.map(download,json.loads((ROOT/'source_manifest.json').read_text())):print('Verified',name,flush=True)
