from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import os,json,hashlib,requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
ROOT=Path(__file__).resolve().parent;WORK=Path(os.environ.get('GSE203612_WORK',ROOT/'work'));INPUT=WORK/'inputs';INPUT.mkdir(parents=True,exist_ok=True)
def download(it):
    p=INPUT/it['name']
    if not p.exists():
        with requests.Session() as sess:
            sess.mount('https://',HTTPAdapter(max_retries=Retry(total=4,backoff_factor=1,status_forcelist=[429,500,502,503,504])))
            r=sess.get(it['url'],timeout=(20,180));r.raise_for_status();tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_bytes(r.content);assert tmp.stat().st_size==it['bytes'];tmp.replace(p)
    assert p.stat().st_size==it['bytes'];sha=hashlib.sha256(p.read_bytes()).hexdigest()
    if 'sha256' in it:assert sha==it['sha256']
    return dict(**it,sha256=sha) if 'sha256' not in it else it
if __name__=='__main__':
    source=ROOT/'source_manifest.json';plan=json.loads((source if source.exists() else ROOT/'download_plan.json').read_text())
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(download,plan))
    source.write_text(json.dumps(rows,indent=2)+'\n');print('Verified',len(rows),'files',flush=True)
