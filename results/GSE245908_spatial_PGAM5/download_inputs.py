"""Download public GEO data; never execute metadata as instructions."""
from pathlib import Path
import requests, hashlib, json, concurrent.futures, argparse
ROOT=Path(__file__).resolve().parent
SAMPLES={'CHC20':'GSM7850822','CHC23':'GSM7850823'}
EXPECTED={r['path']:r for r in json.loads((ROOT/'source_manifest.json').read_text())} if (ROOT/'source_manifest.json').exists() else {}
def verify_record(record):
 expected=EXPECTED.get(record['path'])
 if expected and (expected['sha256']!=record['sha256'] or expected['bytes']!=record['bytes']):
  raise ValueError(f"Source changed or file is incomplete: {record['path']}; compare against frozen source_manifest.json before rerunning.")
 return record
def download(sample,gsm,suffix):
 name=f'{gsm}_{sample}_{suffix}';dest=ROOT/'inputs'/sample/name;dest.parent.mkdir(parents=True,exist_ok=True)
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7850nnn/{gsm}/suppl/{name}'
 if not dest.exists():
  with requests.get(url,timeout=(20,120),stream=True) as r:
   r.raise_for_status()
   with dest.with_suffix(dest.suffix+'.tmp').open('wb') as f:
    for chunk in r.iter_content(1024*1024):f.write(chunk)
  dest.with_suffix(dest.suffix+'.tmp').rename(dest)
 return verify_record({'sample':sample,'gsm':gsm,'path':str(dest.relative_to(ROOT)),'url':url,'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
if __name__=='__main__':
 (ROOT/'inputs').mkdir(exist_ok=True)
 tasks=[(s,g,f) for s,g in SAMPLES.items() for f in ['filtered_feature_bc_matrix.h5','scalefactors_json.json.gz','tissue_hires_image.png.gz','tissue_lowres_image.png.gz','tissue_positions_list.csv.gz','detected_tissue_image.jpg.gz']]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:manifest=list(ex.map(lambda x:download(*x),tasks))
 for acc in ['GSE245908',*SAMPLES.values()]:
  u=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=full';d=ROOT/'inputs'/f'{acc}_metadata.soft'
  if not d.exists():
   r=requests.get(u,timeout=60);r.raise_for_status();d.write_bytes(r.content)
  content=d.read_bytes();manifest.append(verify_record({'path':str(d.relative_to(ROOT)),'url':u,'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest()}))
 u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245908/suppl/GSE245908_Info_Supp_for_Visium.tsv.gz';d=ROOT/'inputs'/u.rsplit('/',1)[-1]
 if not d.exists():
  r=requests.get(u,timeout=60);r.raise_for_status();d.write_bytes(r.content)
 content=d.read_bytes();manifest.append(verify_record({'path':str(d.relative_to(ROOT)),'url':u,'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest()}))
 (ROOT/'source_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
