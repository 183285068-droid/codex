from pathlib import Path
import zipfile,json,hashlib
ROOT=Path(__file__).resolve().parent
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and p.suffix!='.zip' and '__pycache__' not in p.parts and '.venv' not in p.parts and p.name!='package_manifest.json')
rows=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files];(ROOT/'package_manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
archive=ROOT/'GSE281759_PGAM5_plasma_patient4_all_results.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for p in files+[ROOT/'package_manifest.json']:z.write(p,arcname=str(p.relative_to(ROOT)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in rows:assert hashlib.sha256(z.read(r['path'])).hexdigest()==r['sha256']
assert archive.stat().st_size<100*1024**2
print('Verified',len(files)+1,'files; ZIP MB',round(archive.stat().st_size/1e6,2))
