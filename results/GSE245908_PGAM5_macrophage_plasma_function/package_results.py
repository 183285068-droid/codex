from pathlib import Path
import json,hashlib,zipfile
ROOT=Path(__file__).resolve().parent
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and p.suffix!='.zip' and '__pycache__' not in p.parts and '.venv' not in p.parts and p.name!='package_manifest.json')
records=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]
(ROOT/'package_manifest.json').write_text(json.dumps(records,indent=2)+'\n');archive=ROOT/'GSE245908_PGAM5_macrophage_plasma_function_all_results.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
 for p in files+[ROOT/'package_manifest.json']:z.write(p,arcname=str(p.relative_to(ROOT)))
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for row in records:assert hashlib.sha256(z.read(row['path'])).hexdigest()==row['sha256']
print('Verified',len(files)+1,'files; ZIP MB',round(archive.stat().st_size/1e6,2))
