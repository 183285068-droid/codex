from pathlib import Path
import hashlib,json,zipfile
R=Path(__file__).resolve().parent
z=R/'GSE283406_PGAM5_Tcell_expression_all_results.zip'
files=[p for p in sorted(R.rglob('*')) if p.is_file() and p!=z and p.name!='package_manifest.json' and '__pycache__' not in p.parts]
records=[dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]
(R/'package_manifest.json').write_text(json.dumps(records,indent=2))
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
 for p in files+[R/'package_manifest.json']:f.write(p,p.relative_to(R))
with zipfile.ZipFile(z) as f:
 assert f.testzip() is None
 for v in records:assert hashlib.sha256(f.read(v['path'])).hexdigest()==v['sha256']
print('Verified',len(records),'files; ZIP bytes',z.stat().st_size)
