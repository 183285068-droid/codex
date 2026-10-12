from pathlib import Path
import hashlib,json,zipfile
R=Path(__file__).resolve().parent
z=R/'GSE281759_PGAM5_Tcell_extended_all_results.zip'
files=[p for p in sorted(R.rglob('*')) if p.is_file() and p!=z and p.name!='package_manifest.json' and '__pycache__' not in p.parts]
manifest=[dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]
(R/'package_manifest.json').write_text(json.dumps(manifest,indent=2))
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
 for p in files+[R/'package_manifest.json']:f.write(p,p.relative_to(R))
with zipfile.ZipFile(z) as f:
 assert f.testzip() is None
 for r in manifest:assert hashlib.sha256(f.read(r['path'])).hexdigest()==r['sha256']
print('Verified',len(manifest),'files; ZIP bytes',z.stat().st_size)
