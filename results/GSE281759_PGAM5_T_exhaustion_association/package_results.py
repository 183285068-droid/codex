from pathlib import Path
import zipfile,hashlib,json
ROOT=Path(__file__).resolve().parent
name='GSE281759_PGAM5_T_exhaustion_association_results.zip';p=ROOT/name
files=[f for f in ROOT.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.suffix!='.zip' and f.name!='package_manifest.json' and not f.name.startswith('.')]
with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for f in sorted(files):z.write(f,arcname='GSE281759_PGAM5_T_exhaustion_association/'+str(f.relative_to(ROOT)))
assert p.stat().st_size<30*1024**2
with zipfile.ZipFile(p) as z:assert z.testzip() is None
manifest=dict(name=name,bytes=p.stat().st_size,MiB=round(p.stat().st_size/1024**2,2),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),files=len(files),zip_crc_passed=True)
(ROOT/'package_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(manifest)
