from pathlib import Path
import json,hashlib,zipfile
ROOT=Path(__file__).resolve().parent
p=ROOT/'GSE281759_SPP1_macrophage_T_interaction_results.zip'
files=[f for f in ROOT.rglob('*') if f.is_file() and '__pycache__' not in f.parts and 'work' not in f.relative_to(ROOT).parts and f.suffix!='.zip' and f.name!='package_manifest.json' and not f.name.startswith('.')]
with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for f in sorted(files):z.write(f,arcname=ROOT.name+'/'+str(f.relative_to(ROOT)))
with zipfile.ZipFile(p) as z:assert z.testzip() is None
assert p.stat().st_size<30*1024**2
manifest=dict(name=p.name,bytes=p.stat().st_size,MiB=round(p.stat().st_size/1024**2,2),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),files=len(files),zip_crc_passed=True)
(ROOT/'package_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(manifest)
