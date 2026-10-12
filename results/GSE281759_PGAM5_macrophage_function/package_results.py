from pathlib import Path
import zipfile,json,hashlib
ROOT=Path(__file__).resolve().parent
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and p.suffix!='.zip' and '__pycache__' not in p.parts and '.venv' not in p.parts and p.name!='package_manifest.json')
manifest=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
(ROOT/'package_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
with zipfile.ZipFile(ROOT/'GSE281759_PGAM5_macrophage_function_results.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in files+[ROOT/'package_manifest.json']:z.write(p,arcname=str(p.relative_to(ROOT)))
print('Packaged',len(files)+1,'files')
