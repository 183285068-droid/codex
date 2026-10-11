"""Keep downloadable archives below 30 MiB; verify ZIP CRC and record SHA256."""
from pathlib import Path
import hashlib,json,zipfile
ROOT=Path(__file__).resolve().parent
SAMPLES=json.loads((ROOT/'sample_scope.json').read_text());LIMIT=28*1024**2
base=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.zip' and p.name!='package_manifest.json' and not p.name.startswith('.') and 'work' not in p.relative_to(ROOT).parts and 'inputs' not in p.relative_to(ROOT).parts]
assign={};used=set()
for patient in ['A009','B013','C004','C008']:
    groups=[]
    for item in SAMPLES:
        if item['patient']!=patient:continue
        files=[p for p in base if p.name.startswith(item['sample']+'_')];groups.append(files);used.update(files)
    batches=[];batch=[];size=0
    for files in groups:
        s=sum(p.stat().st_size for p in files)
        if size+s>LIMIT and batch:batches.append(batch);batch=[];size=0
        batch.extend(files);size+=s
    if batch:batches.append(batch)
    for k,files in enumerate(batches,1):assign[f'GSE281759_{patient}_results'+(f'_{k:02d}' if len(batches)>1 else '')+'.zip']=files
assign['GSE281759_overview_methods_results.zip']=[p for p in base if p not in used]
readme=ROOT/'README.md';text=readme.read_text().split('\n## 下载结果包')[0]
text+='\n## 下载结果包\n\n所有包各自低于30MiB；打开GitHub文件后点击 **Download raw file**。公共包包含总览、方法和来源，各患者包包含对应逐片PNG/PDF、原始低分辨组织图、注释/标志/候选边表及便携AnnData。\n\n'
for name in ['GSE281759_overview_methods_results.zip']+[n for n in assign if n!='GSE281759_overview_methods_results.zip']:text+=f'- [{name}]({name})\n'
text+='\n[包大小与SHA256](package_manifest.json)。原始完整矩阵/高分辨图请用download_inputs.py重新获取。\n';readme.write_text(text)
manifest=[]
for name,files in assign.items():
    p=ROOT/name
    with zipfile.ZipFile(p,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for f in sorted(files):z.write(f,arcname=str(f.relative_to(ROOT)))
    assert p.stat().st_size<30*1024**2,(name,p.stat().st_size)
    with zipfile.ZipFile(p) as z:assert z.testzip() is None
    manifest.append(dict(name=name,bytes=p.stat().st_size,MiB=round(p.stat().st_size/1024**2,2),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),files=len(files)))
    print(name,manifest[-1]['MiB'],'MiB',flush=True)
(ROOT/'package_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
