from pathlib import Path
import urllib.request,hashlib
R=Path(__file__).resolve().parent;p=R/'GSE151530_matrix.mtx.gz';expected='50dad60e779b3332350e2552ed4260369ccba8db4f73918bbe2a9f6549db6403'
if not p.exists():
 url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_matrix.mtx.gz';print('Downloading official GEO matrix',flush=True)
 with urllib.request.urlopen(url,timeout=120) as response,p.open('wb') as output:
  while chunk:=response.read(1024**2):output.write(chunk)
h=hashlib.sha256()
with p.open('rb') as f:
 while chunk:=f.read(1024**2):h.update(chunk)
assert h.hexdigest()==expected,'Matrix checksum mismatch; do not use this file'
print('Matrix SHA256 verified',p.stat().st_size,flush=True)
