from pathlib import Path
import hashlib,urllib.request,shutil
R=Path(__file__).resolve().parent;dest=R/'GSE202642_matrix.mtx.gz';expected='855789bdd50f129a4131d84e285d0b52a1b481644accc747acc6f99dcc28e2e4'
if not dest.exists():
 temp=R/'GSE202642_matrix.mtx.gz.download'
 with urllib.request.urlopen('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_matrix.mtx.gz',timeout=90) as s,temp.open('wb') as o:shutil.copyfileobj(s,o,4*1024*1024)
 temp.replace(dest)
h=hashlib.sha256()
with dest.open('rb') as f:
 for c in iter(lambda:f.read(8388608),b''):h.update(c)
assert h.hexdigest()==expected,'Source hash mismatch; do not run clustering'
print('Official matrix SHA256 verified')
