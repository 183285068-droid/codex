from pathlib import Path
import urllib.request,json,hashlib
R=Path(__file__).resolve().parent;I=R/'inputs';I.mkdir(exist_ok=True)
for name,expected in json.loads((R/'source_summary.json').read_text())['input_hashes'].items():
 path=I/name
 if not path.exists():
  url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/'+name
  print('Downloading',name,flush=True);urllib.request.urlretrieve(url,path)
 with open(path,'rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==expected
print('Official inputs SHA256 verified')
