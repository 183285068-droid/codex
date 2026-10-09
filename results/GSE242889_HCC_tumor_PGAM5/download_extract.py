from pathlib import Path
import urllib.request,tarfile
R=Path(__file__).resolve().parent
urls=['https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242889/soft/GSE242889_family.soft.gz']
for gsm,name in [(7774399,'1T_C21'),(7774400,'2T_C24'),(7774401,'3T_C25'),(7774402,'4T_C29'),(7774403,'5T_C36')]:
 urls.append(f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7774nnn/GSM{gsm}/suppl/GSM{gsm}_{name}.tar.gz')
for url in urls:
 p=R/url.rsplit('/',1)[1]
 if not p.exists():urllib.request.urlretrieve(url,p)
 if p.name.startswith('GSM'):
  with tarfile.open(p,'r:') as t:
   assert all(not m.issym() and not m.islnk() and not Path(m.name).is_absolute() and '..' not in Path(m.name).parts for m in t.getmembers())
   t.extractall(R/'counts',filter='data')
 print(p.name)
