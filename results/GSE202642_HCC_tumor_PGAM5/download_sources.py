from pathlib import Path
import urllib.request,shutil
R=Path(__file__).resolve().parent
sources={
 'GSE202642_matrix.mtx.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_matrix.mtx.gz',
 'GSE202642_features.tsv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_features.tsv.gz',
 'GSE202642_barcodes.tsv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/suppl/GSE202642_barcodes.tsv.gz',
 'GSE202642_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202642/soft/GSE202642_family.soft.gz'}
for name,url in sources.items():
 dest=R/name
 if dest.exists():
  print('Existing',name,'verify SHA256 against input_sha256.json');continue
 temp=dest.with_suffix(dest.suffix+'.download')
 with urllib.request.urlopen(url,timeout=90) as src,temp.open('wb') as out:shutil.copyfileobj(src,out,4*1024*1024)
 temp.replace(dest);print('Downloaded',name,flush=True)
