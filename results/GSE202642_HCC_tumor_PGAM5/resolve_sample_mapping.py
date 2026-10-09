import gzip,re,json,urllib.request,time,xml.etree.ElementTree as ET,collections,csv
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
P=Path(__file__).resolve().parent
libs=collections.defaultdict(set)
for l in gzip.open(P/'GSE202642_barcodes.tsv.gz','rt'):
 b,x=l.strip().rsplit('-',1);libs[x].add(b)
soft=gzip.open(P/'GSE202642_family.soft.gz','rt').read()
records=[]
for block in soft.split('^SAMPLE = ')[1:]:
 gsm=block.splitlines()[0];title=re.search(r'!Sample_title = (.+)',block)[1];srx=re.search(r'SRA: .*?(SRX\d+)',block)[1];records.append((gsm,title,srx))
def fetch(url,path):
 if path.exists():return path.read_bytes()
 for n in range(4):
  try:
   d=urllib.request.urlopen(url,timeout=40).read();path.write_bytes(d);return d
  except Exception:
   if n==3:raise
   time.sleep(2+n)
def work(rec):
 gsm,title,srx=rec
 d=fetch('https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id='+srx+'&retmode=xml',P/(srx+'.xml'))
 runs=ET.fromstring(d).findall('.//RUN');assert len(runs)>=1
 eligible=[r for r in runs if any('_R1_' in x.attrib.get('filename','') for x in r.findall('.//SRAFile'))];assert eligible
 srr=eligible[0].attrib['accession']
 evidence=[]
 for start in [1,100001]:
  url=f'https://www.ncbi.nlm.nih.gov/Traces/sra-reads-be/read?retmode=json&acc={srr}&retstart={start}&retmax=1000'
  data=json.loads(fetch(url,P/f'{srr}_reads_{start}.json'));spots=data['spots'];assert len(spots)==1000
  bars=[]
  for s in spots:
   # Original 10x R1: the 28-bp technical read (16-bp barcode +12-bp UMI).
   inds=[i for i,n in enumerate(s['READ_LEN']) if n==28 and s['READ_TYPE'][i]==0];assert len(inds)==1
   off=s['READ_START'][inds[0]];bars.append(s['READ'][off:off+16])
  c={x:sum(b in bs for b in bars) for x,bs in libs.items()};best=max(c,key=c.get);runner=sorted(c.values())[-2]
  assert c[best]>=100 and c[best]>=10*max(runner,1),(gsm,c)
  evidence.append((best,c))
 assert evidence[0][0]==evidence[1][0]
 row={'GSM':gsm,'sample_name':title.rsplit('[',1)[1].rstrip(']'),'tissue':'HCC_tumor' if 'carcinoma' in title else 'adjacent_liver','SRX':srx,'SRR':srr,'library_suffix':evidence[0][0],'reads_checked':2000,'matches_batch1':evidence[0][1][evidence[0][0]],'matches_batch2':evidence[1][1][evidence[0][0]],'evidence_counts_json':json.dumps([e[1] for e in evidence],sort_keys=True)}
 print(row,flush=True);return row
with ThreadPoolExecutor(max_workers=2) as ex:rows=list(ex.map(work,records))
assert len(set(r['library_suffix'] for r in rows))==11
with open(P/'sample_mapping_verified.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
barcode_rows=[]
for row in rows:
 for start in [1,100001]:
  data=json.loads((P/f"{row['SRR']}_reads_{start}.json").read_text())
  for spot in data['spots']:
   i=next(i for i,n in enumerate(spot['READ_LEN']) if n==28 and spot['READ_TYPE'][i]==0)
   off=spot['READ_START'][i]
   barcode_rows.append({'GSM':row['GSM'],'SRR':row['SRR'],'batch_start':start,'spot_id':spot['SPOT_ID'],'raw_16bp_cell_barcode':spot['READ'][off:off+16]})
with gzip.open(P/'sample_mapping_read_barcode_evidence.csv.gz','wt') as f:
 w=csv.DictWriter(f,fieldnames=list(barcode_rows[0]));w.writeheader();w.writerows(barcode_rows)
print('All 11 suffixes uniquely resolved',flush=True)
