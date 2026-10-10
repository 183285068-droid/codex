"""Check downloaded payloads before rerunning analyses that overwrite outputs."""
from pathlib import Path
import hashlib
import pandas as pd

R=Path(__file__).resolve().parent

def main():
    table=pd.read_csv(R/'SHA256_manifest.csv')
    for _,row in table.iterrows():
        p=R/row['path'];assert p.is_file(),row['path']
        assert p.stat().st_size==row['bytes'],row['path']
        assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
    print('PASS:',len(table),'payload files verified.')

if __name__=='__main__':main()
