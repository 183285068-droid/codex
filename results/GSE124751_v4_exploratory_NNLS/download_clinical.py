"""Download public official clinical payloads; no credentials needed.

Run after media.githubusercontent.com is allowed. Original LFS pointers are
preserved separately. These files must be inspected before endpoint analysis.
"""
from pathlib import Path
import hashlib
import urllib.request
import urllib.error
import sys

ROOT = Path(__file__).resolve().parent
COMMIT = '18763f8eee86d4886d7210068f745f1e509498ee'
FILES = {
    'patient': ('629c2d2994250b55f069ddd4a7a5900f5ac0f3ab1d9da37b96686f2c8e8729c1',12050),
    'sample': ('5ecf933d8834175841f59efa2f8b9f251531877fd89a58a1da81b8042762ac49',20728),
}


def main():
    for kind,(sha,size) in FILES.items():
        name = f'data_clinical_{kind}.txt'
        url = f'https://media.githubusercontent.com/media/cBioPortal/datahub/{COMMIT}/public/lihc_amc_prv/{name}'
        output = ROOT/'inputs'/f'lihc_amc_prv_{name}'
        if output.exists() and hashlib.sha256(output.read_bytes()).hexdigest()==sha:
            print('Verified cached clinical file:',output.name)
            continue
        try:
            with urllib.request.urlopen(url,timeout=30) as response:
                data = response.read()
        except urllib.error.URLError as error:
            print('Public clinical download failed:', str(error.reason), file=sys.stderr)
            print('Check network access to media.githubusercontent.com; no survival calculation was performed.',file=sys.stderr)
            raise SystemExit(2)
        if len(data)!=size or hashlib.sha256(data).hexdigest()!=sha:
            raise RuntimeError('Clinical payload does not match pinned official Git LFS object: '+name)
        output.write_bytes(data)
        print('Downloaded and SHA256 verified:',output.name)
    print('Next: inspect endpoint definitions, exact patient ID linkage, age agreement and missingness before survival analysis.')


if __name__ == '__main__':
    main()
