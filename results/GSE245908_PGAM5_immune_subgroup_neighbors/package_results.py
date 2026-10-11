"""Package reviewed results with a file manifest and verify ZIP integrity."""
from pathlib import Path
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parent
NAME = 'GSE245908_PGAM5_immune_subgroup_neighbors_results.zip'


def main():
    excluded = {'__pycache__', 'work', '.venv', '.git'}
    files = sorted(p for p in ROOT.rglob('*') if p.is_file()
                   and not excluded.intersection(p.relative_to(ROOT).parts)
                   and not any(part.startswith('.') for part in p.relative_to(ROOT).parts)
                   and p.name not in {NAME, 'package_manifest.json'})
    records = [{'path': str(p.relative_to(ROOT)), 'bytes': p.stat().st_size,
                'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
    archive = ROOT / NAME
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in files:
            z.write(p, str(p.relative_to(ROOT)))
        z.writestr('bundle_file_manifest.json', json.dumps(records, indent=2) + '\n')
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for record in records:
            data = z.read(record['path'])
            assert len(data) == record['bytes']
            assert hashlib.sha256(data).hexdigest() == record['sha256']
    assert archive.stat().st_size < 30 * 1024**2
    result = {'archive': NAME, 'bytes': archive.stat().st_size,
              'MiB': archive.stat().st_size / 1024**2,
              'sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
              'ZIP_CRC_verified': True, 'files_verified': len(records), 'files': records}
    (ROOT / 'package_manifest.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'files'}, indent=2))


if __name__ == '__main__':
    main()
