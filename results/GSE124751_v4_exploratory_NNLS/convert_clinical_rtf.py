"""Decode these two ASCII macOS RTF tables; require exact official SHA256.

This is a constrained reader for the supplied documents, not a general RTF
parser. Reject unexpected controls; never execute document content.
"""
from pathlib import Path
import hashlib, json

ROOT = Path(__file__).resolve().parent
EXPECTED = {
    'patient':('629c2d2994250b55f069ddd4a7a5900f5ac0f3ab1d9da37b96686f2c8e8729c1',12050),
    'sample':('5ecf933d8834175841f59efa2f8b9f251531877fd89a58a1da81b8042762ac49',20728),
}


def main():
    audit=[]
    for kind,(expected,size) in EXPECTED.items():
        source=ROOT/'inputs'/f'uploaded_data_clinical_{kind}.rtf'
        raw=source.read_bytes()
        text=raw.decode('ascii')
        assert text.startswith('{\\rtf1') and text.endswith('}')
        body=text[text.index('#Patient Identifier'):-1]
        body=body.replace('\\\n','\n')
        assert not any(c in body for c in ['\\','{','}']), 'Unexpected RTF body controls'
        decoded=body.encode('ascii')
        actual=hashlib.sha256(decoded).hexdigest()
        assert actual==expected and len(decoded)==size, 'Not the pinned official clinical payload'
        output=ROOT/'inputs'/f'lihc_amc_prv_data_clinical_{kind}.txt'
        output.write_bytes(decoded)
        audit.append({'kind':kind,'uploaded_RTF_SHA256':hashlib.sha256(raw).hexdigest(),
                      'decoded_SHA256':actual,'decoded_bytes':len(decoded),
                      'exact_official_LFS_payload_match':True})
    (ROOT/'clinical_RTF_conversion_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('Both decoded clinical tables exactly match official pinned Git LFS payloads.')


if __name__=='__main__':
    main()
