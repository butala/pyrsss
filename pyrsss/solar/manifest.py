"""SHA256SUMS manifests for fetched data, the way SphericalCT pins fixtures.

A fetched archive is only as trustworthy as its record of what was fetched:
this module writes and checks ``SHA256SUMS`` files in the style of
``shasum -a 256 -c`` -- one ``<digest>  <name>`` line per file -- so a
campaign's data can be re-verified later and a re-fetch can be told apart
from a silently different product.

Pure standard library: no dependency on the fetchers themselves.
"""

import hashlib
from pathlib import Path

MANIFEST_NAME = 'SHA256SUMS'
_CHUNK = 1 << 20


def sha256(path):
    """
    Return the hex SHA-256 digest of *path*, streamed so a 100 MB FITS file
    does not land in memory.
    """
    h = hashlib.sha256()
    with open(path, 'rb') as fid:
        while chunk := fid.read(_CHUNK):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(directory, paths=None):
    """
    Write ``SHA256SUMS`` in *directory* covering *paths* (default: every
    regular file there except the manifest itself). Returns the manifest
    path.
    """
    directory = Path(directory)
    if paths is None:
        paths = sorted(p for p in directory.iterdir()
                       if p.is_file() and p.name != MANIFEST_NAME)
    manifest = directory / MANIFEST_NAME
    with open(manifest, 'w') as fid:
        for path in paths:
            fid.write(f'{sha256(path)}  {Path(path).name}\n')
    return manifest


def read_manifest(manifest):
    """
    Parse a ``SHA256SUMS`` file into ``{name: digest}``.
    """
    out = {}
    with open(manifest) as fid:
        for line in fid:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            digest, _, name = line.partition('  ')
            out[Path(name).name] = digest
    return out


def verify(directory, manifest=None):
    """
    Check every file listed in the manifest of *directory*.

    Returns ``(ok, missing, changed)`` -- the names that verify, the names
    listed but absent, and the names present but digested differently. A
    caller that wants a hard failure raises on non-empty ``missing`` or
    ``changed``; this function is quiet so a survey can report.
    """
    directory = Path(directory)
    manifest = Path(manifest) if manifest is not None \
        else directory / MANIFEST_NAME
    expected = read_manifest(manifest)
    ok, missing, changed = [], [], []
    for name, digest in sorted(expected.items()):
        path = directory / name
        if not path.is_file():
            missing.append(name)
        elif sha256(path) == digest:
            ok.append(name)
        else:
            changed.append(name)
    return ok, missing, changed
