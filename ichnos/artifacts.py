"""No-overwrite publication, input fingerprints and source provenance."""
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib
import importlib.metadata
import json
import os
import subprocess


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def file_fingerprints(paths):
    return {str(Path(p).resolve()): sha256_file(p) for p in sorted(set(map(str, paths)))}


def code_provenance():
    root = Path(__file__).resolve().parents[1]
    files = sorted(p for folder in ('ichnos', 'ichnos_image', 'scripts')
                   for p in (root / folder).rglob('*.py'))
    hashes = {str(p.relative_to(root)): sha256_file(p) for p in files}
    def git(*args):
        try:
            result = subprocess.run(['git', '-C', str(root), *args], capture_output=True, text=True)
        except FileNotFoundError:
            return None
        return result.stdout.strip() if result.returncode == 0 else None
    try:
        version = importlib.metadata.version('ichnos-tool')
    except importlib.metadata.PackageNotFoundError:
        version = 'uninstalled'
    runtime = {}
    for package in ('numpy', 'scipy', 'pandas', 'scikit-image', 'pillow', 'libroadrunner', 'tellurium'):
        try:
            runtime[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            runtime[package] = None
    return dict(package_version=version, runtime_packages=runtime, git_commit=git('rev-parse', 'HEAD'),
                git_status=git('status', '--porcelain'), source_sha256=hashes,
                source_tree_sha256=hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest())


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


@contextmanager
def output_lock(destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock = destination.with_name(destination.name + '.lock')
    with lock.open('x'):
        try:
            if destination.exists():
                raise FileExistsError(f'output path already exists: {destination}')
            yield
        finally:
            lock.unlink()


def publish_csv_bundle(destination, writer):
    """Publish CSV + manifest, never replacing either; manifest is commit marker.

    A process crash between links can leave a CSV without a manifest; it is
    incomplete, never an accepted run. Ordinary exceptions remove owned outputs.
    """
    destination = Path(destination)
    manifest = destination.with_suffix(destination.suffix + '.manifest.json')
    with output_lock(destination):
        if manifest.exists():
            raise FileExistsError(f'output manifest already exists: {manifest}')
        with TemporaryDirectory(prefix='.ichnos-', dir=destination.parent) as directory:
            csv = Path(directory) / 'cells.csv'
            meta = Path(directory) / 'manifest.json'
            metadata = writer(csv)
            metadata['output_csv_sha256'] = sha256_file(csv)
            write_json(meta, metadata)
            os.link(csv, destination)  # atomic no-clobber, even with unrelated writers
            try:
                os.link(meta, manifest)
            except BaseException:
                destination.unlink()
                raise
    return destination


@contextmanager
def new_output_directory(destination):
    """Build a complete bundle off-path; publish after successful completion."""
    destination = Path(destination)
    with output_lock(destination):
        with TemporaryDirectory(prefix='.ichnos-', dir=destination.parent) as temporary:
            staged = Path(temporary) / 'run'
            staged.mkdir()
            yield staged
            if destination.exists():
                raise FileExistsError(f'output path already exists: {destination}')
            staged.rename(destination)
