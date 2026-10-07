import json
import pytest
from ichnos.artifacts import publish_csv_bundle, sha256_file
from ichnos_image.pipeline import process_experiment


def test_existing_results_survive_even_when_inputs_are_invalid(tmp_path):
    path = tmp_path / 'cells.csv'; path.write_text('valuable data')
    with pytest.raises((FileExistsError, AttributeError)):
        process_experiment([None], path, bleed_green_to_red=0)
    assert path.read_text() == 'valuable data'


def test_failure_does_not_publish_partial_csv_or_manifest(tmp_path):
    def fail(path):
        path.write_text('partial')
        raise RuntimeError('second image failed')
    with pytest.raises(RuntimeError):
        publish_csv_bundle(tmp_path / 'cells.csv', fail)
    assert list(tmp_path.iterdir()) == []


def test_manifest_hash_and_no_overwrite(tmp_path):
    def writer(path):
        path.write_text('complete')
        return {'schema_version': 2}
    path = tmp_path / 'cells.csv'
    publish_csv_bundle(path, writer)
    manifest = json.loads((tmp_path / 'cells.csv.manifest.json').read_text())
    assert manifest['output_csv_sha256'] == sha256_file(path)
    with pytest.raises(FileExistsError):
        publish_csv_bundle(path, writer)
    assert path.read_text() == 'complete'


def test_direct_export_also_refuses_overwrite(tmp_path):
    from ichnos_image.export import export_csv
    path = tmp_path/'cells.csv'
    export_csv([],path)
    original = path.read_bytes()
    with pytest.raises(FileExistsError):
        export_csv([],path)
    assert path.read_bytes() == original


def test_manifest_publication_failure_rolls_back_owned_csv(monkeypatch,tmp_path):
    import ichnos.artifacts as artifacts
    real_link = artifacts.os.link
    def link(source,target):
        if str(target).endswith('.manifest.json'):
            raise OSError('simulated manifest publication failure')
        return real_link(source,target)
    monkeypatch.setattr(artifacts.os,'link',link)
    def writer(path):
        path.write_text('complete')
        return {}
    with pytest.raises(OSError):
        publish_csv_bundle(tmp_path/'cells.csv',writer)
    assert list(tmp_path.iterdir()) == []
