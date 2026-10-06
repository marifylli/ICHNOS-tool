import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]


def test_end_to_end_cli_report_and_no_overwrite(tmp_path):
    output=tmp_path/'verification'
    command=[sys.executable,str(ROOT/'scripts/verify_end_to_end.py'),'--out-dir',str(output)]
    result=subprocess.run(command,capture_output=True,text=True)
    assert result.returncode==0,result.stderr+result.stdout
    report=json.loads((output/'report.json').read_text())
    assert report['passed'] is True
    assert len(report['checks'])==8
    assert all(check['passed'] for check in report['checks'])
    assert report['experimentally_validated'] is False
    assert report['assumptions']['noise']=='none'
    assert report['assumptions']['known_session_scale']==1.5
    for name,digest in report['files_sha256'].items():
        assert hashlib.sha256((output/name).read_bytes()).hexdigest()==digest
    calibration=json.loads((output/'calibration.json').read_text())
    assert abs(calibration['c_session']-1.5)<1e-10
    joint=json.loads((output/'joint.json').read_text())
    assert joint['dose_estimate_uM']==75
    assert joint['elapsed_time_estimate_hours']==.5
    assert joint['sample_linkage']['calibration_applied_again'] is False
    ambiguous=json.loads((output/'ambiguous.json').read_text())
    assert len(ambiguous['candidate_pairs'])==6
    assert ambiguous['dose_estimate_uM'] is None
    assert json.loads((output/'low-green.json').read_text())['rejected'] is True
    before=(output/'report.json').read_bytes()
    again=subprocess.run(command,capture_output=True,text=True)
    assert again.returncode!=0
    assert (output/'report.json').read_bytes()==before


def test_visual_report_has_fixed_scales_counts_and_integrity_guard(tmp_path):
    import runpy
    import numpy as np
    from PIL import Image
    output=tmp_path/'run'
    run=subprocess.run([sys.executable,str(ROOT/'scripts/verify_end_to_end.py'),'--out-dir',str(output)],capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    render=runpy.run_path(str(ROOT/'scripts/render_e2e_report.py'))['render']
    visual=tmp_path/'visual'
    index=render(output,visual)
    html=index.read_text()
    assert 'no experimental validation' in html
    assert 'time-course.png' in html and 'dose-response.png' in html
    assert 'known synthetic generation doses' in html
    for chart in ('time-course.png','dose-response.png'):
        assert Image.open(visual/chart).width > 1000
    assert 'insufficient usable cells' in html
    metadata=json.loads((visual/'preview_metadata.json').read_text())
    assert len(metadata['images'])==12
    assert all(row['extracted_cells']==3 for row in metadata['images'])
    assert metadata['display_maxima']['green']==1000
    strong=np.asarray(Image.open(visual/'dose75-0-green.png'))
    weak=np.asarray(Image.open(visual/'low-green-0-green.png'))
    assert strong[:,:,1].max()==255
    assert 20<weak[:,:,1].max()<30
    import pytest
    with pytest.raises(FileExistsError):
        render(output,visual)
    cells=output/'dose75-cells.csv'
    cells.write_text(cells.read_text()+'\n')
    with pytest.raises(ValueError,match='differs from verification report'):
        render(output,tmp_path/'changed')
    assert not (tmp_path/'changed').exists()
