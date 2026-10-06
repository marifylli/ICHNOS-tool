"""Synthetic image check: two observables recover a late oxidative snapshot."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from ichnos.calibration import fit_session_calibration, model_reference_id
from ichnos.snapshot_decoder import build_observable_table, build_snapshot_calibration, decode_snapshot, save_artifact
from ichnos_image.pipeline import ImageSet, process_image_set
from ichnos_image.export import export_csv

ROOT=Path(__file__).resolve().parents[1]
ACQUISITION=dict(objective='40X',green_units='corrected_synthetic_image_units',
    extraction_green='synthetic scalar plane',extraction_red='synthetic scalar plane',
    gain_setting='fixed synthetic gain',exposure_ms_green=100,exposure_ms_red=100)


def image_summary(out,name,model_green,model_ratio):
    yy,xx=np.mgrid[:128,:128]
    mask=np.zeros((128,128),dtype=bool)
    for cy,cx in ((32,32),(64,96),(96,32)):
        mask |= (yy-cy)**2+(xx-cx)**2<100
    green=np.where(mask,20.+100.*model_green,0.)
    red=1.5*model_ratio*green
    bright=np.where(mask,.2,.5)
    np.savez_compressed(out/f'{name}.npz',green=green,red=red,bright_field=bright)
    records=process_image_set(ImageSet(green=green,red=red,bright_field=bright,
        session_id='snapshot-synthetic',sample_id=name,condition_id=name,timepoint=0,acquisition_order=0,
        exposure_ms_green=100,exposure_ms_red=100,nd_filter_green=0,nd_filter_red=0,
        objective='40X',burner_hours=1,lamp_warmup_minutes=30),bleed_green_to_red=0)
    export_csv(records,out/f'{name}-cells.csv')
    if len(records)!=3 or not all(row.qc_pass for row in records):
        raise ValueError('synthetic scene did not produce three QC-passing cells')
    return float(np.median([row.corrected_mean_green for row in records])),float(np.median([row.ratio_red_green for row in records]))


def run(out):
    out=Path(out)
    if out.exists():
        raise FileExistsError(f'output directory already exists: {out}')
    out.mkdir(parents=True)
    reference_dir=out/'ratio-reference'
    result=subprocess.run([sys.executable,str(ROOT/'scripts/run_protocol.py'),
        '--variant','ox','--dose','25','--dose-units','uM','--times-hours','.5','1','8',
        '--initialization','equilibrium','--out-dir',str(reference_dir)],capture_output=True,text=True)
    if result.returncode:
        raise ValueError(result.stderr)
    reference=json.loads((reference_dir/'metadata.json').read_text())
    fluorescence=pd.read_csv(reference_dir/'fluorescence.csv')
    measured=[image_summary(out,f'ratio-reference-{i}',row.green,row.ratio_red_green)[1]
              for i,row in enumerate(fluorescence.itertuples())]
    ratio=fit_session_calibration(fluorescence.ratio_red_green.to_numpy(),measured,
        session_id='snapshot-synthetic',reference_id=model_reference_id(reference),source_kind='synthetic',
        source='separate model-generated reference images at 25 uM; not independent biological validation')
    table=build_observable_table(variant='ox',doses_uM=[0,10,200,400,800],times_hours=[.5,1,8],initialization='equilibrium')
    save_artifact(table.to_dict(),out/'table.json')
    low,_=image_summary(out,'green-low',table.green[0,2],table.ratio_table.ratios[0,2])
    high,_=image_summary(out,'green-high',table.green[-1,2],table.ratio_table.ratios[-1,2])
    config=dict(session_id='snapshot-synthetic',source_kind='synthetic',acquisition=ACQUISITION,
        low_reference=dict(specimen_id='low',biological_replicate_id='synthetic-low-culture',
            source='synthetic unstressed reference image',dose_uM=0,time_hours=8,corrected_green=low),
        high_reference=dict(specimen_id='high',biological_replicate_id='synthetic-high-culture',
            source='synthetic high-dose reference image; saturation not established',dose_uM=800,time_hours=8,corrected_green=high),
        min_image_span=100.,min_model_span=1.,ratio_calibration=ratio.to_dict(),ratio_reference_metadata=reference,
        ratio_reference_specimen_ids=[f'ratio-reference-{i}' for i in range(3)],
        ratio_reference_replicate_ids=['synthetic-ratio-culture'])
    save_artifact(config,out/'calibration-config.json')
    calibration=build_snapshot_calibration(table,**config)
    save_artifact(calibration,out/'calibration.json')
    truth=build_observable_table(variant='ox',doses_uM=[400,801],times_hours=[8],initialization='equilibrium')
    target_green,target_ratio=image_summary(out,'target',truth.green[0,0],truth.ratio_table.ratios[0,0])
    observation=dict(session_id='snapshot-synthetic',specimen_id='target',biological_replicate_id='synthetic-target-culture',
        data_kind='synthetic',acquisition=ACQUISITION,ratio_scale='image',
        image_ratio_red_green=target_ratio,corrected_green=target_green)
    save_artifact(observation,out/'observation.json')
    result=decode_snapshot(table,calibration,observation,ratio_tolerance=.005,green_tolerance=.002,green_floor=100.)
    save_artifact(result,out/'decoded.json')
    passed=(result['status']=='unique_grid_pair' and result['dose_estimate_uM']==400 and
            result['elapsed_time_estimate_hours']==8 and result['ratio_only_candidate_count']>1)
    report=dict(passed=passed,expected=dict(dose_uM=400,elapsed_time_hours=8),
        actual=dict(dose_uM=result['dose_estimate_uM'],elapsed_time_hours=result['elapsed_time_estimate_hours']),
        ratio_only_candidate_count=result['ratio_only_candidate_count'],joint_candidate_count=len(result['candidate_pairs']),
        known_green_equation='corrected_image_green = 20 + 100 * Observed_Green',
        known_ratio_equation='image_ratio = 1.5 * Measured_Ratio_RG',
        green_reference_times_hours=[8,8],green_reference_doses_uM=[0,800],
        ratio_tolerance=.005,green_norm_tolerance=.002,noise='none',
        biological_replicate_ids_are_synthetic_labels=True,experimentally_validated=False,posterior_computed=False,
        limitations='same-model software verification; no biological noise, independent experimental calibration, or domain-wide identifiability claim')
    save_artifact(report,out/'verification.json')
    print(json.dumps(report,indent=2))
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir',type=Path,required=True)
    args=parser.parse_args()
    try:
        report=run(args.out_dir)
    except (ValueError,KeyError,OSError) as exc:
        parser.exit(1,f'Error: {exc}\n')
    if not report['passed']:
        parser.exit(1,'Snapshot recovery check failed\n')


if __name__=='__main__':
    main()
