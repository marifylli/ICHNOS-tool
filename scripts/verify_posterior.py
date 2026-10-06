"""Synthetic replicate-noise and posterior check using a snapshot verification run."""
import argparse
import json
from pathlib import Path

import numpy as np
from ichnos.posterior import fit_replicate_noise, decode_posterior, prior_grid
from ichnos.snapshot_decoder import load_observable_table, save_artifact, validate_snapshot_calibration


def run(snapshot_dir,out_dir,*,seed=20261006,heldout_count=24):
    snapshot_dir,out_dir=Path(snapshot_dir),Path(out_dir)
    if out_dir.exists():
        raise FileExistsError(f'output directory already exists: {out_dir}')
    if isinstance(heldout_count,bool) or heldout_count<1:
        raise ValueError('heldout_count must be positive')
    table=load_observable_table(snapshot_dir/'table.json')
    calibration=json.loads((snapshot_dir/'calibration.json').read_text())
    validate_snapshot_calibration(table,calibration)
    if calibration['source_kind']!='synthetic':
        raise ValueError('this synthetic verifier requires synthetic calibration artifacts')
    generator=np.random.default_rng(seed)
    known_covariance=np.array([[.02**2,-.35*.02*.035],[-.35*.02*.035,.035**2]])
    ilow,ihigh=calibration['image_green_anchors']
    mlow,mhigh=calibration['model_green_anchors']
    means=np.stack([np.log(table.ratio_table.ratios),(table.green-mlow)/(mhigh-mlow)],axis=-1)
    def generate(i,j,name):
        z=generator.multivariate_normal(means[i,j],known_covariance)
        return dict(session_id=calibration['session_id'],specimen_id=name,biological_replicate_id=name+'-culture',
            data_kind='synthetic',ratio_scale='image',acquisition=calibration['acquisition'],
            observation_unit='biological_replicate_summary',
            image_ratio_red_green=float(np.exp(z[0])*calibration['ratio_calibration']['c_session']),
            corrected_green=float(ilow+z[1]*(ihigh-ilow)),
            source='independent synthetic biological-summary label; not an experimental replicate')
    training=[]
    # Three distinct conditions, twelve separately drawn synthetic cultures each.
    for group,(i,j) in enumerate([(0,0),(len(table.green)//2,1),(len(table.green)-1,len(table.green[0])-1)]):
        for k in range(12):
            row=generate(i,j,f'train-{group}-{k}')
            row.update(condition_id=f'condition-{group}',dose_uM=float(table.ratio_table.doses_uM[i]),
                       time_hours=float(table.ratio_table.times_hours[j]))
            training.append(row)
    noise=fit_replicate_noise(table,calibration,training,green_floor=100,min_replicates_per_condition=3)
    prior={'kind':'uniform_grid'}
    i=np.flatnonzero(np.isclose(table.ratio_table.doses_uM,400,rtol=0,atol=1e-12))
    j=np.flatnonzero(np.isclose(table.ratio_table.times_hours,8,rtol=0,atol=1e-12))
    if len(i)!=1 or len(j)!=1:
        raise ValueError('example verification requires the 400 uM / 8 h model point')
    target=generate(int(i[0]),int(j[0]),'heldout-illustration')
    result=decode_posterior(table,calibration,noise,target,prior=prior,max_mahalanobis_squared=25,credible_mass=.95)
    probabilities=prior_grid(table,prior).ravel()
    experiments=[]
    heldout_rows=[]
    hits=0
    for k in range(heldout_count):
        index=generator.choice(len(probabilities),p=probabilities)
        i,j=np.unravel_index(index,table.green.shape)
        row=generate(int(i),int(j),f'heldout-{k}')
        inferred=decode_posterior(table,calibration,noise,row,prior=prior,max_mahalanobis_squared=25,credible_mass=.95)
        truth=dict(dose_uM=float(table.ratio_table.doses_uM[i]),elapsed_time_hours=float(table.ratio_table.times_hours[j]))
        included=any(candidate['dose_uM']==truth['dose_uM'] and candidate['elapsed_time_hours']==truth['elapsed_time_hours']
                     for candidate in inferred['credible_set'])
        hits+=int(included)
        experiments.append(dict(**truth,in_credible_set=included,status=inferred['status'],
                                credible_mass_achieved=inferred.get('credible_mass_achieved')))
        heldout_rows.append(row)
    passed=(np.isclose(np.sum(result['posterior_mass']),1) and result['credible_mass_achieved']>=.95)
    report=dict(computational_checks_passed=bool(passed),source_kind='synthetic',seed=seed,
        n_training_cultures=len(training),n_conditions=3,known_covariance=known_covariance.tolist(),
        estimated_covariance=noise['covariance'],illustration_true_pair=dict(dose_uM=400,elapsed_time_hours=8),
        illustration_status=result['status'],illustration_map_grid_pairs=result['map_grid_pairs'],
        credible_mass_requested=.95,n_heldout=heldout_count,n_true_pairs_in_credible_set=hits,
        heldout_inclusion_fraction=hits/heldout_count,
        heldout_prior=prior,heldout_results=experiments,
        calibration_uncertainty_propagated=False,covariance_uncertainty_propagated=False,
        experimentally_validated=False,
        limitations='small synthetic inclusion experiment under the generating model/prior; not experimental coverage validation; plug-in covariance and fixed calibration')
    out_dir.mkdir(parents=True)
    for name,value in [('table',table.to_dict()),('calibration',calibration),('replicates',training),('noise',noise),
                       ('prior',prior),('observation',target),('posterior',result),('heldout-observations',heldout_rows),('verification',report)]:
        save_artifact(value,out_dir/f'{name}.json')
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    fig,axis=plt.subplots(figsize=(8,5),constrained_layout=True)
    raster=axis.imshow(result['posterior_mass'],aspect='auto',origin='lower',vmin=0,cmap='viridis')
    axis.set_xticks(range(len(table.ratio_table.times_hours)),[f'{x:g}' for x in table.ratio_table.times_hours])
    axis.set_yticks(range(len(table.ratio_table.doses_uM)),[f'{x:g}' for x in table.ratio_table.doses_uM])
    for point in result['credible_set']:
        col=int(np.flatnonzero(table.ratio_table.times_hours==point['elapsed_time_hours'])[0])
        row=int(np.flatnonzero(table.ratio_table.doses_uM==point['dose_uM'])[0])
        axis.plot(col,row,'s',markerfacecolor='none',markeredgecolor='white',markersize=22)
    axis.set(xlabel='Elapsed time grid node (hours)',ylabel='Dose grid node (μM)',
             title='Conditional synthetic posterior mass\nWhite squares: 95% discrete credible set')
    fig.colorbar(raster,ax=axis,label='Probability mass per grid node')
    fig.savefig(out_dir/'posterior.png',dpi=160)
    plt.close(fig)
    print(json.dumps({k:v for k,v in report.items() if k!='heldout_results'},indent=2))
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot-dir',type=Path,required=True)
    parser.add_argument('--out-dir',type=Path,required=True)
    parser.add_argument('--seed',type=int,default=20261006)
    parser.add_argument('--heldout-count',type=int,default=24)
    args=parser.parse_args()
    try:
        report=run(args.snapshot_dir,args.out_dir,seed=args.seed,heldout_count=args.heldout_count)
    except (ValueError,KeyError,OSError) as exc:
        parser.exit(1,f'Error: {exc}\n')
    if not report['computational_checks_passed']:
        parser.exit(1,'Posterior normalization/credible-mass check failed\n')


if __name__=='__main__':
    main()
