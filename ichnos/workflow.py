"""One configuration -> image QC, paired summaries and conditional inference."""
from collections import Counter
from dataclasses import asdict
from html import escape
import json
from pathlib import Path

import pandas as pd

from .artifacts import code_provenance, file_fingerprints, new_output_directory, sha256_file, write_json
from .population import summarize_cells
from .snapshot_adapter import decode_summary
from .snapshot_decoder import load_observable_table, validate_snapshot_calibration
from .posterior import validate_noise, prior_grid


def run_workflow(config_path, out_dir):
    from ichnos_image.manifest import _build_image_sets, _calibrate_bleed_per_session
    from ichnos_image.pipeline import process_experiment
    from ichnos_image.focus import FocusPolicy
    from ichnos_image.saturation_audit import saturation_impact
    config_path, out_dir = Path(config_path).resolve(), Path(out_dir)
    if out_dir.exists():
        raise FileExistsError(f'output path already exists: {out_dir}')
    config = json.loads(config_path.read_text())
    allowed = {'image_manifest', 'controls', 'bleed_default', 'data_kind', 'image_options',
               'table', 'calibration', 'noise', 'prior', 'max_mahalanobis_squared',
               'credible_mass', 'min_cells', 'green_floor', 'ratio_tolerance', 'green_tolerance'}
    if set(config) - allowed:
        raise ValueError(f'unknown configuration fields: {sorted(set(config)-allowed)}')
    if config.get('data_kind') not in {'synthetic', 'experimental'}:
        raise ValueError('data_kind must explicitly be synthetic or experimental')
    def path(key):
        return (config_path.parent / config[key]).resolve()
    def read(key):
        return json.loads(path(key).read_text())
    table, calibration = load_observable_table(path('table')), read('calibration')
    validate_snapshot_calibration(table, calibration)
    noise = read('noise') if config.get('noise') else None
    prior = read('prior') if noise is not None else None
    if noise is not None:
        validate_noise(table, calibration, noise)
        prior_grid(table, prior)
        if 'max_mahalanobis_squared' not in config:
            raise ValueError('max_mahalanobis_squared is required for posterior inference')
    elif not {'ratio_tolerance','green_tolerance','green_floor'}.issubset(config):
        raise ValueError('snapshot compatibility requires tolerances and green_floor')
    floor = config.get('green_floor', noise['green_floor'] if noise is not None else None)
    if noise is not None and floor != noise['green_floor']:
        raise ValueError('green_floor differs from noise training')
    image_options = dict(config.get('image_options', {}))
    policy = FocusPolicy(**image_options.pop('focus_policy', {}))
    if config['data_kind'] == 'experimental':
        if policy.mode != 'enforce':
            raise ValueError('experimental workflow requires explicit enforced focus thresholds')
        if not config.get('controls') or config.get('bleed_default') is not None:
            raise ValueError('experimental workflow requires controls and forbids bleed_default')
    extraction = {k:image_options.pop(k) for k in ('green_extraction','red_extraction','saturation_value') if k in image_options}
    images = _build_image_sets(path('image_manifest'), **extraction)
    if not images:
        raise ValueError('empty image manifest')
    for item in images:
        if not item.sample_id:
            raise ValueError('every image requires an explicit sample_id')
    controls = path('controls') if config.get('controls') else None
    bleed = _calibrate_bleed_per_session(controls, {im.session_id for im in images},
        config.get('bleed_default'), **{k:v for k,v in extraction.items() if k != 'saturation_value'})
    sources = [config_path, path('image_manifest'), path('table'), path('calibration')]
    if noise is not None:
        sources += [path('noise'), path('prior')]
    if controls is not None:
        sources += [controls] + [controls.parent / p for p in pd.read_csv(controls)[
            ['control_green_path','control_red_path']].to_numpy().ravel()]
    sources += [p for item in images for p in item.source_paths]
    fingerprints = file_fingerprints(sources)
    with new_output_directory(out_dir) as out:
        cells_path = process_experiment(images, out/'cells.csv', input_paths=sources,
            bleed_green_to_red=bleed, focus_policy=policy, **image_options)
        cells = pd.read_csv(cells_path, dtype={k:str for k in
            ('session_id','sample_id','condition_id','specimen_id','biological_replicate_id','acquisition_json')})
        summaries = summarize_cells(cells, data_kind=config['data_kind'],
            min_cells=config.get('min_cells',3), green_floor=floor)
        summaries.to_csv(out/'samples.csv', index=False)
        write_json(out/'saturation-impact.json', saturation_impact(cells))
        write_json(out/'table.json', table.to_dict())
        write_json(out/'calibration.json', calibration)
        if noise is not None:
            write_json(out/'noise.json', noise); write_json(out/'prior.json', prior)
        results, seen = [], set()
        for index, row in summaries.iterrows():
            summary = row.to_dict()
            key = (summary['session_id'], summary['sample_id'], int(summary['timepoint']))
            seen.add(key)
            identity = dict(session_id=key[0], sample_id=key[1], timepoint=key[2])
            try:
                result = decode_summary(table,calibration,summary,
                    source=f"cells.csv sha256:{sha256_file(cells_path)}", noise=noise,prior=prior,
                    max_mahalanobis_squared=config.get('max_mahalanobis_squared'),
                    credible_mass=config.get('credible_mass',.95),green_floor=floor,
                    ratio_tolerance=config.get('ratio_tolerance'),green_tolerance=config.get('green_tolerance'))
            except (ValueError, KeyError, TypeError) as exc:
                result = dict(status='rejected',reason=str(exc),posterior_computed=False)
            filename = f'decoded-{index:04d}.json'
            write_json(out/filename, result)
            status = result['status']
            results.append({**identity, 'status':status, 'result_file':filename,
                'reason':result.get('reason', status if status in {'poor_model_fit','prior_only','ambiguous',
                    'below_detection_floor','out_of_response_domain','no_grid_pair_match'} else ''),
                'credible_set_size':len(result.get('credible_set',[])),
                'map_grid_pairs':result.get('map_grid_pairs',[])})
        for item in images:
            key = (item.session_id,item.sample_id,item.timepoint)
            if key not in seen:
                seen.add(key)
                results.append(dict(session_id=key[0],sample_id=key[1],timepoint=key[2],
                    status='rejected',reason='no_segmented_cells',result_file=None,credible_set_size=0,map_grid_pairs=[]))
        counts = Counter(reason for value in cells.qc_reasons.fillna('') for reason in value.split(';') if reason)
        report = dict(schema_version=1,data_kind=config['data_kind'],experimentally_validated=False,
            n_cells=len(cells),n_cells_qc_pass=int(cells.qc_pass.sum()),qc_rejection_counts=dict(counts),
            focus_policy=asdict(policy),focus_status_counts=cells.focus_status.value_counts().to_dict(),
            samples=results,status_counts=dict(Counter(r['status'] for r in results)),
            limitations=['conditional on model, calibration, grid, prior and noise',
                         'credible sets may be disconnected; a MAP is not proof of uniqueness',
                         'QC thresholds require acquisition-specific validation'])
        write_json(out/'qc-report.json', report)
        render_report(report, out/'report.html')
        if file_fingerprints(sources) != fingerprints:
            raise ValueError('input files changed during workflow')
        write_json(out/'run-manifest.json', dict(schema_version=1,config=config,inputs=fingerprints,
            code=code_provenance(),outputs={p.name:sha256_file(p) for p in sorted(out.iterdir()) if p.is_file()}))
    return report


def render_report(report, path):
    rows = []
    for item in report['samples']:
        link = (f'<a href="{escape(item["result_file"])}">Full result / posterior</a>'
                if item['result_file'] else '')
        maps = ', '.join(f'{p["dose_uM"]:g} µM / {p["elapsed_time_hours"]:g} h' for p in item['map_grid_pairs'])
        rows.append('<tr>'+''.join(f'<td>{escape(str(item[k]))}</td>' for k in
            ('session_id','sample_id','timepoint','status','reason','credible_set_size'))+
            f'<td>{escape(maps)}</td><td>{link}</td></tr>')
    html = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>ICHNOS run report</title>
<style>body{font:16px system-ui;max-width:1400px;margin:40px auto;padding:20px;color:#193142}table{border-collapse:collapse;width:100%}th,td{padding:10px;text-align:left;border-bottom:1px solid #ccd}pre{white-space:pre-wrap}</style>
<h1>ICHNOS run report</h1><p>Conditional model inference. Experimental accuracy has not been established.</p>'''
    html += f'<p>{report["n_cells_qc_pass"]} / {report["n_cells"]} cells passed QC.</p>'
    html += '<pre>'+escape(json.dumps(report['qc_rejection_counts'],indent=2))+'</pre>'
    html += '<table><tr>'+''.join(f'<th>{s}</th>' for s in
        ('Session','Sample','Point','Status','Reason','Credible nodes','MAP grid pairs','Artifact'))+'</tr>'+''.join(rows)+'</table>'
    html += '<p>All tied MAP points and disconnected credible sets are retained in the JSON results. A poor-fit posterior is diagnostic, not an accepted estimate.</p></html>'
    Path(path).write_text(html,encoding='utf-8')
