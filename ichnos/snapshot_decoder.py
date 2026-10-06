"""Two-observable discrete snapshot compatibility with traceable green anchors.

Explicit tolerances define compatibility, not a likelihood or posterior.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from .calibration import _nonempty, model_reference_id, session_calibration_from_dict
from .decoder import DoseTable, _build_response_table, _digest, _vector
from .sample_decoder import check_reference_compatibility


@dataclass
class ObservableTable:
    ratio_table: DoseTable
    green: np.ndarray

    def __post_init__(self):
        values = np.asarray(self.green, dtype=object)
        if values.shape != self.ratio_table.ratios.shape:
            raise ValueError('green must have the same dose/time shape as ratios')
        self.green = _vector(values.ravel(), 'Observed_Green').reshape(values.shape)

    def to_dict(self):
        data = dict(schema_version=1, kind='two_observable_table',
                    ratio_table=self.ratio_table.to_dict(), green=self.green.tolist(),
                    green_observable='Observed_Green')
        return {**data, 'artifact_id': _digest(data)}


def build_observable_table(**kwargs):
    ratio_table, green = _build_response_table(**kwargs)
    return ObservableTable(ratio_table, green)


def save_artifact(data, path):
    encoded = json.dumps(data, indent=2, allow_nan=False) + '\n'
    with Path(path).open('x', encoding='utf-8') as handle:
        handle.write(encoded)


def _check_digest(data, kind):
    if data.get('schema_version') != 1 or data.get('kind') != kind:
        raise ValueError(f'expected {kind} schema 1')
    unsigned = {key:value for key,value in data.items() if key != 'artifact_id'}
    if data.get('artifact_id') != _digest(unsigned):
        raise ValueError(f'modified {kind} artifact')


def load_observable_table(path):
    data = json.loads(Path(path).read_text())
    _check_digest(data, 'two_observable_table')
    ratio = data['ratio_table']
    if ratio.get('schema_version') != 1 or ratio.get('artifact_id') != _digest({k:v for k,v in ratio.items() if k != 'artifact_id'}):
        raise ValueError('modified ratio table')
    if data['green_observable'] != 'Observed_Green':
        raise ValueError('unsupported green observable')
    return ObservableTable(DoseTable(ratio['doses_uM'], ratio['times_hours'], ratio['ratios'], ratio['provenance']), data['green'])


def _number(value, name, *, positive=False):
    return float(_vector([value], name, positive=positive)[0])


def _ids(values, name):
    if not isinstance(values, list) or not values:
        raise ValueError(f'{name} requires distinct nonempty identifiers')
    checked = [_nonempty(value, name) for value in values]
    if len(set(checked)) != len(checked):
        raise ValueError(f'{name} requires distinct nonempty identifiers')
    return checked


def _acquisition(value):
    if not isinstance(value, dict):
        raise ValueError('acquisition must be a settings dictionary')
    for name in ('objective', 'green_units', 'extraction_green', 'extraction_red', 'gain_setting'):
        _nonempty(value.get(name), name)
    for name in ('exposure_ms_green', 'exposure_ms_red'):
        _number(value.get(name), name, positive=True)
    json.dumps(value, allow_nan=False)
    return value


def _anchor(table, reference):
    for name in ('specimen_id', 'biological_replicate_id', 'source'):
        _nonempty(reference.get(name), name)
    dose = _number(reference['dose_uM'], 'reference dose')
    time = _number(reference['time_hours'], 'reference time')
    image_green = _number(reference['corrected_green'], 'reference corrected green')
    rows = np.flatnonzero(np.isclose(table.ratio_table.doses_uM, dose, rtol=0, atol=1e-12))
    columns = np.flatnonzero(np.isclose(table.ratio_table.times_hours, time, rtol=0, atol=1e-12))
    if len(rows) != 1 or len(columns) != 1:
        raise ValueError('reference dose/time must identify exactly one stored model point')
    return image_green, float(table.green[rows[0], columns[0]])


def build_snapshot_calibration(
    table, *, session_id, source_kind, acquisition, low_reference, high_reference,
    min_image_span, min_model_span, ratio_calibration, ratio_reference_metadata,
    ratio_reference_specimen_ids, ratio_reference_replicate_ids,
):
    """Bind two fixed-time green anchors and an existing ratio calibration.

    IDs declare independent reference specimens/cultures; software cannot
    verify their experimental provenance. Anchors are not uncertainty estimates.
    """
    _nonempty(session_id, 'session_id')
    if source_kind not in {'synthetic', 'experimental_reference'}:
        raise ValueError('source_kind must be synthetic or experimental_reference')
    _acquisition(acquisition)
    ratio = session_calibration_from_dict(ratio_calibration)
    if ratio.session_id != session_id or ratio.source_kind != source_kind:
        raise ValueError('ratio calibration session/source differs')
    if model_reference_id(ratio_reference_metadata) != ratio.reference_id:
        raise ValueError('ratio reference identity differs')
    check_reference_compatibility(table.ratio_table, ratio_reference_metadata)
    _ids(ratio_reference_specimen_ids, 'ratio_reference_specimen_ids')
    _ids(ratio_reference_replicate_ids, 'ratio_reference_replicate_ids')
    image_low, model_low = _anchor(table, low_reference)
    image_high, model_high = _anchor(table, high_reference)
    if low_reference['dose_uM'] != 0:
        raise ValueError('low reference must be unstressed (dose zero)')
    if high_reference['dose_uM'] <= 0:
        raise ValueError('high reference must have a positive declared dose')
    if low_reference['specimen_id'] == high_reference['specimen_id']:
        raise ValueError('green references must be distinct specimens')
    min_image_span = _number(min_image_span, 'min_image_span', positive=True)
    min_model_span = _number(min_model_span, 'min_model_span', positive=True)
    if image_high-image_low < min_image_span or model_high-model_low < min_model_span:
        raise ValueError('green reference span is too small or reversed')
    data = dict(
        schema_version=1, kind='snapshot_calibration', table_artifact_id=table.to_dict()['artifact_id'],
        session_id=session_id, source_kind=source_kind, acquisition=acquisition,
        low_reference=low_reference, high_reference=high_reference,
        min_image_span=min_image_span, min_model_span=min_model_span,
        ratio_calibration=ratio_calibration, ratio_reference_metadata=ratio_reference_metadata,
        ratio_reference_specimen_ids=ratio_reference_specimen_ids,
        ratio_reference_replicate_ids=ratio_reference_replicate_ids,
        image_green_anchors=[image_low,image_high], model_green_anchors=[model_low,model_high],
        green_equation='(green - low_reference) / (high_reference - low_reference)',
        high_reference_saturation_validated=False, reference_independence='declared IDs only',
        experimentally_validated=False, uncertainty_estimated=False,
    )
    # Detach mutable caller dictionaries before signing.
    data = json.loads(json.dumps(data, allow_nan=False))
    return {**data, 'artifact_id':_digest(data)}


def validate_snapshot_calibration(table, calibration):
    _check_digest(calibration, 'snapshot_calibration')
    keys = ('session_id','source_kind','acquisition','low_reference','high_reference',
            'min_image_span','min_model_span','ratio_calibration','ratio_reference_metadata',
            'ratio_reference_specimen_ids','ratio_reference_replicate_ids')
    expected = build_snapshot_calibration(table, **{key:calibration[key] for key in keys})
    if expected != calibration:
        raise ValueError('calibration differs from table or recorded reference mapping')


def prepare_snapshot(table, calibration, observation, *, green_floor):
    """Validate identity/units and normalize a single snapshot for inference."""
    validate_snapshot_calibration(table, calibration)
    green_floor = _number(green_floor, 'green_floor')
    for name in ('session_id','specimen_id','biological_replicate_id'):
        _nonempty(observation.get(name), name)
    if observation['session_id'] != calibration['session_id']:
        raise ValueError('observation session differs')
    expected_kind = 'synthetic' if calibration['source_kind']=='synthetic' else 'experimental'
    if observation.get('data_kind') != expected_kind:
        raise ValueError('observation data kind differs from calibration source')
    if observation.get('ratio_scale') != 'image':
        raise ValueError('snapshot requires original image ratio, ratio_scale=image')
    if _acquisition(observation.get('acquisition')) != calibration['acquisition']:
        raise ValueError('observation acquisition settings differ from references')
    references = [calibration['low_reference'], calibration['high_reference']]
    specimen_ids = calibration['ratio_reference_specimen_ids'] + [r['specimen_id'] for r in references]
    culture_ids = calibration['ratio_reference_replicate_ids'] + [r['biological_replicate_id'] for r in references]
    if observation['specimen_id'] in specimen_ids or observation['biological_replicate_id'] in culture_ids:
        raise ValueError('target overlaps reference specimen or biological replicate IDs')
    image_ratio = _number(observation['image_ratio_red_green'], 'image_ratio_red_green')
    image_green = _number(observation['corrected_green'], 'corrected_green')
    ratio = image_ratio / calibration['ratio_calibration']['c_session']
    ilow, ihigh = calibration['image_green_anchors']
    mlow, mhigh = calibration['model_green_anchors']
    green = (image_green-ilow)/(ihigh-ilow)
    predicted_green = (table.green-mlow)/(mhigh-mlow)
    if not np.isfinite([ratio,green]).all() or not np.isfinite(predicted_green).all():
        raise ValueError('normalized observations/predictions must be finite')
    result = dict(
        schema_version=1, table_artifact_id=table.to_dict()['artifact_id'],
        calibration_artifact_id=calibration['artifact_id'],
        observation=json.loads(json.dumps(observation,allow_nan=False)),
        ratio_red_green=ratio, green_norm=green, green_norm_outside_reference_interval=bool(green<0 or green>1),
        green_was_clipped=False,
        green_floor_image_units=green_floor, candidate_pairs=[],
        dose_estimate_uM=None, elapsed_time_estimate_hours=None,
        inference='two-observable rectangular compatibility on a discrete grid',
        uncertainty_estimated=False, posterior_computed=False, experimentally_validated=False,
        reference_independence='declared IDs checked; experimental provenance not independently verified',
        decoder_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    result['status'] = 'below_detection_floor' if image_green <= green_floor else 'ready'
    return result, predicted_green


def decode_snapshot(table, calibration, observation, *, ratio_tolerance, green_tolerance, green_floor):
    """Decode one snapshot using rectangular two-observable compatibility."""
    ratio_tolerance = _number(ratio_tolerance, 'ratio_tolerance', positive=True)
    green_tolerance = _number(green_tolerance, 'green_tolerance', positive=True)
    result, predicted_green = prepare_snapshot(table, calibration, observation, green_floor=green_floor)
    result.update(ratio_tolerance=ratio_tolerance, green_tolerance=green_tolerance)
    if result['status'] == 'below_detection_floor':
        return result
    ratio, green = result['ratio_red_green'], result['green_norm']
    ratio_ok = np.abs(table.ratio_table.ratios-ratio) <= ratio_tolerance
    green_ok = np.abs(predicted_green-green) <= green_tolerance
    result['ratio_only_candidate_count'] = int(ratio_ok.sum())
    result['green_only_candidate_count'] = int(green_ok.sum())
    matches = np.argwhere(ratio_ok & green_ok)
    result['candidate_pairs'] = [dict(dose_uM=float(table.ratio_table.doses_uM[i]),
        elapsed_time_hours=float(table.ratio_table.times_hours[j]),
        predicted_ratio=float(table.ratio_table.ratios[i,j]), predicted_green_norm=float(predicted_green[i,j]))
        for i,j in matches]
    if len(matches)==1:
        match=result['candidate_pairs'][0]
        return {**result,'status':'unique_grid_pair','dose_estimate_uM':match['dose_uM'],
                'elapsed_time_estimate_hours':match['elapsed_time_hours']}
    if len(matches)>1:
        return {**result,'status':'ambiguous'}
    outside = (ratio < table.ratio_table.ratios.min()-ratio_tolerance or
               ratio > table.ratio_table.ratios.max()+ratio_tolerance or
               green < predicted_green.min()-green_tolerance or green > predicted_green.max()+green_tolerance)
    return {**result,'status':'out_of_response_domain' if outside else 'no_grid_pair_match'}
