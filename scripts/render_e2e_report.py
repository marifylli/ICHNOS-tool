"""Render an existing verified synthetic run as local HTML and PNG previews."""
import argparse
import hashlib
from html import escape
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from skimage.segmentation import find_boundaries
from ichnos_image.segment import segment_cells


def plot_curves(input_dir, out_dir):
    """Plot saved model points and calibrated synthetic observations."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from ichnos.decoder import load_dose_table

    reference = pd.read_csv(input_dir / "reference/fluorescence.csv")
    cells = pd.read_csv(input_dir / "calibration-reference-cells.csv")
    coefficient = json.loads((input_dir / "calibration.json").read_text())["c_session"]
    image_ratios = cells.groupby("timepoint").ratio_red_green.median().to_numpy() / coefficient
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for axis, field, title in zip(axes, ["green", "red", "ratio_red_green"],
                                 ["Observed green", "Mature red reporter", "mCherry/GFP model ratio"]):
        axis.plot(reference.time_hours, reference[field], "o-", label="Stored model points")
        axis.set(xlabel="Hours after exposure onset", ylabel=title, title=title)
        axis.grid(alpha=.2)
    for positions, label, marker in [(slice(None,None,2), "Calibrated image ratios: fit", "s"),
                                      (slice(1,None,2), "Calibrated image ratios: holdout", "x")]:
        axes[2].scatter(reference.time_hours.to_numpy()[positions],image_ratios[positions],
                        marker=marker,s=60,label=label,zorder=3)
    axes[2].legend(fontsize=8)
    fig.suptitle("Synthetic reference: 25 μM — model time course")
    fig.savefig(out_dir / "time-course.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
    for axis, filename, sample, true_dose, title in [
        (axes[0], "known-table.json", "dose75", 75., "Discrete dose grid"),
        (axes[1], "continuous-table.json", "dose37.5", 37.5, "Checked interpolation grid")]:
        table = load_dose_table(input_dir / filename)
        summary = pd.read_csv(input_dir / f"{sample}-summary/samples.csv")
        for column, time in enumerate(table.times_hours):
            line, = axis.plot(table.doses_uM, table.ratios[:,column], "o-", label=f"Model grid: {time:g} h")
            observed = summary.loc[np.isclose(summary.sampling_time_hours,time),"calibrated_ratio_red_green_median"].iloc[0]
            axis.scatter([true_dose],[observed],marker="x",s=90,color=line.get_color(),zorder=3,
                         label=f"Synthetic observation: {time:g} h")
        axis.set(xlabel="Dose (μM)",ylabel="Calibrated mCherry/GFP ratio",title=title)
        axis.grid(alpha=.2)
        axis.legend(fontsize=8)
    fig.suptitle("Synthetic dose–response: markers placed at known generation doses")
    fig.savefig(out_dir / "dose-response.png", dpi=160)
    plt.close(fig)


def render(input_dir, out_dir):
    input_dir, out_dir = Path(input_dir), Path(out_dir)
    if out_dir.exists():
        raise FileExistsError(f"output directory already exists: {out_dir}")
    report = json.loads((input_dir/'report.json').read_text())
    if report.get('data_kind') != 'synthetic':
        raise ValueError('this renderer requires a synthetic verification report')
    hashes = report.get('files_sha256', {})
    if not hashes:
        raise ValueError('report lacks artifact fingerprints')
    for name, digest in hashes.items():
        path = (input_dir/name).resolve()
        if not path.is_relative_to(input_dir.resolve()):
            raise ValueError('artifact path escapes input directory')
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'artifact differs from verification report: {name}')
    paths = sorted(input_dir.glob('*.npz'))
    if not paths:
        raise ValueError('no synthetic images found')
    arrays = []
    maxima = dict(green=0., red=0., bright_field=0.)
    for path in paths:
        if path.name not in hashes:
            raise ValueError('image is not fingerprinted in the verification report')
        with np.load(path, allow_pickle=False) as data:
            channels = {name:data[name].copy() for name in maxima}
        shape = channels['green'].shape
        if len(shape) != 2 or any(a.shape != shape or not np.isfinite(a).all() or (a<0).any() for a in channels.values()):
            raise ValueError('image channels must be matching finite nonnegative 2D arrays')
        for name, array in channels.items():
            maxima[name] = max(maxima[name],float(array.max()))
        arrays.append((path,channels))
    out_dir.mkdir(parents=True)
    parts = ['<!doctype html><html lang="en"><meta charset="utf-8">',
             '<title>ICHNOS synthetic verification</title>',
             '<style>body{font:16px system-ui;max-width:1100px;margin:35px auto;padding:0 20px;color:#183043}table{border-collapse:collapse;width:100%;margin:20px 0}td,th{border:1px solid #ccd6df;padding:9px;text-align:left}figure{display:inline-block;margin:8px}img{width:220px;image-rendering:pixelated}section{border-top:1px solid #ccd6df;margin-top:30px}small{color:#43556a}</style>',
             '<h1>ICHNOS — synthetic verification</h1>',
             '<p><strong>Synthetic images • software verification • no experimental validation</strong></p>',
             f'<p>Overall checks passed: {escape(str(report["passed"]))}</p>',
             '<p>Green: GFP. Red: mCherry. Bright field: input used to identify cells. Yellow outlines: segmentation boundaries. Colors show intensity, not dose.</p>',
             '<p>Each channel uses one fixed scale across every image in this run; no per-image contrast stretching. Units are arbitrary synthetic image units, not concentrations. Each PNG is a display preview; NPZ retains the numerical data.</p>',
             '<table><tr><th>Check</th><th>Passed</th><th>Expected</th><th>Actual</th></tr>']
    labels = dict(reference_cell_extraction='Reference cells detected', session_scale='Session calibration scale', calibration_holdout='Calibration at held-out times', known_time_dose='Dose at known times', joint_dose_time='Dose and elapsed time', loose_tolerance_ambiguity='Ambiguous response', continuous_dose='Continuous dose estimate', low_green_rejection='Insufficient green signal')
    for row in report['checks']:
        row = {**row, 'name': labels.get(row['name'], row['name'])}
        parts.append('<tr>'+''.join('<td>'+escape(json.dumps(row[key]))+'</td>' for key in ('name','passed','expected','actual'))+'</tr>')
    parts.append('</table>')
    plot_curves(input_dir, out_dir)
    parts.extend(['<h2>Model curves and synthetic observations</h2>',
                  '<p>Lines connect stored simulation points. They do not add new simulated times or doses. Observation markers are calibrated image ratios; dose coordinates are known synthetic generation doses, not decoder estimates. Absolute model green/red outputs are not camera intensity units.</p>',
                  '<img style="width:100%;image-rendering:auto" src="time-course.png" alt="Reference model time courses with calibrated image ratio markers">',
                  '<img style="width:100%;image-rendering:auto" src="dose-response.png" alt="Model dose response grids and calibrated synthetic observations">',
                  '<h2>Images and extracted measurements</h2>'])
    manifest = []
    for path, channels in arrays:
        sample,index = path.stem.rsplit('-',1)
        index = int(index)
        cells = pd.read_csv(input_dir/f'{sample}-cells.csv')
        selected = cells[cells.timepoint==index]
        times = selected.sampling_time_hours.unique()
        if len(times) != 1:
            raise ValueError('image time cannot be resolved from cell export')
        time = float(times[0])
        labels = segment_cells(channels['bright_field'],method='otsu')
        label_count = int(np.unique(labels[labels>0]).size)
        if label_count != len(selected):
            raise ValueError('recomputed segmentation count differs from exported cells')
        parts.append(f'<section><h3>{escape(sample)} — {time:g} h after onset</h3>')
        parts.append(f'<p>Extracted cells: {len(selected)}. QC-passing cells: {int(selected.qc_pass.sum())}. Original median mCherry/GFP: {selected.ratio_red_green.median():.6g}.</p>')
        summary_path = input_dir/f'{sample}-summary/samples.csv'
        if summary_path.exists():
            summary = pd.read_csv(summary_path)
            row = summary[summary.timepoint==index].iloc[0]
            calibrated = row.calibrated_ratio_red_green_median
            value = 'unavailable: insufficient usable cells' if pd.isna(calibrated) else f'{calibrated:.6g}'
            parts.append(f'<p>Cells used after green-floor filtering: {int(row.n_cells_used)}. Calibrated median mCherry/GFP: {value}.</p>')
        else:
            parts.append(f'<p>Calibration reference: {"fit" if index%2==0 else "holdout"} time. Known synthetic reference dose: 25 μM.</p>')
        for name,title in [('green','GFP'),('red','mCherry'),('bright_field','Bright field'),('segmentation','Segmentation')]:
            source = 'bright_field' if name=='segmentation' else name
            upper = maxima[source]
            normalized = np.zeros(channels[source].shape) if upper==0 else channels[source]/upper
            gray = np.rint(normalized*255).astype(np.uint8)
            rgb = np.zeros((*gray.shape,3),dtype=np.uint8)
            if name=='green': rgb[:,:,1]=gray
            elif name=='red': rgb[:,:,0]=gray
            else: rgb[:]=gray[:,:,None]
            if name=='segmentation': rgb[find_boundaries(labels,mode='inner')]=[255,220,0]
            filename = f'{path.stem}-{name}.png'
            Image.fromarray(rgb).save(out_dir/filename)
            parts.append(f'<figure><img src="{escape(filename,quote=True)}" alt="{escape(title)} for {escape(sample)}"><figcaption>{title}<br><small>Display range: 0–{upper:.6g}</small></figcaption></figure>')
        parts.append('</section>')
        manifest.append(dict(image=path.name,sample=sample,sampling_time_hours=time,extracted_cells=len(selected)))
    parts.append('<p>Segmentation is recomputed using the same deterministic default Otsu function as the verification pipeline. The cell count is checked against the exported cell CSV. This viewer supports these idealized synthetic runs; it is not a general microscope-image viewer.</p></html>')
    (out_dir/'index.html').write_text('\n'.join(parts),encoding='utf-8')
    (out_dir/'preview_metadata.json').write_text(json.dumps(dict(
        synthetic=True,display_maxima=maxima,images=manifest,
        source_report_sha256=hashlib.sha256((input_dir/'report.json').read_bytes()).hexdigest(),
        segmentation='recomputed default Otsu; count checked against exported cells'),indent=2)+'\n')
    return out_dir/'index.html'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--out-dir',type=Path,required=True)
    args=parser.parse_args()
    try:
        print(render(args.input_dir,args.out_dir))
    except (ValueError,KeyError,OSError) as exc:
        parser.exit(1,f'Error: {exc}\n')


if __name__=='__main__':
    main()
