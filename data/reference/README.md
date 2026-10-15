# Reference cell table (October 2026)

`cells_all_per_s.csv.gz` is the exposure-normalised cell table of the team's
October 2026 imaging experiment. The decoders are trained on it, and every
analysis in the README and on the wiki can be rerun from it.

- 67,580 cells, 47 columns; pandas reads the compressed file directly
  (`pd.read_csv("data/reference/cells_all_per_s.csv.gz")`).
- SHA256 of the uncompressed CSV:
  `a7586dbb84537a9495863f83a3f3c2754c91954461775c50fffc0dead591a94f`
- Column definitions: [docs/data_contract.md](../../docs/data_contract.md).
- How it was produced (bright-field alignment, Cellpose segmentation,
  exposure normalisation), with every command: [docs/run_log.md](../../docs/run_log.md).

| Session | Stressor | Cultures | Cells |
| --- | --- | --- | ---: |
| `20261005_h2o2_tr` | H2O2, 0–600 µM | transfected | 2,230 |
| `20261006_h2o2_tr` | H2O2, 0–600 µM | transfected | 1,737 |
| `20261007_h2o2_tr` | H2O2, 0–600 µM | transfected | 3,328 |
| `20261006_dtt_tr` | DTT, 0–200 µM | transfected | 11,138 |
| `20261007_dtt_tr` | DTT, 0–200 µM | transfected | 11,629 |
| `20261008_h2o2_nt` | H2O2 | non-transfected (autofluorescence) | 15,490 |
| `20261009_dtt_nt` | DTT | non-transfected (autofluorescence) | 22,028 |

The raw microscope images are not in this repository because of their size.
