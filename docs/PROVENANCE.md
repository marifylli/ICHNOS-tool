# Provenance

Every file in this repository that was taken from another repository is listed
here with the exact source commit and the SHA256 of the source file as it was
at that commit. Nothing is copied without an entry here.

## Source repositories (frozen commits)

| Repository | Commit | Role |
| --- | --- | --- |
| [DryLabTool](https://github.com/YannisBalasis/DryLabTool/tree/2692ea89c5cab24eb35ec628d6ccde2bac68f2ca) | `2692ea89` | image pipeline and shared schema |
| [Ichnos_PULSE](https://github.com/ioanna888/Ichnos_PULSE/tree/e66de65cb7ca7c170bfff96c181577f5e9df4db3) | `e66de65c` | SBML builder and the four source models (Step 3) |
| [ICHNOS-ablation](https://github.com/marifylli/ICHNOS-ablation/tree/186a50a3bc7b5855667e0f6b46835c5ea95d18e0) | `186a50a3` | reusable Jacobian / population-noise functions (Step 7) |
| [ichnos-fisher](https://github.com/vicdo18/ichnos-fisher/tree/caa971d196ce502938a0f285a06b9125bf6f338e) | `caa971d1` | noise methods and documented limits (Step 7) |

## Step 1 — migrated from DryLabTool `2692ea89`

| Source path | Path here | Change |
| --- | --- | --- |
| `ichnos/__init__.py` | `ichnos/__init__.py` | verbatim |
| `ichnos/schema.py` | `ichnos/schema.py` | verbatim |
| `ichnos/config.py` | `ichnos/config.py` | verbatim |
| `ichnos_image/__init__.py` | `ichnos_image/__init__.py` | verbatim |
| `ichnos_image/schema.py` | `ichnos_image/schema.py` | verbatim (re-export shim) |
| `ichnos_image/segment.py` | `ichnos_image/segment.py` | verbatim |
| `ichnos_image/correct.py` | `ichnos_image/correct.py` | verbatim |
| `ichnos_image/extract.py` | `ichnos_image/extract.py` | verbatim |
| `ichnos_image/export.py` | `ichnos_image/export.py` | verbatim |
| `ichnos_image/metadata.py` | `ichnos_image/metadata.py` | verbatim |
| `ichnos_image/pipeline.py` | `ichnos_image/pipeline.py` | verbatim |
| `ichnos_image/synthesize.py` | `ichnos_image/synthesize.py` | verbatim |
| `scripts/run_pipeline.py` | `scripts/run_pipeline.py` | verbatim |
| `scripts/calibrate_focus_threshold.py` | `tools/calibrate_focus_threshold.py` | verbatim |
| `scripts/calibrate_registration_threshold.py` | `tools/calibrate_registration_threshold.py` | verbatim |
| `scripts/calibrate_rolling_ball_radius.py` | `tools/calibrate_rolling_ball_radius.py` | verbatim |
| `scripts/validate_segmentation_quantitative.py` | `tests/helpers/segmentation_scenes.py` | extracted: _make_scene + _match_and_score only, bodies unchanged |
| `tests/test_correct.py` | `tests/image/test_correct.py` | verbatim |
| `tests/test_export.py` | `tests/image/test_export.py` | verbatim |
| `tests/test_extract.py` | `tests/image/test_extract.py` | verbatim |
| `tests/test_focus_score.py` | `tests/image/test_focus_score.py` | verbatim |
| `tests/test_ichnos_core.py` | `tests/image/test_ichnos_core.py` | verbatim |
| `tests/test_metadata.py` | `tests/image/test_metadata.py` | verbatim |
| `tests/test_pipeline.py` | `tests/image/test_pipeline.py` | verbatim |
| `tests/test_pipeline_orchestration.py` | `tests/image/test_pipeline_orchestration.py` | verbatim |
| `tests/test_segment_cellpose.py` | `tests/image/test_segment_cellpose.py` | verbatim |
| `tests/test_segment_resize.py` | `tests/image/test_segment_resize.py` | verbatim |
| `tests/test_segmentation_quality.py` | `tests/image/test_segmentation_quality.py` | import rewired to tests.helpers.segmentation_scenes |
| `tests/test_synthesize.py` | `tests/image/test_synthesize.py` | verbatim |

### SHA256 of the source files (at `2692ea89`)

```
75d9f053e4be53c4feabfca1ccb88bb41cf03c0df611389beb73c0c2498ab45c  ichnos/__init__.py
5a3e75e2f27d3d50be57487cd8cd85b6ceba262cddb3a39ea3aeed535162d222  ichnos/schema.py
4909bf23342441e13b468b380ee5038cf25625e526e684d8343182ac9a89c614  ichnos/config.py
3b85c06bf889c83e5a7c8318ad6409431f0daae2ecfb69bb9ca54f16d6b23b19  ichnos_image/__init__.py
8c2ec827a1c1146c7ec959a9c738797e347fd7e1432a6b18f1803b2b79072cc7  ichnos_image/schema.py
4c7224498633b0281e53bce6d108c3b4a06a157166ca576ea1d4840d3a330d0a  ichnos_image/segment.py
b0f68973e0f0d74186979ff57d57a33c5c99a4d1835ad81d41cda942ed11d93f  ichnos_image/correct.py
9b2a1bf4d0f27103c678257a0b93e8c40b9601719f600d3e3792bd8c08fe2e0b  ichnos_image/extract.py
327d34d9df91d5cdfd00250734ac0503cd9e8c2db799afcc7db1069e3311e29b  ichnos_image/export.py
65ae32c9cf999d6340e11e7be550d965ccd56cd87f029b6e0ebfdbfa70281642  ichnos_image/metadata.py
9db98de1459c3a35a4e291ca204590209fe383564eac1a7e5a172ae9c06c84e7  ichnos_image/pipeline.py
6c4071af994b9a8c1e67a31495a0ca7d33a68dcba0730e01393964c9d50cd2cc  ichnos_image/synthesize.py
82ebd09c5860468ef75c87dc3b92856c42766647f48ebc2bdb2d861d77c91c62  scripts/run_pipeline.py
dac60dcd4aac6ee30124b25497717265d16b691829b78cded70546214ecf337d  scripts/calibrate_focus_threshold.py
033a91ff717c641bfe20f6b8e433c0a836d0398d95366acbe7e403a5c45ea4d1  scripts/calibrate_registration_threshold.py
0960baa7d87164128dd3ed2a46585ec9985cfe8cc6bd91ef1b8255130d1a6307  scripts/calibrate_rolling_ball_radius.py
6120d41a1196da41ae974dad1055dcb93b133fc3d7efa93bee919647e1aa484a  scripts/validate_segmentation_quantitative.py
df94197e1f0b79dbdcba4caba09126ee19616f5f495a65be171d5e6a381acca6  tests/test_correct.py
8b146322f34d4e6cf3b8058f8356119e43a897a5ddf12888197b4d7fb87d49c1  tests/test_export.py
1951670d309c69f1e8d47be4651792f19d86112f098189e603d9511467acd475  tests/test_extract.py
86d64a9706d3debe4cd9fc329b60352d24d262a8a28e6c6110a7ae2b24cde721  tests/test_focus_score.py
3787e6cf5c4d4d69b4cd97d55bb7dd74cd0f1f48c611efec6ee2cb9c491981b2  tests/test_ichnos_core.py
bc9e96a89fdab136efbac5a1148cf708bef67636857e63df1857e4a3da850e3f  tests/test_metadata.py
52da2e5e34cbe386254ef9b159451a90fe141c4cf827515be2efd44736a4c941  tests/test_pipeline.py
6ad14ebf1b28706f1b986fdaa3bfef1c4a782cba7acef4319fa60d517e520a28  tests/test_pipeline_orchestration.py
5a280b392ff4dcb4b0e5cb497ce0fcb1991da92961abe812079ca436c0cb7e73  tests/test_segment_cellpose.py
bb45db6074a2950d08eb598e32accc46a693bc22382aa6e3267110ec1a519556  tests/test_segment_resize.py
afb931c657dc3b9ea61198aeafb4cdd45eef841f61660dc99d651cabd4dd6fcd  tests/test_segmentation_quality.py
f17af34bd2dadc3ff8ba71ab9aac35ef90ca9b334ae80496c1e31ce2561dc9db  tests/test_synthesize.py
```

## Licensing

`ICHNOS-ablation` is MIT. `DryLabTool`, `Ichnos_PULSE` and `ichnos-fisher` carry
no licence file at the commits above. Written permission from each repository
owner, or a licence added upstream, is required before this repository is made
public or released. Tracked as an open item; see `docs/limitations.md`.

