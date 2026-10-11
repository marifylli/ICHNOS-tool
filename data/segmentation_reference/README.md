# Segmentation reference: hand-marked cells

Seven bright-field crops (700 × 700 px, 0.1725 µm/px), one per session.
`crops.csv` says which image and where in it each crop was cut
(`y0`, `x0` = top-left corner in the full frame). The crops are
contrast-stretched for viewing only.

The marks are the reference that bright-field masks are scored against
(recall: marked cells that got a mask; precision: masks that hold a mark).
Mark what you see, not what you think the software will find.

## How to mark (Fiji / ImageJ, ~10 min per crop)

1. Open the PNG. `Analyze → Set Scale…` → *Click to Remove Scale* (so X, Y are pixels).
2. `Analyze → Set Measurements…` → tick only *Centroid*.
3. Choose the **Multi-point** tool (right-click the point tool). Click once in the
   centre of every cell.
4. `Analyze → Measure` (Ctrl/Cmd + M). Save the Results table as
   `<crop name>.points.csv` (e.g. `20261007_h2o2_tr.points.csv`) in this folder.

## What counts as a cell

- Every cell with its rim visible, whether bright or dark, round or elongated.
- A bud is a separate cell **only** if a clear neck separates it from the mother;
  otherwise mark the pair once.
- Cells cut by the crop edge: mark them if their centre is inside the crop.
- Do not mark debris, dust, the air bubble in `20261009_dtt_nt`, or blurred
  shapes with no rim.
- If unsure, mark it and add a second file `<crop name>.unsure.csv` with those
  points; they will be left out of the score.
