import csv

import pytest

from ichnos_image.acquisition_log import COLUMNS, parse_log_text, write_rows

# Excerpt in the shape pdftotext -layout produces, including a form feed page break.
LOG = """Control σε θρεπτικό
Μέρα    Brightfield    Exp    Red    Green
5/10    17:30          375    17:32  17:32

Οξειδωτικό στρες (H2O2) — χρονοσειρά 5–6/10
Για την μετρηση στις 05/10 το H2O2 μπηκε 13:41.
Για τις μετρησεις στις 06/10 το Η2Ο2 μπηκε 11:55.
0h     5/10   0     14:28   400   14:32   14:33
3h     5/10   300   17:15   400   17:13   17:14
              (cluster)
\f6h     6/10   600   18:09   300   18:07   18:07

ER στρες (DTT) — χρονοσειρά 6/10
To DTT μπηκε μια φορά τις 12:49
1h   13:49*   100 (cluster) 15:14   300   15:13   15:13
1h   13:49*   100           15:15   250   15:14   15:14
              (normal)
6h   19:49*   0             19:09   400   19:08   19:08

ER στρες (DTT) — τριπλέτες 7/10
Το DTT προστέθηκε 12:50.
Στρώση 16:50
Exp = 200 σε όλες τις λήψεις.
0      R1     17:34     17:31     17:34
0      R1     17:34     17:31     17:34
"""


@pytest.fixture
def rows():
    return {(r["section"], r["tp_label"], r["dose_uM"], r["rep"]): r for r in parse_log_text(LOG).rows}


def test_series_rows_get_their_own_days_addition_time(rows):
    first = rows[("h2o2_series", "0h", "0", "")]
    late = rows[("h2o2_series", "6h", "600", "")]
    assert (first["added"], first["measurement_time_hours"]) == ("13:41", "0.850")
    assert (late["added"], late["measurement_time_hours"]) == ("11:55", "6.200")


def test_row_after_form_feed_is_still_parsed(rows):
    assert ("h2o2_series", "6h", "600", "") in rows


def test_cluster_and_normal_tags_on_the_same_or_next_line(rows):
    assert rows[("h2o2_series", "3h", "300", "cluster")]["red"] == "17:13"
    assert rows[("dtt_0610", "1h", "100", "cluster")]["green"] == "15:13"
    assert rows[("dtt_0610", "1h", "100", "normal")]["green"] == "15:14"


def test_listed_exposure_is_used_and_unlisted_is_derived_and_flagged(rows):
    listed = rows[("h2o2_series", "0h", "0", "")]
    derived = rows[("dtt_0610", "1h", "100", "normal")]
    assert listed["exposure_acq_ms"] == "1600" and "acq_derived_4x" not in listed["flags"]
    assert derived["exposure_acq_ms"] == "1000" and "acq_derived_4x" in derived["flags"]


def test_exposure_from_section_note_is_flagged(rows):
    row = rows[("dtt_0710_s1650", "4h", "0", "R1")]
    assert row["exp"] == "200" and "exp_from_section_note" in row["flags"]


def test_impossible_plating_time_is_corrected_and_flagged(rows):
    row = rows[("dtt_0610", "6h", "0", "")]
    assert row["plating"] == "18:49" and row["sampling_time_hours"] == "6.000"
    assert "plating_19:49_read_as_18:49" in row["flags"]


def test_duplicate_log_row_is_flagged_not_dropped():
    result = parse_log_text(LOG)
    dup = [r for r in result.rows if r["section"] == "dtt_0710_s1650"]
    assert len(dup) == 2 and "duplicate_row_in_log" in dup[1]["flags"]


def test_misfiled_2h_folders_follow_the_log():
    text = """.η μετασχηματισμένα (DTT) — τριπλέτες 9/10
Το DTT μπήκε 12:20
2 ΩΡΕΣ
100   R1   15:34   15:33   15:33
50    R1   15:40   15:39   15:39
200   R1   15:46   15:45   15:45
"""
    rows = {r["dose_uM"]: r for r in parse_log_text(text).rows}
    assert rows["100"]["folder"] == "D:igem DTT non transfected/200/2h"
    assert rows["200"]["folder"] == "D:igem DTT non transfected/100/2h"
    assert rows["50"]["folder"] == "D:igem DTT non transfected/50/2h"
    assert "log_authoritative" in rows["100"]["flags"]


def test_doses_are_micromolar(rows):
    assert rows[("h2o2_series", "0h", "0", "")]["dose_units"] == "uM"


def test_missing_headings_are_reported_not_ignored():
    result = parse_log_text(LOG)
    assert any("nt_h2o2_0810" in w for w in result.warnings)


def test_write_refuses_to_overwrite(tmp_path):
    out = tmp_path / "log.csv"
    write_rows(parse_log_text(LOG).rows, out)
    with out.open(encoding="utf-8") as handle:
        assert next(csv.reader(handle)) == COLUMNS
    with pytest.raises(FileExistsError):
        write_rows([], out)
