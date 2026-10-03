"""Turning SBML ids into readable names, for display.

Migrated from Ichnos_PULSE python/ichnos_io.py @ e66de65c, function bodies
unchanged. Kept separate from io.py because these are presentation helpers,
not artifact I/O.

These are DISPLAY helpers only. The merge logic and any simulator call must
keep using real SBML ids. The strict name-to-id resolver that errors on a
missing or ambiguous name is a Step 4 item and is deliberately not here yet:
_find_id_by_name() below still returns None when a name is absent, which is
the migrated behaviour.
"""
import libsbml

def _build_id_to_name_map(sbml_str):
    """Maps every species/parameter SBML id to its human-readable Name (same
    GUID-vs-Name situation as find_param_value_anywhere). Used purely for
    DISPLAY (plot legends, printed values) — never for the merge logic
    itself, which must keep using real ids."""
    doc = libsbml.readSBMLFromString(sbml_str)
    m = doc.getModel()
    id_to_name = {}
    for s in m.getListOfSpecies():
        id_to_name[s.getId()] = s.getName() or s.getId()
    for p in m.getListOfParameters():
        id_to_name[p.getId()] = p.getName() or p.getId()
    return id_to_name


def _relabel_result_columns(result, id_to_name):
    """Rewrites a roadrunner simulate() result's column headers from raw
    SBML ids (e.g. '[mw7303685c_10d1_49f1_b0ab_a24defe78516]') to their
    human-readable Names (e.g. '[TetR_active]'), so plot legends and printed
    'Final values' are actually readable. roadrunner's NamedArray.colnames is
    directly settable — this doesn't touch the underlying data, only labels."""
    new_colnames = []
    for col in result.colnames:
        if col == "time":
            new_colnames.append(col)
            continue
        inner = col.strip("[]")
        new_colnames.append(f"[{id_to_name.get(inner, inner)}]")
    result.colnames = new_colnames
    return result


def _find_id_by_name(sbml_str, target_name):
    """Reverse of _build_id_to_name_map: given a human-readable Name, finds
    the actual SBML id to use with roadrunner's r[id] getter/setter (which
    operates on real ids, not Names). Returns None if not found."""
    doc = libsbml.readSBMLFromString(sbml_str)
    m = doc.getModel()
    for p in m.getListOfParameters():
        if (p.getName() or p.getId()) == target_name:
            return p.getId()
    for s in m.getListOfSpecies():
        if (s.getName() or s.getId()) == target_name:
            return s.getId()
    return None
