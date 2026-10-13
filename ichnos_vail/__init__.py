"""ichnos_vail — the ICHNOS application: a local Streamlit front end over the command-line scripts.

ichnos_vail computes nothing itself. Every action builds the same argument list a
user would type (``commands``), runs it as a background job that survives a
closed browser tab (``jobs``), and reads the files the script wrote
(``results``). The scripts and their tests stay the single source of truth.
"""
