# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  pypsa-kz-data authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Modifies the line capacities of buses of interest (boi)
"""
import logging
import os
import pypsa
import pandas as pd
import sys
sys.path.append("scripts/")
from _helpers import configure_logging


logger = logging.getLogger(__name__)


def limit_line_capacities(n, boi, max_limit):
    """
    Limits s_nom of lines attached to buses of interest (boi) by max_limit.

    Parameters
    ----------
    n: pypsa network
    
    boi: list of str
        Buses of interest to which lines being limited are attached, 
        e.g. ["KZ.14.1_AC", "KZ.12.1_AC"]

    max_limit: float
        Maximum line capacity s_nom to which lines are limited. The unit of 
        measure is MVA, e.g. 700.0

    Returns
    -------
    n : pypsa network
        Network with updated line capacities (s_nom)
    """
    logger.info(f"Limiting s_nom by {max_limit} MVA for lines connected to \
                {boi} buses")

    for line in n.lines.index:
        if n.lines.loc[line].bus0 in boi or n.lines.loc[line].bus1 in boi:
            if n.lines.loc[line, 's_nom'] >= max_limit:
                n.lines.loc[line, 's_nom'] = max_limit

def add_line(n):
    extra_line = n.lines.loc[["2"]].copy()
    extra_line.bus0 = "KZ.3_1_AC"
    extra_line.bus1 = "KZ.4_1_AC"
    extra_line.index = ["21"]
    n.import_components_from_dataframe(extra_line, "Line")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("modify_lines")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
        
    boi = snakemake.config["lines"]["modify_lines"]["bus_of_interest"]
    max_limit = snakemake.config["lines"]["modify_lines"]["max_limit"]

    limit_line_capacities(n, boi, max_limit)

    if snakemake.config["lines"]["modify_lines"].get('new_line_kz', False):
        add_line(n)

    # Snakemake output
    n.export_to_netcdf(snakemake.output.network)
