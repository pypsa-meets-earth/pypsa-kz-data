# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Description
"""
import logging
import os
import pypsa
import pandas as pd
import sys
sys.path.append("scripts/")
from _helpers import configure_logging

logger = logging.getLogger(__name__)


def rescale_load(n, gadm_demand):
    """
    Rescales demand profiles of each GADM region based on the official 
    national report.

    Parameters
    ----------
    n: pypsa network

    gadm_demand: str
        Path to csv file containing officcial consumption by GADM regions, 
        e.g. "data/custom_data/official_demand.csv"
        Regional demands is reported in GWh.

    Returns
    -------
    n : pypsa network
        Now updated with rescaled demand profiles for each GADM zone
    """
    logger.info(f"Rescaling GADM demand profile based on national report")

    official_demand = pd.read_csv(gadm_demand, index_col=0)
    
    pypsa_demand = n.loads_t.p_set.sum(axis=0)
    
    scale_factor = (official_demand["Consumption (GWh)"] * 1e3) / pypsa_demand
    
    n.loads_t.p_set.loc[:,official_demand.index] = (
        n.loads_t.p_set.loc[:,official_demand.index] * scale_factor.loc[official_demand.index]
    )

    if True in n.generators.query("carrier=='coal'").p_nom_extendable.unique():
        n.generators.loc[n.generators.carrier=="coal","p_nom_max"] = n.generators.loc[n.generators.carrier=="coal","p_nom"]
        n.generators.loc[n.generators.carrier=="coal", "p_nom"] = 0
        n.generators.loc[n.generators.carrier=="coal", "p_nom_min"] = 0


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("modify_demand")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
    
    # Snakemake imports:
    gadm_demand = snakemake.input["gadm_demand_data"]
    
    rescale_load(n, gadm_demand)

    # Snakemake output
    n.export_to_netcdf(snakemake.output.network)
