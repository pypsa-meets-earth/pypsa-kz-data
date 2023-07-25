# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Adds external nodes and assigns the demand profiles to external buses.
"""
import logging
import os
import pypsa
import pandas as pd
from calendar import monthrange
import sys
sys.path.append("scripts/")
from _helpers import configure_logging


logger = logging.getLogger(__name__)


def add_external_nodes(n, data_imp_exp):
    """
    Add external nodes with demand profile build from data_imp_exp.

    Parameters
    ----------
    n: pypsa network
    
    data_imp_exp: str
        Path to csv file containing import and exports data from neighboring
        countries, e.g. "data/custom_data/electricity_exp-imp.csv"
        Demands are reported in GWh.

    Returns
    -------
    n : pypsa network
        Network with added external buses
    """
    logger.info(f"Adding external buses")

    imp_exp_elec = pd.read_csv(data_imp_exp, sep=",", decimal=",")

    # Trimming data to a single year exp-exports, imp-imports
    exp_elec = imp_exp_elec.iloc[:12,1:].iloc[:,1::2]
    imp_elec = imp_exp_elec.iloc[:12,1:].iloc[:,0::2]

    # Renaming the columns
    exp_elec.columns = ["Russia", "Kyrgyzstan", "Uzbekistan"]
    imp_elec.columns = ["Russia", "Kyrgyzstan", "Uzbekistan"]

    # Export import difference 
    external_loads = (exp_elec - imp_elec) * 1e3 # convert to MWh

    # initialize external load profile
    ex_load_profile = pd.DataFrame()
    ex_load_profile.index = n.snapshots

    for i in external_loads.index:
        year = n.snapshots.year[0]
        _, days = monthrange(int(year), int(i)+1)
        hourly_load = external_loads.loc[i]/(days*24)
        decimal = [0 if (int(i) + 1) < 10 else ""][0]
        ex_load_profile.loc[f"{year}-{decimal}{int(i)+1}",["RU", "KG", "UZ"]] = hourly_load.values

    # creating three buses
    y_RU, x_RU = 55.064853, 61.509004
    y_KG, x_KG = 42.721375, 74.544894
    y_UZ, x_UZ = 40.950348, 69.377086
    external_buses = n.buses.loc[["KZ.12_1_AC", "KZ.12_1_AC", "KZ.12_1_AC"]].copy()
    external_buses.index = ["RU", "KG", "UZ"]
    external_buses.x = [x_RU, x_KG, x_UZ]
    external_buses.y = [y_RU, y_KG, y_UZ]

    n.import_components_from_dataframe(external_buses, "Bus")

    # connecting lines
    external_lines = n.lines.loc[["5","5","5"]].copy() # external lines is set to have same property as line 5
    external_lines.bus0 = ["KZ.10_1_AC", "KZ.1_1_AC", "KZ.12_1_AC"]
    external_lines.bus1 = ["RU", "KG", "UZ"]
    external_lines.index = ["18", "19","20"]

    n.import_components_from_dataframe(external_lines, "Line")

    # attaching loads
    external_loads = n.loads.loc[["KZ.12_1_AC","KZ.12_1_AC","KZ.12_1_AC"],:].copy()
    external_loads.index = ["RU", "KG", "UZ"]
    external_loads.bus = ["RU", "KG", "UZ"]

    n.import_components_from_dataframe(external_loads, "Load")

    n.import_series_from_dataframe(ex_load_profile, "Load", "p_set")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("add_exports")
    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)
    
    # Snakemake input
    data_imp_exp = snakemake.input.imp_exp

    add_external_nodes(n, data_imp_exp)

    # Snakemake output
    n.export_to_netcdf(snakemake.output.network)
