<!--
SPDX-FileCopyrightText:  pypsa-kz-data authors

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# About the Project

<img src="https://github.com/pypsa-meets-earth/pypsa-kz-data/assets/53824825/ca7893de-26e2-47ad-a3e4-d91cd6716652" alt="Open Energy Transition Logo" width="280" height="100" align="right">
<img src="https://github.com/pypsa-meets-earth/pypsa-kz-data/assets/53824825/63bd0250-c54a-4ce1-8df3-eb116baac01b" alt="Agora Energiewende Logo" width="240" height="100">
<br>

Agora Energiewende aims to model the Kazakh power system, incorporating a substantial increase in variable generation, such as solar and wind, surpassing the current official mid-term policy goal of 15% of all renewable energy sources (RES) in generation by 2030. This endeavor has received support from Open Energy Transition on the modeling side.

## Development status: Active and Stable

[![CI-Linux](https://github.com/pypsa-meets-earth/pypsa-kz-data/actions/workflows/ci-linux.yml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-kz-data/actions/workflows/ci-linux.yml)
![Size](https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-kz-data?label=Repo%20size)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![REUSE status](https://api.reuse.software/badge/github.com/pypsa-meets-earth/pypsa-kz-data)](https://api.reuse.software/info/github.com/pypsa-meets-earth/pypsa-kz-data)

# pypsa-kz-data
Extra data for the Kazakhstan model that will be used as input for PyPSA-Earth.
Repo design oriented on: https://github.com/pypsa-meets-earth/pypsa-zm-data

![image](https://user-images.githubusercontent.com/61968949/231397315-bc490876-abb6-45c4-bf01-e26f90c9db93.png)

# Data

Contains openly available data for Kazakhstan

## Demand data

Monthly electricity demand data with monthly aggregation provided by Kazakhstan operator of the electric energy and power market

`kz_demand_validation.csv`
(`demand_valid_clean.R` contains for details on how data were extracted and a simple vizualization; needs to work any R terminal + installation of two libraries)

# Model Execution

The provided workflow builds on [PyPSA-Earth](https://github.com/pypsa-meets-earth/pypsa-earth). Therefore, first, the PyPSA-Earth repository must be forked and the fork should then be cloned. A fork can be created by navigating to the [PyPSA-Earth](https://github.com/pypsa-meets-earth/pypsa-earth) website. By clicking on the fork-symbol in the upper right corner, a fork is created and linked to the specific user.

Next, we also need to fork the [pypsa-kz-data](https://github.com/pypsa-meets-earth/pypsa-kz-data) repository. A fork can be created by navigating to the [pypsa-kz-data](https://github.com/pypsa-meets-earth/pypsa-kz-data) website and clicking the fork symbol in the upper right corner.

In order to clone both forks to the correct locations on a local machine, the following commands can be used using the local machines shell:
```bash
git clone https://github.com/<user-name>/pypsa-earth
```
`<user-name>` must be replaced with the personal github-username.

After that, one must change to the freshly created pypsa-earth repository.
```bash
cd pypsa-earth/
```
and repeat the cloning, this time for the pypsa-kz-data repository.
```bash
git clone https://github.com/<user-name>/pypsa-kz-data
```
Again, `<user-name>` must be replaced with the personal github-username.

In order to install the pypsa-earth environment, instructions are provided in the pypsa-earth [documentation](https://pypsa-earth.readthedocs.io/en/latest/installation.html), see `Install dependencies` and `Python dependencies`.
After installing the environment, activate it using
```bash
conda activate pypsa-earth
```
Before the whole workflow can be executed, the databundle must be retrieved. This can be done via:
```bash
snakemake -j 1 retrieve_databundle_light
```
This step can optionally be skipped if the `data/` folder with all relevant subfolders already exists.

To adapt the overall workflow for kz, only three further changes are necessary:
1. replace the `pypsa-earth/data/custom_powerplants.csv` with the provided `pypsa-kz-data/data/custom_powerplants.csv`. This can be done using the command:
```bash
cp pypsa-kz-data/data/custom_powerplants.csv data/custom_powerplants.csv
```

2. Open the Snakefile (in `pypsa-earth/`) and navigate to line 25, which should read `configfile: "config.yaml"`. Replace this line with `configfile: "pypsa-kz-data/config_kz_default.yaml"`. To run the default scenario for different weather years (2011, 2013, 2018), add a new line below (line 26) with `configfile: "pypsa-kz-data/config_kz_<year>.yaml"`, where `<year>` represents an existing weather year (2011, 2013, 2018).

The whole workflow can be reproduced by executing
```bash
snakemake -j 1 solve_everything
```

3. Only for scenarios: To run a certain scenario, make sure to update the *2nd* config file, i.e. navigate to line 26 in the Snakefile, which now should read: `configfile: "pypsa-kz-data/config_kz_<year>.yaml"` and replace this with the scenario you want to execute: `configfile: "pypsa-kz-data/config_kz_<year>_discount<p>.yaml"`. `<year>` and `<p>` must be replaced with existing years (2011, 2013, 2018) and discount rates (10 for optimistic scenario, 15 for BAU and 20 for pessimistic scenario). To run the coal exit scenario, the corresponding config file also must be referred in line 26 of the `Snakemake`. For this setting, the line must read `configfile: "pypsa-kz-data/config_kz_<year>_discount<p>_coalexit.yaml"`, where `<year>` and `<p>` are the weather year (2011, 2013, 2018) and discount rates (10, 15, and 20).

Again, the whole workflow can be reproduced by executing the same command as above:
```bash
snakemake -j 1 solve_everything
```

Results are generated and locally saved in `pypsa-earth/results/<scenario_folder>/networks/`.

## Potential errors:
- A rule is killed. In this case, open the `Snakefile` in `pypsa-earth` or open `kz.smk` in `pypsa-kz-data` (depending on the rule which is killed), navigate to the rule that is being killed in the workflow and increase the memory assignment (for example, add a 0 at the end).

- The workflow runs into an error during the `build_powerplants` rule. In this case, try to repeat step 1. of the workflow using the command
```bash
cp pypsa-kz-data/data/custom_powerplants.csv data/custom_powerplants.csv
```

- Unusual error arising from either Snakemake or the `Snakefile` and proving to be challenging to comprehend: Inspect all indentation. Ensure there is no tab spacing; employ only spaces, i.e., ` `. It is probable that the indentations before `configfile: 'pypsa-kz-data/config_kz_<...>.yaml'` are tabs instead of four spaces.

## Comes in handy:
After all cutouts were generated (i.e. the three files `asia-<year>-era5.nc` exist in the folder `pypsa-earth/cutouts/`, where `<year>` is 2011, 2013, and 2018, navigate to `pypsa-earth/pypsa-kz-data`, open the default config file, navigate to line 36, which should read `build_cutout: True`, and set it to `build_cutout: false`. This will save you a lot of time when (re-)runnig scenarios. But remember to set it back to `true` in case one of the cutouts was deleted!

## Acknowledgement
Code development and testing:<td align="center">
    <a href="https://openenergytransition.org/about-us.html#team">
        <b>Open Energy Transition</b>
    </a>
</td>

Model assumptions: <td align="center">
    <a href="https://www.agora-energiewende.de/en/">
        <b>Agora Energiewende</b>
    </a>
</td>
