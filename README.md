<!--
SPDX-FileCopyrightText:  pypsa-kz-data authors

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# About the Project

<img src="https://github.com/pypsa-meets-earth/pypsa-kz-data/assets/53824825/ca7893de-26e2-47ad-a3e4-d91cd6716652" alt="Open Energy Transition Logo" width="280" height="100" align="right">
<img src="https://github.com/pypsa-meets-earth/pypsa-kz-data/assets/53824825/63bd0250-c54a-4ce1-8df3-eb116baac01b" alt="Agora Energiewende Logo" width="240" height="100">
<br>

Agora Energiewende aims to model the Kazakh power system, incorporating a substantial increase in variable generation, such as solar and wind, surpassing the current official mid-term policy goal of 15% of all renewable energy sources (RES) in generation by 2030. This endeavor has received support from Open Energy Transition on the modeling side.

Used for [2024 study](https://www.agora-energiewende.org/publications/modernising-kazakhstans-coal-dependent-power-sector-through-renewables), "Modernising Kazakhstan’s coal-dependent power sector through renewables - Challenges, solutions and scenarios up to 2030 and beyond"

## Development status: Active and Stable

[![CI-Linux](https://github.com/pypsa-meets-earth/pypsa-kz-data/actions/workflows/ci-linux.yml/badge.svg?branch=main)](https://github.com/pypsa-meets-earth/pypsa-kz-data/actions/workflows/ci-linux.yml)
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

## Setting up the general repositories

The provided workflow builds on [PyPSA-Earth](https://github.com/pypsa-meets-earth/pypsa-earth). Therefore, first, the PyPSA-Earth repository must be forked and the fork should then be cloned. A fork can be created by navigating to the [PyPSA-Earth](https://github.com/pypsa-meets-earth/pypsa-earth) website. By clicking on the fork-symbol in the upper right corner, a fork is created and linked to the specific user. While making a fork, unclick `Copy the main branch only` option to fork all branches and tags of `pypsa-earth` repository.

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
Switch to the stable `v0.4.0` version of `pypsa-earth` that is compatible with `pypsa-kz-data` repository:
```bash
git checkout tags/v0.4.0
```
Repeat the cloning, this time for the pypsa-kz-data repository.
```bash
git clone https://github.com/<user-name>/pypsa-kz-data
```
Again, `<user-name>` must be replaced with the personal github-username.

In order to install the pypsa-earth environment, instructions are provided in the pypsa-earth [documentation](https://pypsa-earth.readthedocs.io/en/latest/installation.html), see `Install dependencies` and `Python dependencies`.
After installing the environment, activate it using
```bash
conda activate pypsa-earth
```

## Modeling adaptations for KZ study

To adapt the overall workflow for kz, only two further changes are necessary.

Firstly, open the Snakefile (in `pypsa-earth/`) and navigate to line [1071-1072](https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/Snakefile#L1057-L1058), which should read
```bash
os.system("snakemake -j all solve_all_networks --rerun-incomplete")
os.system("snakemake -j1 make_statistics --force")
```
and replace these two lines with
```bash
os.system(f"snakemake -j1 networks/{wildcards.scenario_name}/base.nc")
os.system("cp pypsa-kz-data/data/custom_powerplants.csv data/custom_powerplants.csv")
os.system("snakemake -j1 solve_everything --rerun-incomplete")
```

Secondly, copy the default configuration file to the pypsa-earth folder using:
```bash
cp pypsa-kz-data/config.kz_default.yaml config.default.yaml
```
In case you already have a custom config file, make sure to replace it as well, using
```bash
cp pypsa-kz-data/config.kz_default.yaml config.yaml
```
**Note!** Run two aforementioned commands in `pypsa-earth` directory.
You are now all set to run all scenarios!

## Running KZ scenarios

To prepare running all scenarios, execute
```bash
snakemake -j1 prepare_kz_scenarios
```
Optionally, to save time for future runs, you can now set `enable: retrieve_databundle: True` in the `config.yaml` to `False`. If you already have build all cutouts for 2011, 2013 and 2018, you can also set `enable: build_cutout: True` to `False`.
Finally, to run all scenarios, execute
```bash
snakemake -j1 run_all_scenarios
```

After all scenarios have executed successfully, all results are generated and locally saved in `pypsa-earth/results/<scenario_folder>/networks/`.

# Potential errors
- A rule is killed. In this case, open the `Snakefile` in `pypsa-earth` or open `kz.smk` in `pypsa-kz-data` (depending on the rule which is killed), navigate to the rule that is being killed in the workflow and increase the memory assignment (for example, add a 0 at the end).

- The workflow runs into an error during the `build_powerplants` rule. In this case, try to repeat step 1. of the workflow using the command
```bash
cp pypsa-kz-data/data/custom_powerplants.csv data/custom_powerplants.csv
```

- Unusual error arising from either Snakemake or the `Snakefile` and proving to be challenging to comprehend: Inspect all indentation. Ensure there is no tab spacing; employ only spaces, i.e., ` `. It is probable that the indentations before
```bash
os.system(f"snakemake -j1 networks/{wildcards.scenario_name}/base.nc")
os.system("cp pypsa-kz-data/data/custom_powerplants.csv data/custom_powerplants.csv")
os.system("snakemake -j1 solve_everything --rerun-incomplete")
```
are tabs instead of four spaces.

- Missing `data/` folder or some relevant subfolders. This should normally be executed automatically when executing the rule `prepare_kz_scenarios`, however might be missing due to incorrect execution. The databundle can be also retrieved manually via:
```bash
snakemake -j 1 retrieve_databundle_light
```

- The rule `retrieve_databundle_light` always executes with an error. To avoid this, try setting `enable: build_cutout: False` to `True`.

# Comes in handy
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
