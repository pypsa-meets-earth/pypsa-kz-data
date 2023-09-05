# About the Project

<img src="https://github.com/pypsa-meets-earth/pypsa-kz-data/assets/53824825/ca7893de-26e2-47ad-a3e4-d91cd6716652" alt="Open Energy Transition Logo" width="280" height="100" align="right">
<img src="https://github.com/pypsa-meets-earth/pypsa-kz-data/assets/53824825/63bd0250-c54a-4ce1-8df3-eb116baac01b" alt="Agora Energiewende Logo" width="240" height="100">
<br>

Agora Energiewende aims to model the Kazakh power system, incorporating a substantial increase in variable generation, such as solar and wind, surpassing the current official mid-term policy goal of 15% of all renewable energy sources (RES) in generation by 2030. This endeavor has received support from Open Energy Transition on the modeling side.

# pypsa-kz-data
Extra data for Kazakhstan model that will be used as input for PyPSA-Earth.
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

## Modeling adaptations for KZ study

To adapt the overall workflow for kz, only two further changes are necessary.

Firstly, open the Snakefile (in `pypsa-earth/`) and navigate to line [1057-1058](https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/Snakefile#L1057-L1058), which should read
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
cp pypsa-kz-data/config.kz_default.yaml config.yaml
```

Done!

## Running KZ scenarios

To run all scenarios for all considered weather years (2011, 2013, 2018), execute first
```bash
snakemake -j1 prepare_kz_scenarios
```
followed by
```bash
snakemake -j1 run_all_scenarios
```
All results are generated and locally saved in `pypsa-earth/results/<scenario_folder>/networks/`.

## Potential errors:
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

## Comes in handy:
After all cutouts were generated (i.e. the three files `asia-<year>-era5.nc` exist in the folder `pypsa-earth/cutouts/`, where `<year>` is 2011, 2013, and 2018, navigate to `pypsa-earth/pypsa-kz-data`, open the default config file, navigate to line 36, which should read `build_cutout: True`, and set it to `build_cutout: false`. This will save you a lot of time when (re-)runnig scenarios. But remember to set it back to `true` in case one of the cutouts was deleted!
