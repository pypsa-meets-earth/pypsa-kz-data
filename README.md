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

2. Open the Snakefile (in `pypsa-earth/`) and navigate to line 25, which should read `configfile: "config.yaml"`. Replace this line with `configfile: "pypsa-kz-data/config_kz_default.yaml"`.

The whole workflow can be reproduced by executing
```bash
snakemake -j 1 solve_everything
```

3. Only for scenarios: To run a certain scenario, make sure to update the config file. I.e. navigate to line 25 in the Snakefile, which now should read: `configfile: "pypsa-kz-data/config_kz_default.yaml"` and add a new line below with `configfile: "pypsa-kz-data/config_kz_<year>_discount<p>.yaml"`. `<year>` and `<p>` must be replaced with existing years (2011, 2013, 2018) and discount rates (10 for optmisitic scenario, 15 for BAU and 20 for pessimistic scenario).

Again, the whole workflow can be reproduced by executing the same command as above:
```bash
snakemake -j 1 solve_everything
```

Results are generated and locally saved in `pypsa-earth/results/<scenario_folder>/networks/`.

# Potential errors:
A rule is killed. In this case, open the `Snakefile` in `pypsa-earth` or open `kz.smk` in `pypsa-kz-data` (depending on the rule which is killed), navigate to the rule that is being killed in the workflow and increase the memory assignment (for example, add a 0 at the end).
