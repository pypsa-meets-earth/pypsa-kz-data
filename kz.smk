
if config['load_options'].get("rescale_demand", True):
    rule modify_demand:
        input:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            gadm_demand_data="pypsa-kz-data/data/official_demand.csv",
        output:
            demand_dummy="logs/" + RDIR + "dummy_log/modify_demand_{simpl}_{clusters}_ec_l{ll}_{opts}.txt",
        priority: 3
        script:
            "scripts/modify_demand.py"


if config["lines"]["modify_lines"].get("limit_line_capacities", True):
    rule modify_lines:
        input:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
        output:
            lines_dummy="logs/" + RDIR + "dummy_log/modify_lines_{simpl}_{clusters}_ec_l{ll}_{opts}.txt",
        priority: 2
        script:
            "scripts/modify_lines.py"


if config['load_options'].get("external_loads", True):
    rule add_exports:
        input:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            imp_exp="pypsa-kz-data/data/electricity_exp-imp.csv",
        output:
            exports_dummy="logs/" + RDIR + "dummy_log/add_exports_{simpl}_{clusters}_ec_l{ll}_{opts}.txt",
        priority: 1
        script:
            "scripts/add_exports.py"


rule solve_everything:
    input:
        expand(
            ["results/" + RDIR + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
             "logs/" + RDIR + "dummy_log/add_exports_{simpl}_{clusters}_ec_l{ll}_{opts}.txt",
             "logs/" + RDIR + "dummy_log/modify_lines_{simpl}_{clusters}_ec_l{ll}_{opts}.txt",
             "logs/" + RDIR + "dummy_log/modify_demand_{simpl}_{clusters}_ec_l{ll}_{opts}.txt"],
            **config["scenario"]
        ),