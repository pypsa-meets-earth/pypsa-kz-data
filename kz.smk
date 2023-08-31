
if config['load_options'].get("rescale_demand", True):
    rule modify_demand:
        input:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            gadm_demand_data="pypsa-kz-data/data/official_demand.csv",
        output:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_d.nc",
        script:
            "scripts/modify_demand.py"


if config["lines"]["modify_lines"].get("limit_line_capacities", True):
    rule modify_lines:
        input:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_d.nc",
        output:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dl.nc",
        script:
            "scripts/modify_lines.py"


if config['load_options'].get("external_loads", True):
    rule add_exports:
        input:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dl.nc",
            imp_exp="pypsa-kz-data/data/electricity_exp-imp.csv",
        output:
            network="networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dle.nc",
        script:
            "scripts/add_exports.py"


def memory(w):
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith("m"):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


rule solve_network_dle:
    params:
        solving=config["solving"],
        augmented_line_connection=config["augmented_line_connection"],
    input:
        "networks/" + RDIR + "elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dle.nc",
    output:
        "results/" + RDIR + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dle.nc",
    log:
        solver=normpath(
            "logs/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dle_solver.log"
        ),
        python="logs/"
        + RDIR
        + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dle_python.log",
        memory="logs/"
        + RDIR
        + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dle_memory.log",
    benchmark:
        (
            "benchmarks/"
            + RDIR
            + "solve_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dle"
        )
    threads: 20
    resources:
        mem=memory,
    shadow:
        "shallow"
    script:
        "../scripts/solve_network.py"


rule solve_everything:
    input:
        expand(
            ["results/" + RDIR + "networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_dle.nc"],
            **config["scenario"]
        ),
