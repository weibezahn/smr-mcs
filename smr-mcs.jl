##### activate environment #####
using Pkg
Pkg.activate(pwd())
using Statistics

##### load functions #####
@info("Loading functions")
include("functions.jl");

##### load project data #####
@info("loading data")
include("data.jl");

##### further simulation data #####

    # number of Monte Carlo runs
    n = Int64(1e5)

    # wholesale electricity price [USD/MWh], lower and upper bound of rand variable
    electricity_price = [52.22, 95.84]

    # weighted average cost of capital (WACC), lower and upper bound of rand variable
    wacc = [0.04, 0.1]

    # scaling
        # scaling options
        opts_scaling = ["Manufacturer", "Roulstone", "Rothwell", "uniform"]
        # scaling parameter, lower and upper bound of random variable
        scaling = [0.25, 0.85]

    # choose scaling option
    opt_scaling = opts_scaling[2]

##### run simulation #####

# initialize result variables
npv_results = DataFrame();
lcoe_results = DataFrame();

# run simulation for all projects
for p in eachindex(pjs)
    @info("running simulation for", p, name = pjs[p].name)
    # generate random variables
    rand_vars = gen_rand_vars(opt_scaling, n, wacc, electricity_price, pjs[p])
    # run Monte Carlo simulation
    results = investment_simulation(pjs[p], rand_vars)
    # normalize NPV to plant capacity [USD/MW]
    npv_results.res = vec(results[1] / pjs[p].plant_capacity)
    rename!(npv_results,:res => pjs[p].name)
    lcoe_results.res = vec(results[2])
    rename!(lcoe_results,:res => pjs[p].name)
end

CSV.write("mcs-npv_results.csv", npv_results)
CSV.write("mcs-lcoe_results.csv", lcoe_results)

##### sensitivity analysis #####

# initialize results variables
si_npv_results = DataFrame();
si_npv_results.si = ["S", "S", "S", "S", "ST", "ST", "ST", "ST"];
si_npv_results.var = ["wacc", "electricity_price", "loadfactor", "investment", "wacc", "electricity_price", "loadfactor", "investment"];
si_lcoe_results = DataFrame();
si_lcoe_results.si = ["S", "S", "S", "S", "ST", "ST", "ST", "ST"];
si_lcoe_results.var = ["wacc", "electricity_price", "loadfactor", "investment", "wacc", "electricity_price", "loadfactor", "investment"];

# run sensitivity analysis for all projects
for p in eachindex(pjs)
    @info("running sensitivity analysis for", p, name = pjs[p].name)
    si_results = sensitivity_index(opt_scaling, n, wacc, electricity_price, pjs[p])
    si_npv_results.res = vcat(collect(si_results[1]),collect(si_results[2]))
    rename!(si_npv_results,:res => pjs[p].name)
    si_lcoe_results.res = vcat(collect(si_results[3]),collect(si_results[4]))
    rename!(si_lcoe_results,:res => pjs[p].name)
end

CSV.write("si-npv_results.csv", si_npv_results)
CSV.write("si-lcoe_results.csv", si_lcoe_results)