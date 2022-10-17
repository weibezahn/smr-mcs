##### activate environment #####
using Pkg
Pkg.activate(pwd())
using Statistics

##### load functions #####
include("functions.jl");

##### load project data #####
include("data.jl");

##### further simulation data #####

    # number of Monte Carlo runs
    n = Int64(1e6)

    # wholesale electricity price [USD/MWh], lower and upper bound of rand variable
    electricity_price = [52.22, 95.84]

    # weighted average cost of capital (WACC), lower and upper bound of rand variable
    wacc = [0.04, 0.1]

    # scaling
        # scaling options
        opts_scaling = ["Manufacturer", "Roulstone", "Rothwell", "uniform"]
        # scaling parameter, lower and upper bound of random variable
        scaling = [0.25, 0.85]

##### run simulation #####

# initialize result variables
npv_results = DataFrame();
lcoe_results = DataFrame();

# choose scaling option
opt_scaling = opts_scaling[2]

# run simulation for all projects
for p in eachindex(pjs)
    @info("running simulation for", p, name = pjs[p].name)
    results = investment_simulation(opt_scaling, n, wacc, electricity_price, pjs[p])
    # normalize NPV to plant capacity [USD/MW]
    npv_results.res = vec(results[1] / pjs[p].plant_capacity)
    rename!(npv_results,:res => pjs[p].name)
    lcoe_results.res = vec(results[2])
    rename!(lcoe_results,:res => pjs[p].name)
end

##### benchmark function runtime #####

using BenchmarkTools

@btime investment_simulation(opts_scaling[2], n, wacc, electricity_price, pjs[5])