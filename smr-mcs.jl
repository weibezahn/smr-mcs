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

for i in opts_scaling
    investment_simulation(i, n, wacc, electricity_price, pjs[5])
end

result = investment_simulation(opts_scaling[2], n, wacc, electricity_price, pjs[5]);

##### benchmark function runtime #####

using BenchmarkTools

@btime investment_simulation(opts_scaling[2], n, wacc, electricity_price, pjs[5])