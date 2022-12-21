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

##### benchmark function runtime #####

using BenchmarkTools

@btime investment_simulation(opts_scaling[2], n, wacc, electricity_price, pjs[5])

##### sensitivity analysis #####
##### testing stage        #####

using LinearAlgebra

    # first-order sensitivity index
    sensitivity_index(sensi_res_A, sensi_res_C) = (sensi_res_A ⋅ sensi_res_C - mean(sensi_res_A)^2) / (sensi_res_A ⋅ sensi_res_A - mean(sensi_res_A)^2)
    # total-effect sensitivity index
    sensitivity_index(sensi_res_A, sensi_res_B, sensi_res_C) = 1 - (sensi_res_B ⋅ sensi_res_C - mean(sensi_res_A)^2) / (sensi_res_A ⋅ sensi_res_A - mean(sensi_res_A)^2)

# testing for Nuscale and Roulstone

    # generate random variable matrices A, B, C
    rand_vars_A = gen_rand_vars(opt_scaling, n, wacc, electricity_price, pjs[5]);
    rand_vars_B = gen_rand_vars(opt_scaling, n, wacc, electricity_price, pjs[5]);
    rand_vars_C1 = (wacc = rand_vars_A.wacc, electricity_price = rand_vars_B.electricity_price, loadfactor = rand_vars_B.loadfactor, investment = rand_vars_B.investment);
    rand_vars_C2 = (wacc = rand_vars_B.wacc, electricity_price = rand_vars_A.electricity_price, loadfactor = rand_vars_B.loadfactor, investment = rand_vars_B.investment);
    rand_vars_C3 = (wacc = rand_vars_B.wacc, electricity_price = rand_vars_B.electricity_price, loadfactor = rand_vars_A.loadfactor, investment = rand_vars_B.investment);
    rand_vars_C4 = (wacc = rand_vars_B.wacc, electricity_price = rand_vars_B.electricity_price, loadfactor = rand_vars_B.loadfactor, investment = rand_vars_A.investment);

    # run Monte Carlo simulations for A, B, C
    sensi_res_A = investment_simulation(pjs[5], rand_vars_A);
    sensi_res_B = investment_simulation(pjs[5], rand_vars_B);
    sensi_res_C1 = investment_simulation(pjs[5], rand_vars_C1);
    sensi_res_C2 = investment_simulation(pjs[5], rand_vars_C2);
    sensi_res_C3 = investment_simulation(pjs[5], rand_vars_C3);
    sensi_res_C4 = investment_simulation(pjs[5], rand_vars_C4);

    # sensitivities for NPV
    S_NPV = (
        wacc = sensitivity_index(sensi_res_A[1], sensi_res_C1[1]),
        electricity_price = sensitivity_index(sensi_res_A[1], sensi_res_C2[1]),
        loadfactor = sensitivity_index(sensi_res_A[1], sensi_res_C3[1]),
        investment = sensitivity_index(sensi_res_A[1], sensi_res_C4[1])
        )
    ST_NPV = (
        wacc = sensitivity_index(sensi_res_A[1], sensi_res_B[1], sensi_res_C1[1]),
        electricity_price = sensitivity_index(sensi_res_A[1], sensi_res_B[1], sensi_res_C2[1]),
        loadfactor = sensitivity_index(sensi_res_A[1], sensi_res_B[1], sensi_res_C3[1]),
        investment = sensitivity_index(sensi_res_A[1], sensi_res_B[1], sensi_res_C4[1])
        )

    # sensitivity for LCOE
    S_LCOE = (
        wacc = sensitivity_index(sensi_res_A[2], sensi_res_C1[2]),
        electricity_price = sensitivity_index(sensi_res_A[2], sensi_res_C2[2]),
        loadfactor = sensitivity_index(sensi_res_A[2], sensi_res_C3[2]),
        investment = sensitivity_index(sensi_res_A[2], sensi_res_C4[2])
        )
    ST_LCOE = (
        wacc = sensitivity_index(sensi_res_A[2], sensi_res_B[2], sensi_res_C1[2]),
        electricity_price = sensitivity_index(sensi_res_A[2], sensi_res_B[2], sensi_res_C2[2]),
        loadfactor = sensitivity_index(sensi_res_A[2], sensi_res_B[2], sensi_res_C3[2]),
        investment = sensitivity_index(sensi_res_A[2], sensi_res_B[2], sensi_res_C4[2])
        )