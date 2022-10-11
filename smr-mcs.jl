##### activate environment #####
using Pkg
using CSV, DataFrames
using BenchmarkTools
Pkg.activate(pwd())

##### defining structs #####

"""
struct for investment project data
"""
mutable struct project
    name::String                    # investment concept
    type::String                    # investment type
    investment::Float64             # investment estimate by manufacturer [USD/MW]
    plant_capacity::Float64         # plant capacity [MW]
    learning_factor::Float64        # learning factor
    time::AbstractVector            # project time [years] (construction time, plant lifetime)
    loadfactor::AbstractVector      # load factor, lower, and upper bound of rand variable
    operating_cost::AbstractVector  # O&M cost (O&M fix cost [USD/MW], O&M variable cost [USD/MWh], fuel cost [USD/MWh])
    reference_pj::AbstractVector    # reference reactor (investment costs [USD/MW], plant capacity [MW])
end

##### loading project data #####

pjs_dat = CSV.File("project_data.csv") |> DataFrame # read project data from CSV into a dataframe
pjs = []                                            # initialize empty projects vector

# populate array with projects using project data input
for i in 1:nrow(pjs_dat)
    push!(pjs,project(
        pjs_dat.name[i],
        pjs_dat.type[i],
        pjs_dat.investment[i],
        pjs_dat.plant_capacity[i],
        pjs_dat.learning_factor[i],
        [pjs_dat.construction_time[i],pjs_dat.operating_time[i]],
        [pjs_dat.loadfactor_lower[i],pjs_dat.loadfactor_upper[i]],
        [pjs_dat.operating_cost_fix[i],pjs_dat.operating_cost_variable[i],pjs_dat.operating_cost_fuel[i]],
        [pjs_dat.reference_pj_investment[i];pjs_dat.reference_pj_capacity[i]]
        ))
end

##### functions #####

"""
function for Monte Carlo investment simulation
"""
function investment_simulation(opt_scaling::String, n::Int64, wacc::Vector, electricity_price::Vector, pj::project)

    # project data
        # scaling case distinction
            # Note that the scaling parameter here are converted such that Rothwell and Roulstone coincide.
            if opt_scaling == "Manufacturer"
                @info("using manufacturer estimates")
                # deterministic investment cost based on manufacturer estimates [USD]
                rand_investment = pj.investment * pj.plant_capacity * ones(n) * (1-pj.learning_factor)
            elseif opt_scaling == "Roulstone"
                @info("using Roulstone scaling")
                # random investment cost based on Roulstone [USD]
                rand_scaling = scaling[1] .+ (scaling[2] - scaling[1]) .* rand(n,1)
                rand_investment = pj.reference_pj[1] * pj.reference_pj[2] * (1-pj.learning_factor) * (pj.plant_capacity/pj.reference_pj[2]) .^ (rand_scaling)
            elseif opt_scaling == "Rothwell"
                @info("using Rothwell scaling")
                # random investment cost based on Rothwell [USD]
                rand_scaling = (2^(scaling[1]-1) + (2^(scaling[2]-1) - 2^(scaling[1]-1))) .* rand(n,1)
                rand_investment = pj.reference_pj[1] * pj.reference_pj[2] * (1-pj.learning_factor) * (pj.plant_capacity/pj.reference_pj[2]) .^ (1 .+ log.(rand_scaling) ./ log(2))
            elseif opt_scaling == "uniform"
                @info("using uniform scaling")
                # random investment cost uniform [USD]
                investment_scaled = pj.reference_pj[1] * pj.reference_pj[2] * (pj.plant_capacity/pj.reference_pj[2]) .^ (scaling)
                rand_investment = investment_scaled[1] .+ (investment_scaled[2]-investment_scaled[1]) .* rand(n,1) .* (1-pj.learning_factor)
            else
                @error("Option for the scaling method is unknown.")
            end
        # total time of reactor project, i.e., construction and operating time
        total_time = pj.time[1] + pj.time[2]
        # O&M costs
            # fixed O&M costs [USD/year]
            operating_cost_fix = pj.plant_capacity * pj.operating_cost[1]
            # variable O&M costs [USD/MWh]
            operating_cost_variable = pj.operating_cost[2] + pj.operating_cost[3]

    # generation of uniformly distributed random variables
        rand_electricity_price = electricity_price[1] .+ (electricity_price[2]-electricity_price[1]) * rand(n,total_time)
        rand_wacc = wacc[1] .+ (wacc[2]-wacc[1]) * rand(n)
        rand_loadfactor = pj.loadfactor[1] .+ (pj.loadfactor[2]-pj.loadfactor[1]) * rand(n,total_time)

    # interest during construction factor calculation [Rothwell (2016): Economics of Nuclear Power, Eq. (3.3.8)]
        # idc = .5 * rand_wacc * construction_time + 1/6 * (rand_wacc .^ 2) * (construction_time^2)

    # initialize variables
        # cash inflow
        cash_in = zeros(Float64, n, total_time)
        # cash outflow
        cash_out = zeros(Float64, n, total_time)
        # cashflow = inflow - outflow
        cash_net = zeros(Float64, n, total_time)
        # discounted cash outflow
        disc_cash_out = zeros(Float64, n, total_time)
        # NPV in period t
        disc_cash_net = zeros(Float64, n, total_time)
        # generated electricity
        electricity = zeros(Float64, n, total_time)
        # discounted generated electricity
        disc_electricity = zeros(Float64, n, total_time)

    # simulation loop
        for t in 1:total_time
            if t <= pj.time[1]
                # outflow(:,t) = plant_capacity*random_investment./construction_time.*(1+idc);
                cash_out[:,t] = rand_investment ./ pj.time[1]
                cash_net[:,t] = cash_in[:,t] - cash_out[:,t]
            else
                # amount of electricity produced
                electricity[:,t] = pj.plant_capacity .* rand_loadfactor[:,t] .* 8760
                # cash inflow from electricity sales
                cash_in[:,t] = rand_electricity_price[:,t] .* electricity[:,t]
                # cash outflow O&M costs
                cash_out[:,t] = operating_cost_fix .+ operating_cost_variable * electricity[:,t]
                cash_net[:,t] = cash_in[:,t] - cash_out[:,t]
                disc_electricity[:,t] = electricity[:,t] ./ ((1 .+ rand_wacc[:]) .^ (t-1))
            end
            disc_cash_out[:,t] = cash_out[:,t] ./ ((1 .+ rand_wacc[:]) .^ (t-1))
            disc_cash_net[:,t] = cash_net[:,t] ./ ((1 .+ rand_wacc[:]) .^ (t-1))
        end

    # NPV and LCOE calculation
    npv = sum(disc_cash_net[:,:])
    lcoe = sum(disc_cash_out[:,:]) / sum(disc_electricity[:,:])

    # output
    @info("simulation results",pj.name,pj.type,npv,lcoe)
    return(npv,lcoe)
end

##### input data #####

# scaling option
opts_scaling = ["Manufacturer", "Roulstone", "Rothwell", "uniform"]

# number of Monte Carlo runs
n = Int64(1e6)

# wholesale electricity price [USD/MWh], lower and upper bound of rand variable
electricity_price = [52.22, 95.84]

# weighted average cost of capital (WACC), lower and upper bound of rand variable
wacc = [0.04, 0.1]

# scaling parameter, lower and upper bound of random variable
scaling = [0.25, 0.85]

##### run simulation #####

for i in opts_scaling
    investment_simulation(i, n, wacc, electricity_price, pjs[5])
end

##### benchmark function runtime #####

@btime investment_simulation(opts_scaling[2], n, wacc, electricity_price, pjs[5])
