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

##### defining functions #####

"""
function to generate random variables
"""
function gen_rand_var(opt_scaling::String, n::Int64, wacc::Vector, electricity_price::Vector, pj::project)

    @info "generating random variables"

    # total time of reactor project, i.e., construction and operating time
    total_time = pj.time[1] + pj.time[2]

    # generation of uniformly distributed random non-project specific variables
        rand_wacc = wacc[1] .+ (wacc[2]-wacc[1]) * rand(n)
        rand_electricity_price = electricity_price[1] .+ (electricity_price[2]-electricity_price[1]) * rand(n,total_time)
        rand_loadfactor = pj.loadfactor[1] .+ (pj.loadfactor[2]-pj.loadfactor[1]) * rand(n,total_time)
    
    # generation of uniformly distributed random project specific variables
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

    # output
    return(rand_wacc,rand_electricity_price,rand_loadfactor,rand_investment)

end

"""
function for Monte Carlo runs
"""
function mc_run(n::Int64, pj::project, rand_var)

    @info "running Monte Carlo simulation"

    # project data
        # total time of reactor project, i.e., construction and operating time
        total_time = pj.time[1] + pj.time[2]
        # O&M costs
            # fixed O&M costs [USD/year]
            operating_cost_fix = pj.plant_capacity * pj.operating_cost[1]
            # variable O&M costs [USD/MWh]
            operating_cost_variable = pj.operating_cost[2] + pj.operating_cost[3]

    # initialize variables
        # random variables
        rand_wacc = rand_var[1]
        rand_electricity_price = rand_var[2]
        rand_loadfactor = rand_var[3]
        rand_investment = rand_var[4]
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
                # cash_out[:,t] = pj.plant_capacity * rand_investment ./ pj.time[1] .* (1+idc)
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

    # output
    return(disc_cash_out,disc_cash_net,disc_electricity)

end

"""
function to calculate NPV and LCOE
"""
function npv_lcoe(disc_res)
 
    @info "calculating NPV and LCOE"

    disc_cash_out = disc_res[1]
    disc_cash_net = disc_res[2]
    disc_electricity = disc_res[3]

    npv = sum(disc_cash_net,dims=2)
    lcoe = sum(disc_cash_out,dims=2) ./ sum(disc_electricity,dims=2)

    return(npv,lcoe)

end

"""
function for Monte Carlo investment simulation
"""
function investment_simulation(opt_scaling::String, n::Int64, wacc::Vector, electricity_price::Vector, pj::project)

    # generate random variables
    rand_var = gen_rand_var(opt_scaling, n, wacc, electricity_price, pj)

    # interest during construction factor calculation [Rothwell (2016): Economics of Nuclear Power, Eq. (3.3.8)]
        # idc = .5 * rand_wacc * pj.time[1] + 1/6 * (rand_wacc .^ 2) * (pj.time[1]^2)
    
    # run the Monte Carlo simulation
    disc_res= mc_run(n, pj, rand_var)

    # NPV and LCOE calculation
    res = npv_lcoe(disc_res)

    # output
    @info "simulation results" pj.name pj.type NPV=mean(res[1]) LCOE=mean(res[2])
    return(res)

end