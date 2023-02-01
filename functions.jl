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
The gen_rand_vars function generates random variables for a given investment project, pj, based on a specified scaling option, opt_scaling, the number of simulations to run, n, a range of weighted average cost of capital (WACC) values, wacc, and a range of electricity price values, electricity_price.
The total time of the reactor project, which is the sum of the construction time and operating time, is calculated by adding the elements at index 1 and 2 of the pj.time vector.
Next, the function generates uniformly distributed random variables for the WACC, electricity price, and load factor using the rand function. The range of values for each variable is determined by the input ranges of wacc and electricity_price, and the pj.loadfactor attribute.
The function then generates random project-specific variables based on the opt_scaling input. The function has four different branches of execution based on the value of opt_scaling.
    If opt_scaling is "Manufacturer", the function uses manufacturer estimates to calculate a deterministic investment cost, rand_investment, which is the product of pj.investment, pj.plant_capacity, and ones(n) (where ones(n) is an array of n ones) multiplied by (1-pj.learning_factor)
    If opt_scaling is "Roulstone", the function uses Roulstone scaling to calculate a random investment cost, rand_investment, which is the product of pj.reference_pj[1], pj.reference_pj[2], (1-pj.learning_factor), (pj.plant_capacity/pj.reference_pj[2]) raised to the power of a random scaling factor, rand_scaling, within the range specified in the scaling input.
    If opt_scaling is "Rothwell", the function uses Rothwell scaling to calculate a random investment cost, rand_investment, which is the product of pj.reference_pj[1], pj.reference_pj[2], (1-pj.learning_factor), (pj.plant_capacity/pj.reference_pj[2]) raised to the power of 1 + logarithm base 2 of rand_scaling which is calculated by the range specified in the scaling input
    If opt_scaling is "uniform", the function uses a uniform scaling to calculate a random investment cost, rand_investment, which is a random value within the range specified in the scaling input
    If the opt_scaling is not one of the above, the function print an error message, "Option for the scaling method is unknown."
Finally, the function returns the generated random variables in the form of a named tuple with four fields: wacc, electricity_price, loadfactor, and investment and their corresponding values.
"""
function gen_rand_vars(opt_scaling::String, n::Int64, wacc::Vector, electricity_price::Vector, pj::project)

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
            if opt_scaling == "manufacturer"
                @info("using manufacturer estimates")
                # deterministic investment cost based on manufacturer estimates [USD]
                rand_investment = pj.investment * pj.plant_capacity * ones(n) * (1-pj.learning_factor)
            elseif opt_scaling == "roulstone"
                @info("using Roulstone scaling")
                # random investment cost based on Roulstone [USD]
                rand_scaling = scaling[1] .+ (scaling[2] - scaling[1]) .* rand(n,1)
                rand_investment = pj.reference_pj[1] * pj.reference_pj[2] * (1-pj.learning_factor) * (pj.plant_capacity/pj.reference_pj[2]) .^ (rand_scaling)
            elseif opt_scaling == "rothwell"
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
    return(wacc = rand_wacc, electricity_price = rand_electricity_price, loadfactor = rand_loadfactor, investment = rand_investment)

end

"""
The mc_run function performs a Monte Carlo simulation for an investment project, pj, based on a specified number of simulations to run, n, and a set of random variables, rand_vars.
The function first initializes the total time of the reactor project, which is the sum of the construction time and operating time, and the fixed and variable operating and maintenance (O&M) costs.
It then initializes variables for the random variables, cash inflow, cash outflow, cashflow, discounted cash outflow, NPV in period t, generated electricity, and discounted generated electricity.
The function then enters a simulation loop that iterates over the total time of the project. For each time step, the function performs different calculations depending on whether the time step is before or after the construction time.
    If the time step is before the construction time, the function calculates the cash outflow as the investment cost divided by the construction time, and cashflow as the difference between cash inflow and cash outflow.
    If the time step is after the construction time, the function calculates the amount of electricity produced, cash inflow from electricity sales, cash outflow for O&M costs, cashflow, discounted cash outflow, discounted electricity.
Finally, the function returns the results of the simulation in the form of a named tuple with three fields: disc_cash_out, disc_cash_net, and disc_electricity and their corresponding values.
"""
function mc_run(n::Int64, pj::project, rand_vars)

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
        rand_wacc = rand_vars.wacc
        rand_electricity_price = rand_vars.electricity_price
        rand_loadfactor = rand_vars.loadfactor
        rand_investment = rand_vars.investment
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
    return(disc_cash_out = disc_cash_out, disc_cash_net = disc_cash_net, disc_electricity = disc_electricity)

end

"""
The npv_lcoe function calculates the net present value (NPV) and levelized cost of electricity (LCOE) for a given input, disc_res. The input, disc_res, is assumed to be an object with three attributes: disc_cash_out, disc_cash_net, and disc_electricity.
Next, the function creates three local variables disc_cash_out, disc_cash_net, and disc_electricity which are assigned the values of the corresponding attributes of the input disc_res.
The NPV is then calculated by taking the sum of the disc_cash_net variable along the second dimension using the sum function. Similarly, the LCOE is calculated by taking the ratio of the sum of the disc_cash_out variable along the second dimension to the sum of the disc_electricity variable along the second dimension.
Finally, the function returns a named tuple with two fields, npv and lcoe with the calculated values.
"""
function npv_lcoe(disc_res)
 
    @info "calculating NPV and LCOE"

    disc_cash_out = disc_res.disc_cash_out
    disc_cash_net = disc_res.disc_cash_net
    disc_electricity = disc_res.disc_electricity

    npv = sum(disc_cash_net,dims=2)
    lcoe = sum(disc_cash_out,dims=2) ./ sum(disc_electricity,dims=2)

    return(npv = npv, lcoe = lcoe)

end

"""
The investment_simulation function simulates an investment project, given an instance of the project type pj and a set of random variables rand_vars.
The function runs a Monte Carlo simulation by calling the mc_run function and passing in the number of simulations to run, n, the project pj, and the set of random variables rand_vars. The function saves the results of the Monte Carlo simulation in the local variable disc_res.
The function then calculates the net present value (NPV) and the levelized cost of electricity (LCOE) by calling the npv_lcoe function and passing in the disc_res variable as input. The results of this calculation are saved in the local variable res.
The function then returns the results of the simulation in the form of the res variable.
"""
function investment_simulation(pj::project, rand_vars)
   
    # run the Monte Carlo simulation
    disc_res = mc_run(n, pj, rand_vars)

    # NPV and LCOE calculation
    res = npv_lcoe(disc_res)

    # output
    @info "simulation results" pj.name pj.type NPV = mean(res[1]) LCOE = mean(res[2])
    return(res)

end

# first-order sensitivity index
function si_first_order(A,B,AB)
    
    si = mean(B .* (AB .- A)) / var(vcat(A, B), corrected = false)

    if isapprox(si, -0.0, atol = 1e-2)
        return 0.0
    else
        return si
    end

end

# total-effect sensitivity index
function si_total_order(A,B,AB)
    
    si = 0.5 * mean((A .- AB).^2) / var(vcat(A, B), corrected = false)

    if isapprox(si, -0.0, atol = 1e-2)
        return 0.0
    else
        return si
    end

end

function sensitivity_index(opt_scaling::String, n::Int64, wacc::Vector, electricity_price::Vector, pj::project)

    # generate random variable matrices A and B
    @info "generating matrix A"
    rand_vars_A = gen_rand_vars(opt_scaling, n, wacc, electricity_price, pj);
    @info "generating matrix B"
    rand_vars_B = gen_rand_vars(opt_scaling, n, wacc, electricity_price, pj);

    # build random variable matrices AB from A and B for each random variable
    @info "building matrices AB"
    rand_vars_AB1 = (wacc = rand_vars_B.wacc, electricity_price = rand_vars_A.electricity_price, loadfactor = rand_vars_A.loadfactor, investment = rand_vars_A.investment);
    rand_vars_AB2 = (wacc = rand_vars_A.wacc, electricity_price = rand_vars_B.electricity_price, loadfactor = rand_vars_A.loadfactor, investment = rand_vars_A.investment);
    rand_vars_AB3 = (wacc = rand_vars_A.wacc, electricity_price = rand_vars_A.electricity_price, loadfactor = rand_vars_B.loadfactor, investment = rand_vars_A.investment);
    rand_vars_AB4 = (wacc = rand_vars_A.wacc, electricity_price = rand_vars_A.electricity_price, loadfactor = rand_vars_A.loadfactor, investment = rand_vars_B.investment);

    # run Monte Carlo simulations for A, B, and AB
    @info "matrix A:"
    sensi_res_A = investment_simulation(pj, rand_vars_A);
    @info "matrix B:"
    sensi_res_B = investment_simulation(pj, rand_vars_B);
    @info "matrix AB1"
    sensi_res_AB1 = investment_simulation(pj, rand_vars_AB1);
    @info "matrix AB2"
    sensi_res_AB2 = investment_simulation(pj, rand_vars_AB2);
    @info "matrix AB3"
    sensi_res_AB3 = investment_simulation(pj, rand_vars_AB3);
    @info "matrix AB4"
    sensi_res_AB4 = investment_simulation(pj, rand_vars_AB4);

    # sensitivity indices for NPV
    s_npv = (
        wacc = si_first_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB1[1]),
        electricity_price = si_first_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB2[1]),
        loadfactor = si_first_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB3[1]),
        investment = si_first_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB4[1]),
        )
    st_npv = (
        wacc = si_total_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB1[1]),
        electricity_price = si_total_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB2[1]),
        loadfactor = si_total_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB3[1]),
        investment = si_total_order(sensi_res_A[1], sensi_res_B[1], sensi_res_AB4[1])
        )

    # sensitivity indices for LCOE
    s_lcoe = (
        wacc = si_first_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB1[2]),
        electricity_price = si_first_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB2[2]),
        loadfactor = si_first_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB3[2]),
        investment = si_first_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB4[2]),
        )
    st_lcoe = (
        wacc = si_total_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB1[2]),
        electricity_price = si_total_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB2[2]),
        loadfactor = si_total_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB3[2]),
        investment = si_total_order(sensi_res_A[2], sensi_res_B[2], sensi_res_AB4[2])
        )
    
    # output
    @info "sensitivity results" pj.name pj.type S_NPV = s_npv ST_NPV = st_npv S_LCOE = s_lcoe ST_LCOE = st_lcoe
    return(s_npv = s_npv, st_npv = st_npv, s_lcoe = s_lcoe, st_lcoe = st_lcoe)

end

function gen_scaled_investment(scaling::Vector, pj::project)
    
    # generation of project specific scaled investment cost ranges
    # note that the scaling parameter here are converted such that Rothwell and Roulstone coincide.
        # deterministic investment cost based on manufacturer estimates [USD/MW]
        scaled_investment = pj.investment
        # scaled investment cost based on Roulstone [USD/MW]
        scaled_investment = vcat(scaled_investment, pj.reference_pj[1] * pj.reference_pj[2] * (1-pj.learning_factor) * (pj.plant_capacity/pj.reference_pj[2]) .^ scaling / pj.plant_capacity)
        # scaled investment cost based on Rothwell [USD/MW]
        scaled_investment = vcat(scaled_investment, pj.reference_pj[1] * pj.reference_pj[2] * (1-pj.learning_factor) * (pj.plant_capacity/pj.reference_pj[2]) .^ (1 .+ log.(scaling) ./ log(2)) / pj.plant_capacity)

    # output
    return(scaled_investment = scaled_investment)

end