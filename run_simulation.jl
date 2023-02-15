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

# summary statistics
npv_summary = describe(npv_results, :all)
lcoe_summary = describe(lcoe_results, :all)

# output
CSV.write("$outputpath/mcs-npv_results-$opt_scaling.csv", npv_results);
CSV.write("$outputpath/mcs-npv_summary-$opt_scaling.csv", npv_summary[!,1:8]);
CSV.write("$outputpath/mcs-lcoe_results-$opt_scaling.csv", lcoe_results);
CSV.write("$outputpath/mcs-lcoe_summary-$opt_scaling.csv", lcoe_summary[!,1:8]);

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

#output
CSV.write("$outputpath/si-npv_results-$opt_scaling.csv", si_npv_results);
CSV.write("$outputpath/si-lcoe_results-$opt_scaling.csv", si_lcoe_results);