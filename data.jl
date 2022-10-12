##### data #####

using CSV, DataFrames

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