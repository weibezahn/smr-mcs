##### julia job file #####
# main file
# to be started as cluster job

par_job = parse(Int, ARGS[1]) # read-in scaling parameter argument

include("smr-mcs.jl")
include("plots.jl")