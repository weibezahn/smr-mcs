##### activate environment #####
using Pkg
Pkg.activate(pwd())
using Statistics

inputpath = "_input"
outputpath = "_output"

##### load functions #####
@info("Loading functions")
include("functions.jl");

##### load project data #####
@info("loading data")
include("data.jl");

##### further simulation data #####

    # number of Monte Carlo runs
    n = Int64(1e6);

    # wholesale electricity price [USD/MWh], lower and upper bound of rand variable
    electricity_price = [52.2, 95.8];

    # weighted average cost of capital (WACC), lower and upper bound of rand variable
    wacc = [0.04, 0.1];

    # scaling
        # scaling options
        opts_scaling = ["manufacturer", "roulstone", "rothwell", "uniform"];
        # scaling parameter, lower and upper bound of random variable
        scaling = [0.25, 0.85];

    # choose scaling option
    if @isdefined(par_job) == true
        # read scaling option from job script parameter
        opt_scaling = opts_scaling[par_job];
        @info("using scaling option $opt_scaling")
    else
        # define scaling option locally
        opt_scaling = opts_scaling[1];
    end

##### run simulation #####

@info("running simulation")
include("run_simulation.jl")