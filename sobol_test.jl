##### activate environment #####
using Pkg
Pkg.activate(pwd())

# testing our implementation from
# vs. the implementation in GlobalSensitivityAnalysis.jl

using Distributions
using DataStructures
using LinearAlgebra

using GlobalSensitivityAnalysis

# non-linear function to simulate
f(x,y) = x * exp(-y)

# number of runs
n = 100

# define the data (from GlobalSensitivityAnalysis)
data = SobolData(
    params = OrderedDict(
        :x => Uniform(0,1),
        :y => Uniform(0,1)),
    calc_second_order = true,
    N = 100
)

# generate samples using Sobol sequence (from GlobalSensitivityAnalysis)
samples = GlobalSensitivityAnalysis.sample(data);

# run the model
Y = zeros(size(samples, 1), 1);

for i in 1:size(samples, 1)
    Y[i] = f(samples[i,1], samples[i,2])
end

# extract the output in order to use the same results

function split_output(model_output::AbstractArray{<:Number, S}, N, D, calc_second_order::Bool) where S

    if calc_second_order
        stepsize = 2 * D + 2 
    else
        stepsize = D + 2
    end

    A = model_output[1:stepsize:end]
    B = model_output[stepsize:stepsize:end]
    
    #preallocate
    AB = Array{Float64}(undef, N, D)
    if calc_second_order
        BA = Array{Float64}(undef, N, D)
    else
        BA = nothing
    end

    for i in 1:D
        AB[:, i] = model_output[i+1:stepsize:end, :] 
        if calc_second_order
            BA[:, i] = model_output[i + D + 1:stepsize:end, :]
        end
    end

    return A, B, AB, BA
end

# matrix f(A)
A = split_output(Y, data.N, size(samples, 2), true)[1];
# matrix f(B)
B = split_output(Y, data.N, size(samples, 2), true)[2];
# matrix f(A_B), i.e. model results from input matrix A with one random parameter replaced with values from matrix B for each respective column
AB = split_output(Y, data.N, size(samples, 2), true)[3];
# matrix f(B_A), i.e. model results from input matrix B with one random parameter replaced with values from matrix A for each respective column
BA = split_output(Y, data.N, size(samples, 2), true)[4];

# first-order sensitivity index from GlobalSensitivityAnalysis following Saltelli et al. (2010)
# mean(B .* (AB .- A), dims = 1) ./ var(vcat(A, B), dims = 1, corrected = false)
# total-effect sensitivity index from GlobalSensitivityAnalysis following Saltelli et al. (2010)
# 0.5 * mean((A .- AB).^2, dims = 1) ./ var(vcat(A, B), dims = 1, corrected = false)

# perform Sobol Analysis
GlobalSensitivityAnalysis.analyze(data, Y)

# first-order sensitivity index from our implementation following Saltelli (2008)
sensitivity_index(sensi_res_A, sensi_res_C) = (((sensi_res_A ⋅ sensi_res_C) / length(sensi_res_A)) - mean(sensi_res_A)^2) / (((sensi_res_A ⋅ sensi_res_A) / length(sensi_res_A)) - mean(sensi_res_A)^2)

[sensitivity_index(A,BA[:,1]),sensitivity_index(A,BA[:,2])]
# second try: since the GSA implementation uses values from B inserted to A, this is switched around here (A and B exchanged)
[sensitivity_index(B,AB[:,1]),sensitivity_index(B,AB[:,2])]

# total-effect sensitivity index from our implementation following Saltelli (2008)
sensitivity_index(sensi_res_A, sensi_res_B, sensi_res_C) = 1 - (((sensi_res_B ⋅ sensi_res_C) / length(sensi_res_A)) - mean(sensi_res_A)^2) / (((sensi_res_A ⋅ sensi_res_A) / length(sensi_res_A)) - mean(sensi_res_A)^2)

[sensitivity_index(A,B,BA[:,1]),sensitivity_index(A,B,BA[:,2])]
# second try: since the GSA implementation uses values from B inserted to A, this is switched around here (A and B exchanged)
[sensitivity_index(B,A,AB[:,1]),sensitivity_index(B,A,AB[:,2])]

sensitivity_index_2(sensi_res_A, sensi_res_B, sensi_res_C) = (((sensi_res_A ⋅ sensi_res_C) / length(sensi_res_A)) - mean(vcat(sensi_res_A,sensi_res_B))^2) / (((sensi_res_A ⋅ sensi_res_A) / length(sensi_res_A)) - mean(vcat(sensi_res_A,sensi_res_B))^2)
[sensitivity_index_2(A,B,BA[:,1]),sensitivity_index_2(A,B,BA[:,2])]