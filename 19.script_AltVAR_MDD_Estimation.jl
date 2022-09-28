using DataFrames
using QuadGK
using Distributions
using JLD
using Optim
using ForwardDiff
using CSV
using LinearAlgebra
using SpecialFunctions
using SparseArrays
using Random

#clearconsole()

#-------------------------------------------------------------
# include Functions
#-------------------------------------------------------------
#cd("$(pwd())/Dropbox/Heterogeneity/Software/KS_Simulation/")
readDir = "$(pwd())/Empirics_Estimation_v2/Functions/"
include(readDir *"vech.jl");
include(readDir *"VAR_Procedures.jl");
include(readDir *"Loaddata.jl");

#-------------------------------------------------------------
# choose specification files
#-------------------------------------------------------------
nMDDSpec  = "1"
nMCMCSpec = "1"
specDir   = "$(pwd())/Empirics_Estimation_v2/SpecFiles/"
include(specDir * "/AltMDDspec" * nMDDSpec * ".jl")
include(specDir * "/AltMDDMCMCspec" * nMCMCSpec * ".jl")

#-------------------------------------------------------------
# load aggregate data
#-------------------------------------------------------------
juliaversion = 15

agg_data, period_agg, ~ = loadaggdata(SampleStart,SampleEnd,juliaversion)
n_agg = size(agg_data)[2]

if juliaversion == 13
    percentiles_data     = CSV.read("$(pwd())/Empirics_Estimation_v2/Data/percentiles_data.csv", header = true); # run with version 1.3.
    gini_coef_data       = CSV.read("$(pwd())/Empirics_Estimation_v2/Data/gini_coef_data.csv", header = true); # run with version 1.3.
    probMass_below1_data = CSV.read("$(pwd())/Empirics_Estimation_v2/Data/probMass_below1_data.csv", header = true); # run with version 1.3.
else
    percentiles_data     = CSV.read("$(pwd())/Empirics_Estimation_v2/Data/percentiles_data.csv", DataFrame, header = true); # run with version 1.5.
    gini_coef_data       = CSV.read("$(pwd())/Empirics_Estimation_v2/Data/gini_coef_data.csv", DataFrame, header = true); # run with version 1.5.
    probMass_below1_data = CSV.read("$(pwd())/Empirics_Estimation_v2/Data/probMass_below1_data.csv", DataFrame, header = true); # run with version 1.5.
end

percentiles_data = Matrix(percentiles_data)[2:end,2:end] # include 10, 20, 50, 80, 90th percentiles
vec_percs = [0.1; 0.2; 0.5; 0.8; 0.9]
gini_coef_data = Matrix(gini_coef_data)[2:end,2] # include gini coefficient
probMass_below1_data = Matrix(probMass_below1_data)[2:end,2] # include fraction of individuals earning less than per-capita GDP

#-------------------------------------------------------------
# Combine agg data and density coefficients
#-------------------------------------------------------------
if nSpecStr == "Pctl"
    data    = [ agg_data percentiles_data ]
else
    data    = [ agg_data gini_coef_data probMass_below1_data]
end

sNameSaveDir = "AltMDD" * nMDDSpec * "_MCMC" * nMCMCSpec * "_" * nSpecStr
savedir  = "$(pwd())/results/" * sNameSaveDir *"/";
try mkdir(savedir) catch; end

#-------------------------------------------------------------
# Compute MDD over hyperparameter grid
#-------------------------------------------------------------
silent  = 1
sNameFile  = "AltMDD" * nMDDSpec* "_MCMC" * nMCMCSpec * "_" * nSpecStr
LogFile    = open(savedir * sNameFile * "_MDDLogFile.txt", "w")

for ll = 1:hyper_n

    time_init_loop = time_ns()

    lambda1 = hyper_grid[ll,1]
    lambda2 = hyper_grid[ll,2]
    lambda3 = hyper_grid[ll,3]

    PHIpdraw, SIGMAtrpdraw, YYpred, lnpYY = bayesianVAR_wokron(data, Tdrop, nlags, nsim, nburn,
                                            n_agg, lambda1, lambda2, lambda3, tau, silent, seedoffset, iwishdf);
    MDD_vec[ll] = lnpYY

    print("Step $ll out of $hyper_n \n")
    print("Hyperparameter lambda1 = $lambda1 \n")
    print("Hyperparameter lambda2 = $lambda2 \n")
    print("Hyperparameter lambda3 = $lambda3 \n")
    print("Marginal data density  = $lnpYY \n")
    time_loop=signed(time_ns()-time_init_loop)/1000000000
    print("Elapsed time           = $(time_loop) \n")
    print("================================= \n")

    write(LogFile,"Step $ll out of $hyper_n \n")
    write(LogFile,"Hyperparameter lambda1 = $lambda1 \n")
    write(LogFile,"Hyperparameter lambda2 = $lambda2 \n")
    write(LogFile,"Hyperparameter lambda3 = $lambda3 \n")
    write(LogFile,"Marginal data density  = $lnpYY \n")
    write(LogFile,"Elapsed time           = $(time_loop) \n")
    write(LogFile,"================================= \n")

end

close(LogFile)

#-------------------------------------------------------------
# Combine MDD terms, save output
#-------------------------------------------------------------

lambda_MDD = [hyper_grid MDD_vec]
sNameCols = ["Lam1", "Lam2", "Lam3", "MDD"]
sNameFile = "AltMDD" * nMDDSpec* "_MCMC" * nMCMCSpec * "_" * nSpecStr
CSV.write(savedir * sNameFile * "_lambda_MDD.csv", DataFrame(permutedims(sNameCols),:auto), writeheader=false)
CSV.write(savedir * sNameFile * "_lambda_MDD.csv", DataFrame(lambda_MDD,:auto), append=true);
