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
juliaversion = 15 #use 13 for 1.3 and 15 for 1.5

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
nVARSpec    = "2"
nVARMCMCSpec= "1"
specDir   = "$(pwd())/Empirics_Estimation_v2/SpecFiles/"
include(specDir * "/AltVARspec" * nVARSpec * ".jl")
include(specDir * "/AltVARMCMCspec" * nVARMCMCSpec * ".jl")

#-------------------------------------------------------------
# load aggregate data
#-------------------------------------------------------------

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

#-------------------------------------------------------------
# Generate prior
#-------------------------------------------------------------

YY, XX, PHIhat, Shat, SIGMAhat = data_to_YY(data, Tdrop, nlags)
prior = generate_prior(YY, PHIhat, SIGMAhat, n_agg, iwishdf, lambda1, lambda2, lambda3)
T  = size(YY)[1]
ny = size(YY)[2]

#-------------------------------------------------------------
# Run Gibbs Sampler
#-------------------------------------------------------------

if silent == 0
   println("")
   println(" Bayesian estimation of VAR: Gibbs Sampling... ")
   println("")
end

# Initialize matrices for collecting draws from posterior
PHIpdraw     = zeros(nsim, ny*nlags, ny);
SIGMAtrpdraw = zeros(nsim, ny, ny);

let counter = 0
SIGMAj = SIGMAhat;
PHIj   = PHIhat;
time_init_loop = time_ns();

    for jj = 1:nsim

        # VAR estimation with the filtered states (without kronecker structure)
        PHIj, SIGMAj         = bayesianVAR_j_wokron(data, Tdrop, nlags, n_agg, prior, SIGMAj, seedoffset, jj)
        PHIpdraw[jj,:,:]     = PHIj
        SIGMAjtr             = cholesky(SIGMAj)
        SIGMAtrpdraw[jj,:,:] = SIGMAjtr.L

        #output
        counter = counter + 1
        if counter == ncount
           println("")
           println(" Draw number:  $jj")
           println(" Remaining draws:  $(nsim-jj)")
           time_loop=signed(time_ns()-time_init_loop)/1000000000
           println(" Elapsed time = $(time_loop)")
           counter = 0
           time_init_loop = time_ns()
        end

    end

end

PHIpdraw     = PHIpdraw[nburn+1:nsim,:,:];
SIGMAtrpdraw = SIGMAtrpdraw[nburn+1:nsim,:,:];

#-------------------------------------------------------------
# Save draws
#-------------------------------------------------------------
sName    = "AltVAR" * nVARSpec * "_MCMC" * nVARMCMCSpec * "_" * nSpecStr

savedir  = "$(pwd())/results/" * sName *"/";
try mkdir(savedir) catch; end
save(savedir * sName* "_PostDraws.jld", "PHIpdraw", PHIpdraw, "SIGMAtrpdraw", SIGMAtrpdraw)

#-------------------------------------------------------------
# MISC ITEMS
#-------------------------------------------------------------

PHIpmean     = dropdims(mean(PHIpdraw,dims=1),dims=1)
SIGMAtrpmean = dropdims(mean(SIGMAtrpdraw,dims=1),dims=1)
save(savedir * sName* "_PostMeans.jld", "PHIpmean", PHIpmean, "SIGMAtrpmean", SIGMAtrpmean)

# save OLS estiamtes and posterior means as CSV files
CSV.write(savedir * sName * "_PHIols.csv", DataFrame(PHIhat,:auto))
CSV.write(savedir * sName * "_SIGMAols.csv", DataFrame(SIGMAhat,:auto))
CSV.write(savedir * sName * "_PHIpmean.csv", DataFrame(PHIpmean,:auto))
SIGMApmean = SIGMAtrpmean*SIGMAtrpmean'
CSV.write(savedir * sName * "_SIGMApmean.csv", DataFrame(SIGMApmean,:auto))
CSV.write(savedir * sName * "_SIGMAtrpmean.csv", DataFrame(SIGMAtrpmean,:auto))
