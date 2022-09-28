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

clearconsole()
juliaversion = 15 #use 13 for 1.3 and 15 for 1.5

#-------------------------------------------------------------
# include Functions
#-------------------------------------------------------------
#cd("$(pwd())/Dropbox/Heterogeneity/Software/KS_Simulation/")
readDir = "$(pwd())/Functions/"
include(readDir *"vech.jl");
include(readDir *"VAR_Procedures.jl");
include(readDir *"Loaddata.jl");


#-------------------------------------------------------------
# choose specification files
#-------------------------------------------------------------
nfVARSpec = "10tc"
nVARSpec  = "1"
nVARMCMCSpec= "1"
specDir   = "$(pwd())/SpecFiles/"
include(specDir * "/fVARspec" * nfVARSpec * ".jl")
include(specDir * "/VARspec" * nVARSpec * ".jl")
include(specDir * "/VARMCMCspec" * nVARMCMCSpec * ".jl")


#-------------------------------------------------------------
# load aggregate data
#-------------------------------------------------------------

agg_data, period_agg, ~ = loadaggdata(SampleStart,SampleEnd,juliaversion)
n_agg = size(agg_data)[2]

#-------------------------------------------------------------
# Load coefficients from density estimation
#-------------------------------------------------------------
sNameLoadDir = "fVAR" * nfVARSpec
loaddir  = "$(pwd())/results/" * sNameLoadDir *"/";

PhatDensCoef_factor, MDD_term1, VinvLam_all, period_Dens, ~, ~, ~ = loaddensdata(SampleStart,SampleEnd,K,nfVARSpec,juliaversion)
n_cross = size(PhatDensCoef_factor)[2]

#-------------------------------------------------------------
# Combine agg data and density coefficients
#-------------------------------------------------------------

data    = [ agg_data PhatDensCoef_factor ]

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
   println("Bayesian estimation of VAR: Gibbs Sampling... ")
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
           println("Draw number:  $jj")
           println("Remaining draws:  $(nsim-jj)")
           time_loop=signed(time_ns()-time_init_loop)/1000000000
           println("Elapsed time = $(time_loop)")
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

sName    = "fVAR" * nfVARSpec * "_VAR" * nVARSpec*"_MCMC" * nVARMCMCSpec
savedir  = "$(pwd())/results/" * sName *"/";
try mkdir(savedir) catch; end
save(savedir * sName* "_PostDraws.jld", "PHIpdraw", PHIpdraw, "SIGMAtrpdraw", SIGMAtrpdraw)
#save(savedir * sName* "_StateDraws.jld", "statepdraw", statepdraw)

#-------------------------------------------------------------
# MISC ITEMS
#-------------------------------------------------------------

PHIpmean     = dropdims(mean(PHIpdraw,dims=1),dims=1)
SIGMAtrpmean = dropdims(mean(SIGMAtrpdraw,dims=1),dims=1)
#statepmean   = dropdims(mean(statepdraw,dims=1),dims=1)
save(savedir * sName* "_PostMeans.jld", "PHIpmean", PHIpmean, "SIGMAtrpmean", SIGMAtrpmean)
#save(savedir * sName* "_StateMeans.jld", "statepmean", statepmean)

# save OLS estiamtes and posterior means as CSV files
CSV.write(savedir * sName * "_PHIols.csv", DataFrame(PHIhat))
CSV.write(savedir * sName * "_SIGMAols.csv", DataFrame(SIGMAhat))
CSV.write(savedir * sName * "_PHIpmean.csv", DataFrame(PHIpmean))
SIGMApmean = SIGMAtrpmean*SIGMAtrpmean'
CSV.write(savedir * sName * "_SIGMApmean.csv", DataFrame(SIGMApmean))
CSV.write(savedir * sName * "_SIGMAtrpmean.csv", DataFrame(SIGMAtrpmean))
