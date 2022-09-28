using DataFrames
using QuadGK
using Distributions
using CSV
using LinearAlgebra
using JLD

#clearconsole()
juliaversion = 15 #use 13 for 1.3 and 15 for 1.5

#-------------------------------------------------------------
# IRF Configuration
#-------------------------------------------------------------

H = 60 # maximum horizon
n_every = 10 # use every n_every'th draw from the posterior
sh_size = 3   # shock size, in multiples of standard deviations
sh_id = 1 # sh_id = 1: TFP shock, sh_id = 2: GDP shock, sh_id = 3: Employment shock

#-------------------------------------------------------------
# include Functions
#-------------------------------------------------------------
#cd("$(pwd())/Dropbox/Heterogeneity/Software/KS_Simulation/")
readDir = "$(pwd())/Empirics_Estimation_v2/Functions/"
include(readDir *"vech.jl");
include(readDir *"logSpline_Procedures.jl");
include(readDir *"VAR_Procedures.jl");
include(readDir *"Loaddata.jl");

#-------------------------------------------------------------
# choose specification files
#-------------------------------------------------------------
nVARSpec    = "2" # with 1= percentiles, and  2= gini, frac
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

n_cross = size(data)[2]-n_agg

#-------------------------------------------------------------
# Load posterior draws
#-------------------------------------------------------------
sName    = "AltVAR" * nVARSpec * "_MCMC" * nVARMCMCSpec * "_" * nSpecStr
loadDir  = "$(pwd())/results/" * sName *"/";

PHIpdraw     = load(loadDir * sName * "_PostDraws.jld", "PHIpdraw")
SIGMAtrpdraw = load(loadDir * sName * "_PostDraws.jld", "SIGMAtrpdraw")
PHIpmean     = load(loadDir * sName * "_PostMeans.jld", "PHIpmean")
SIGMAtrpmean = load(loadDir * sName * "_PostMeans.jld", "SIGMAtrpmean")

#-------------------------------------------------------------
# Generate IRFs at posterior mean of (Phi,Sigmatr)
#-------------------------------------------------------------

println("")
println("Generating IRFs at posterior mean... ")
println("")

# implemented for VAR(1)
YY_IRF         = zeros(H+1,n_agg+n_cross)

# hh=1 is steady state
# hh=2 is period of impact

# impact
YY_IRF[1+1,:]  = sh_size*SIGMAtrpmean[:,sh_id]' # TFP shock
#YY_IRF[1+1,:]  = sh_size*SIGMAtrpmean[:,2]' # GDP shock
#YY_IRF[1+1,:]  = sh_size*SIGMAtrpmean[:,3]' # employment shock

# future periods
for hh = 1+2:H+1
    YY_IRF[hh,:]  = YY_IRF[hh-1,:]'*PHIpmean
end

savedir = "$(pwd())/results/" * sName *"/";
try mkdir(savedir) catch; end
CSV.write(savedir * sName * "_IRF_YY_AggSh"*string(sh_id)*"_pmean.csv", DataFrame(YY_IRF,:auto))

#-------------------------------------------------------------
# Generate IRFs for a subset of posterior draws
#-------------------------------------------------------------

println("")
println("Generating IRFs for posterior draws... ")
println("")

n_subseq                 = floor(Int,size(PHIpdraw)[1]/n_every)
YY_IRF_uncertainty       = zeros(H+1,n_agg+n_cross,n_subseq)

# hh=1 is steady state
# hh=2 is period of impact

for pp = 1:n_subseq

    time_init_loop = time_ns();
    println(" Draw number:  $pp")
    println(" Remaining draws:  $(n_subseq-pp)")
    println(" Posterior draw number: $(pp*n_every)")

    # impact
    YY_IRF_uncertainty[1+1,:,pp] = sh_size*(SIGMAtrpdraw[pp*n_every,:,1])'

    # future periods    z
    for hh = 1+2:H+1
        YY_IRF_uncertainty[hh,:,pp]  = YY_IRF_uncertainty[hh-1,:,pp]'*PHIpdraw[pp*n_every,:,:]
    end

    CSV.write(savedir * sName * "_IRF_YY_AggSh" * string(sh_id) * "_" * string(pp) * ".csv", DataFrame(YY_IRF_uncertainty[:,:,pp],:auto))

    time_loop=signed(time_ns()-time_init_loop)/1000000000
    println(" Elapsed time = $(time_loop)")
    println("")

end

#-------------------------------------------------------------
# No interaction between aggregate and cross-section
#-------------------------------------------------------------

PHIpmean_nointeract = copy(PHIpmean)
PHIpmean_nointeract[n_agg+1:end, 1:n_agg] = zeros(n_cross, n_agg)
PHIpmean_nointeract[1:n_agg, n_agg+1:end] = zeros(n_agg, n_cross)
PHIpmean_nointeract # no interaction between aggregate and cross-section

#-------------------------------------------------------------
# Generate IRFs at posterior mean of (Phi,Sigmatr)
#-------------------------------------------------------------

println("")
println("Generating IRFs at posterior mean... ")
println("")

# implemented for VAR(1)
YY_IRF         = zeros(H+1,n_agg+n_cross)

# hh=1 is steady state
# hh=2 is period of impact

# impact
YY_IRF[1+1,:]  = sh_size*SIGMAtrpmean[:,sh_id]' # TFP shock
#YY_IRF[1+1,:]  = sh_size*SIGMAtrpmean[:,2]' # GDP shock
#YY_IRF[1+1,:]  = sh_size*SIGMAtrpmean[:,3]' # employment shock

# future periods
for hh = 1+2:H+1
    YY_IRF[hh,:]  = YY_IRF[hh-1,:]'*PHIpmean_nointeract
end

savedir = "$(pwd())/results/" * sName *"/";
try mkdir(savedir) catch; end
CSV.write(savedir * sName * "_IRF_nointeract_YY_AggSh"*string(sh_id)*"_pmean.csv", DataFrame(YY_IRF,:auto))
