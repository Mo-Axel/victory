using DataFrames
using QuadGK
using Distributions
using CSV
using LinearAlgebra
using JLD
using Interpolations

clearconsole()

#-------------------------------------------------------------
# IRF Configuration
#-------------------------------------------------------------

H = 40 # maximum horizon
n_every = 10 # use every n_every'th draw from the posterior
sh_size = 3   # shock size, in multiples of standard deviations
sh_id = 3 # sh_id = 1: TFP shock, sh_id = 2: GDP shock, sh_id = 3: Employment shock


#-------------------------------------------------------------
# include Functions
#-------------------------------------------------------------
#cd("$(pwd())/Dropbox/Heterogeneity/Software/KS_Simulation/")
readDir = "$(pwd())/Functions/"
include(readDir *"vech.jl");
include(readDir *"logSpline_Procedures.jl");
include(readDir *"VAR_Procedures.jl");
include(readDir *"Loaddata.jl");
include(readDir *"EmpPercentiles_Procedures.jl")
include(readDir *"GiniFracBelowCutoff_Procedures.jl")
include(readDir *"IRF_Procedures.jl")

#-------------------------------------------------------------
# choose specification files
#-------------------------------------------------------------
nfVARSpec = "10tc"
nModSpec  = "4"
nMCMCSpec = "1"
modName   = "SS"  # VAR or SS

specDir   = "$(pwd())/SpecFiles/"
include(specDir * "/fVARspec" * nfVARSpec * ".jl")
include(specDir * "/" * modName * "spec" * nModSpec * ".jl")
include(specDir * "/" * modName * "MCMCspec" * nMCMCSpec * ".jl")

#-------------------------------------------------------------
# load aggregate data
#-------------------------------------------------------------
juliaversion = 15;
agg_data, period_agg, mean_unrate = loadaggdata(SampleStart,SampleEnd,juliaversion)
n_agg = size(agg_data)[2]

#-------------------------------------------------------------
# Load coefficients from density estimation
#-------------------------------------------------------------
sNameLoadDir = "fVAR" * nfVARSpec
loaddir  = "$(pwd())/results/" * sNameLoadDir *"/";

knots_all = CSV.read(loaddir * sNameLoadDir * "_knots_all.csv", DataFrame, header = true);
knots_all = convert(Array, knots_all)'

ii=getindex.(findall(K_vec.-K.==0),[1 2])[1] # find index ii where K==K_vec
knots = knots_all[quant_sel[ii,:].==1]

PhatDensCoef_factor, MDD_term1, VinvLam_all, period_Dens, PhatDensCoef_lambda, PhatDensCoef_mean, PhatDensCoef_mean_allt = loaddensdata(SampleStart,SampleEnd,K,nfVARSpec,juliaversion)
n_cross = size(PhatDensCoef_factor)[2]

#-------------------------------------------------------------
# Load posterior mean and draws
#-------------------------------------------------------------

sName    = "fVAR" * nfVARSpec * "_" * modName * nModSpec * "_" * "MCMC" * nMCMCSpec;
loadDir  = "$(pwd())/results/" * sName *"/";

PHIpdraw     = load(loadDir * sName * "_PostDraws.jld", "PHIpdraw")
SIGMAtrpdraw = load(loadDir * sName * "_PostDraws.jld", "SIGMAtrpdraw")
PHIpmean     = load(loadDir * sName * "_PostMeans.jld", "PHIpmean")
SIGMAtrpmean = load(loadDir * sName * "_PostMeans.jld", "SIGMAtrpmean")

#-------------------------------------------------------------
# Generate IRFs at posterior mean of (Phi,Sigmatr)
# hh=1 is steady state
# hh=2 is period of impact
#-------------------------------------------------------------

qstar = zeros(n_agg+n_cross)
qstar[sh_id] = 1

println("")
println("Generating IRFs at posterior mean... ")
println("")

YY_IRF,PhatDens_IRF,~ = IRF_qSh(PHIpmean, SIGMAtrpmean, qstar, sh_size, H, xgrid)

savedir = "$(pwd())/results/" * sName *"/";
try mkdir(savedir) catch; end
CSV.write(savedir * sName * "_IRF_PhatDens_AggSh"*string(sh_id)*"_pmean.csv", DataFrame(PhatDens_IRF))
CSV.write(savedir * sName * "_IRF_YY_AggSh"*string(sh_id)*"_pmean.csv", DataFrame(YY_IRF))

#-------------------------------------------------------------
# Generate IRFs for a subset of posterior draws
#-------------------------------------------------------------

println("")
println("Generating IRFs for subset of posterior draws... ")
println("")

n_subseq                 = floor(Int,size(PHIpdraw)[1]/n_every)

for pp = 1:n_subseq

    time_init_loop = time_ns();
    println("Draw number:  $pp")
    println("Remaining draws:  $(n_subseq-pp)")
    println("Posterior draw number: $(pp*n_every)")

    YY_IRF,PhatDens_IRF,~ = IRF_qSh(PHIpdraw[pp*n_every,:,:], SIGMAtrpdraw[pp*n_every,:,:], qstar, sh_size, H, xgrid)

    CSV.write(savedir * sName * "_IRF_PhatDens_AggSh" * string(sh_id) * "_" * string(pp) * ".csv", DataFrame(PhatDens_IRF))
    CSV.write(savedir * sName * "_IRF_YY_AggSh" * string(sh_id) * "_" * string(pp) * ".csv", DataFrame(YY_IRF))

    time_loop=signed(time_ns()-time_init_loop)/1000000000
    println("Elapsed time = $(time_loop)")
    println("")

end
