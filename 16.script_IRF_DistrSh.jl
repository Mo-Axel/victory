using DataFrames
using QuadGK
using Distributions
using CSV
using LinearAlgebra
using JLD
using Interpolations

#clearconsole()

#-------------------------------------------------------------
# IRF Configuration
#-------------------------------------------------------------

H  = 15 # maximum horizon
Hq = 1 # max horizon needed to determine q
n_every = 10 # use every n_every'th draw from the posterior
sh_size = 3   # shock size, in multiples of standard deviations


#-------------------------------------------------------------
# include Functions
#-------------------------------------------------------------
#cd("$(pwd())/Dropbox/Heterogeneity/Software/KS_Simulation/")
readDir = "$(pwd())/CB-fVAR/OVERALL/Functions/"
include(readDir *"vech.jl");
include(readDir *"logSpline_Procedures.jl");
include(readDir *"VAR_Procedures.jl");
include(readDir *"Loaddata.jl");
include(readDir *"GiniFracBelowCutoff_Procedures.jl")
include(readDir *"IRF_Procedures.jl")



#-------------------------------------------------------------
# choose specification files
#-------------------------------------------------------------
nfVARSpec = "10tc"
nModSpec  = "4"
nMCMCSpec = "1"
modName   = "SS"  # VAR or SS

specDir   = "$(pwd())/CB-fVAR/OVERALL/SpecFiles/"
include(specDir * "/fVARspec" * nfVARSpec * ".jl")
include(specDir * "/" * modName * "spec" * nModSpec * ".jl")
include(specDir * "/" * modName * "MCMCspec" * nMCMCSpec * ".jl")

#-------------------------------------------------------------
# load aggregate data
#-------------------------------------------------------------
juliaversion = 15
agg_data, period_agg = loadaggdata(SampleStart,SampleEnd,juliaversion)
n_agg = size(agg_data)[2]

#-------------------------------------------------------------
# Load coefficients from density estimation
#-------------------------------------------------------------
sNameLoadDir = "fVAR" * nfVARSpec
loaddir  = "$(pwd())/CB-fVAR/OVERALL/results/" * sNameLoadDir *"/";
knots_all = CSV.read(loaddir * sNameLoadDir * "_knots_all.csv", DataFrame, header = true);

knots_all = Array(knots_all)'
ii=getindex.(findall(K_vec.-K.==0),[1 2])[1] # find index ii where K==K_vec
knots = knots_all[quant_sel[ii,:].==1]

PhatDensCoef_factor, MDD_term1, VinvLam_all, period_Dens, PhatDensCoef_lambda, PhatDensCoef_mean, PhatDensCoef_mean_allt = loaddensdata(SampleStart,SampleEnd,K,nfVARSpec,juliaversion)
n_cross = size(PhatDensCoef_factor)[2]

#-------------------------------------------------------------
# Load posterior draws
#-------------------------------------------------------------

sName    = "fVAR" * nfVARSpec * "_" * modName * nModSpec * "_" * "MCMC" * nMCMCSpec;
loadDir  = "$(pwd())/CB-fVAR/OVERALL/results/" * sName *"/";

PHIpdraw     = load(loadDir * sName * "_PostDraws.jld", "PHIpdraw")
SIGMAtrpdraw = load(loadDir * sName * "_PostDraws.jld", "SIGMAtrpdraw")
PHIpmean     = load(loadDir * sName * "_PostMeans.jld", "PHIpmean")
SIGMAtrpmean = load(loadDir * sName * "_PostMeans.jld", "SIGMAtrpmean")

#-------------------------------------------------------------
# Create a grid for the rotation q
#-------------------------------------------------------------

lenq = 5000 # out of 500 rotations, find the one that maximizes the Gini coefficients
#Random.seed!(100000*seedindx+10*t+3+seedoffset)
qgrid = randn(lenq,n_cross)

for qq = 1:lenq
   qgrid[qq,:] = qgrid[qq,:]./norm(qgrid[qq,:])
end
qgrid = [zeros(lenq,n_agg) qgrid]; # save candidate shocks

#-------------------------------------------------------------
# Generate IRFs at posterior mean of (Phi,Sigmatr)
# hh=1 is steady state
# hh=2 is period of impact
# implemented for VAR(1)
#-------------------------------------------------------------
println("")
println("Generating IRFs at posterior mean... ")
println("")

# find the optimal q
Gini_vec = zeros(lenq,1)
for qq = 1:lenq
    YY_IRF,~,Gini_IRF = IRF_qSh(PHIpmean, SIGMAtrpmean, qgrid[qq,:], sh_size, Hq, xgrid)
#     return(YY_IRF)
    Gini_vec[qq] = Gini_IRF[2]
end
maxGini_qstar = qgrid[getindex.(findall(Gini_vec.-maximum(Gini_vec).==0),[1 2])[1],:]

# Recompute the IRFs
YY_IRF,PhatDens_IRF,~ = IRF_qSh(PHIpmean, SIGMAtrpmean, maxGini_qstar, sh_size, H, xgrid)
savedir = "$(pwd())/CB-fVAR/OVERALL/results/" * sName *"/";
try mkdir(savedir) catch; end
CSV.write(savedir * sName * "_IRF_PhatDens_DistrSh_pmean.csv", DataFrame(PhatDens_IRF,:auto))
CSV.write(savedir * sName * "_IRF_YY_DistrSh_pmean.csv", DataFrame(YY_IRF,:auto))

#-------------------------------------------------------------
# Generate IRFs for a subset of posterior draws
#-------------------------------------------------------------
println("")
println("Generating IRFs for subset of posterior draws... ")
println("")

n_subseq                 = floor(Int,size(PHIpdraw)[1]/n_every)
errort = 0
for pp = 1:n_subseq

    time_init_loop = time_ns();
    println("Draw number:  $pp")
    println("Remaining draws:  $(n_subseq-pp)")
    println("Posterior draw number: $(pp*n_every)")

    # find the optimal q
    Gini_vec = zeros(lenq,1)
    for qq = 1:lenq
        try
            ~,~,Gini_IRF = IRF_qSh(PHIpdraw[pp*n_every,:,:], SIGMAtrpdraw[pp*n_every,:,:], qgrid[qq,:], sh_size, Hq, xgrid)
            Gini_vec[qq] = Gini_IRF[1+1]
        catch
            ~,~,Gini_IRF = IRF_qSh(PHIpdraw[pp*n_every,:,:], SIGMAtrpdraw[pp*n_every,:,:], qgrid[(qq-1),:], sh_size, Hq, xgrid)
            Gini_vec[qq] = Gini_IRF[1+1]
            errort = errort + 1
        end
    end
    maxGini_qstar = qgrid[getindex.(findall(Gini_vec.-maximum(Gini_vec).==0),[1 2])[1],:]
    # Recompute the IRFs
    try
        YY_IRF,PhatDens_IRF,~ = IRF_qSh(PHIpdraw[pp*n_every,:,:], SIGMAtrpdraw[pp*n_every,:,:], maxGini_qstar, sh_size, H, xgrid)
    catch
        YY_IRF,PhatDens_IRF,~ = IRF_qSh(PHIpdraw[(pp-1)*n_every,:,:], SIGMAtrpdraw[(pp-1)*n_every,:,:], maxGini_qstar, sh_size, H, xgrid)
        errort = errort + 1
    end

    CSV.write(savedir * sName * "_IRF_PhatDens_DistrSh" * "_" * string(pp) * ".csv", DataFrame(PhatDens_IRF,:auto))
    CSV.write(savedir * sName * "_IRF_YY_DistrSh" * "_" * string(pp) * ".csv", DataFrame(YY_IRF,:auto))

    time_loop=signed(time_ns()-time_init_loop)/1000000000
    println("Elapsed time = $(time_loop)")
    println("")

end
