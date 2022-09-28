using DataFrames
using QuadGK
using Distributions
using CSV
using LinearAlgebra
using JLD
using Interpolations

clearconsole()

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

#-------------------------------------------------------------
# choose specification files
#-------------------------------------------------------------
nfVARSpec = "10tc"
nModSpec  = "1"
nMCMCSpec = "1"
modName   = "SS"  # VAR or SS

specDir   = "$(pwd())/SpecFiles/"
include(specDir * "/fVARspec" * nfVARSpec * ".jl")
include(specDir * "/" * modName * "spec" * nModSpec * ".jl")
include(specDir * "/" * modName * "MCMCspec" * nMCMCSpec * ".jl")

#-------------------------------------------------------------
# load aggregate data to obtain unemployment rate
#-------------------------------------------------------------
juliaversion = 15
agg_data, period_agg, mean_unrate = loadaggdata(SampleStart,SampleEnd,juliaversion)
n_agg = size(agg_data)[2]

#-------------------------------------------------------------
# Load knots from density estimation
#-------------------------------------------------------------
sNameLoadDir = "fVAR" * nfVARSpec
loaddir  = "$(pwd())/results/" * sNameLoadDir *"/";

knots_all = CSV.read(loaddir * sNameLoadDir * "_knots_all.csv", DataFrame, header = true);
knots_all = convert(Array, knots_all)'

#-------------------------------------------------------------
# Generate Lambda matrix and vector of means to convert a factors
#-------------------------------------------------------------
PhatDensCoef_factor, MDD_term1, VinvLam_all, period_Dens, PhatDensCoef_lambda, PhatDensCoef_mean, PhatDensCoef_mean_allt = loaddensdata(SampleStart,SampleEnd,K,nfVARSpec,juliaversion)

# later, undo the transformation using inverse hyperbolic sine transform with theta
theta = 1.0

#-------------------------------------------------------------
# Load state's posterior mean
#-------------------------------------------------------------

sName    = "fVAR" * nfVARSpec * "_" * modName * nModSpec * "_" * "MCMC" * nMCMCSpec;
loadDir  = "$(pwd())/results/" * sName *"/";

statepmean = load(loadDir * sName * "_StateMeans.jld", "statepmean")
T = size(statepmean)[1]

ii=getindex.(findall(K_vec.-K.==0),[1 2])[1] # find index ii where K==K_vec
knots = knots_all[quant_sel[ii,:].==1]

# compute empirical CDF for a grid of values
# ngrid     = 20
# grid_temp = range(0.2, stop=1.8, length=ngrid);
ngrid     = 50
grid_temp = range(0.05, stop=2.5, length=ngrid);
grid_temp = [0.0; grid_temp]

# initalize the matrices for the time series of percentiles
vec_percs = [0.1; 0.2; 0.5; 0.8; 0.9]
emp_percs_all = zeros(T,length(vec_percs))

PhatDensCoef_smoothed = zeros(T,K)
Ktilde = size(statepmean)[2]-n_agg
PhatDensCoef_factor_smoothed = zeros(T,Ktilde)

# use numerical integration to get empirical CDF increments evaluated at grid_temp
# iterate over time period

for tt = 1:T

    print("CDF and Percentiles for Period = $tt \n")

    PhatDensCoef_factor = statepmean[tt,n_agg+1:end]
    PhatDensCoef        = coefRecover(PhatDensCoef_factor',PhatDensCoef_lambda, PhatDensCoef_mean_allt[tt,:]')
    unrate              = agg_data[tt,3] + mean_unrate
    emp_percs           = DensPercentiles(PhatDensCoef, knots, unrate, xgrid, grid_temp, vec_percs)
    emp_percs_all[tt,:] = emp_percs'
    PhatDensCoef_smoothed[tt,:] = PhatDensCoef
    PhatDensCoef_factor_smoothed[tt,:] = PhatDensCoef_factor

end

savedir  = "$(pwd())/results/" * sName *"/";
CSV.write(savedir * sName * "_PredPctl_Smoothed.csv", DataFrame(emp_percs_all));
CSV.write(savedir * sName * "_PhatDensCoef_Smoothed.csv", DataFrame(PhatDensCoef_smoothed));
CSV.write(savedir * sName * "_PhatDensCoef_factor_Smoothed.csv", DataFrame(PhatDensCoef_factor_smoothed));
