using DataFrames
using QuadGK
using Distributions
using CSV
using LinearAlgebra
using JLD
using Interpolations
import Statistics
#clearconsole()

#-------------------------------------------------------------
# include Functions
#-------------------------------------------------------------
#cd("$(pwd())/Dropbox/Heterogeneity/Software/KS_Simulation/")
readDir = "$(pwd())/CB-fVAR/OVERALL/Functions/"

include(readDir *"vech.jl");
include(readDir *"logSpline_Procedures.jl");
include(readDir *"VAR_Procedures.jl");
include(readDir *"Loaddata.jl");
include(readDir *"EmpPercentiles_Procedures.jl")
dataDir = "$(pwd())/CB-fVAR/OVERALL/Data/"
#-------------------------------------------------------------
# choose specification files
#-------------------------------------------------------------
nfVARSpec   = "10tc"
K           = 6
SampleStart = 1998
SampleEnd   = 2019

nKSpec    = "K$(K)_"

specDir   = "$(pwd())/CB-fVAR/OVERALL/SpecFiles/"
include(specDir * "/fVARspec" * nfVARSpec * ".jl")

#-------------------------------------------------------------
# load aggregate data
#-------------------------------------------------------------
juliaversion = 15
agg_data, period_agg = loadaggdata(SampleStart,SampleEnd,juliaversion)
T     = size(agg_data)[1]
n_agg = size(agg_data)[2]

#-------------------------------------------------------------
# Load coefficients from density estimation
#-------------------------------------------------------------
sName   = "fVAR" * nfVARSpec
loaddir  = "$(pwd())/CB-fVAR/OVERALL/results/" * sName *"/";

knots_all    = CSV.read(loaddir * sName * "_knots_all.csv", DataFrame, header = true);
PhatDensCoef = CSV.read(loaddir * nKSpec * sName * "_PhatDensCoef.csv", DataFrame, header = true);
period_Dens  = CSV.read(loaddir * nKSpec * sName * "_DensityPeriod.csv", DataFrame, header = true);

knots_all           = Matrix(knots_all)'
PhatDensCoef        = Matrix(PhatDensCoef)
period_Dens         = Matrix(period_Dens)
period_Dens_ind     = dropdims((SampleStart .<= period_Dens .<= SampleEnd),dims=2)
period_Dens         = period_Dens[period_Dens_ind]
PhatDensCoef        = PhatDensCoef[period_Dens_ind,:]

# to undo the transformation using inverse hyperbolic sine transform, use
theta = 1.0

#-------------------------------------------------------------
# Select Knots
#-------------------------------------------------------------

ii=getindex.(findall(K_vec.-K.==0),[1 2])[1] # find index ii where K==K_vec
knots = knots_all[quant_sel[ii,:].==1]

#-------------------------------------------------------------
# Initializations
#-------------------------------------------------------------

# compute empirical CDF for a grid of values
ngrid     = 20
grid_temp = range(0.2, stop=1.8, length=ngrid);
ngrid     = 50
grid_temp = range(0.05, stop=2.5, length=ngrid);
grid_temp = [0.0; grid_temp]

# initalize the matrices for the time series of percentiles
vec_percs = [0.2;0.4; 0.5; 0.6; 0.8]
emp_percs_all = zeros(T,length(vec_percs))

#-------------------------------------------------------------
# Loop over time periods to compute percentiles
#-------------------------------------------------------------

for tt = 1:T

    print("CDF and Percentiles for Period = $tt \n")

    unrate    = Matrix(CSV.read(dataDir * "No_superrate.csv", DataFrame, header = true))[tt,2]
    @time emp_percs = DensPercentiles(PhatDensCoef[tt,:]', knots, xgrid, grid_temp, vec_percs)

    emp_percs_all[tt,:] = emp_percs'

end

savedir  = "$(pwd())/CB-fVAR/OVERALL/results/" * sName *"/";
CSV.write(savedir * nKSpec * sName * "_PredPctL_MLE.csv", DataFrame(emp_percs_all,:auto));
