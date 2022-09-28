using DataFrames
using QuadGK
using Distributions
using CSV
using LinearAlgebra
using JLD
using Interpolations

#clearconsole()
juliaversion  = 15
#-------------------------------------------------------------
# IRF Configuration
#-------------------------------------------------------------

H = 15 # maximum horizon
sh_size = 3   # shock size, in multiples of standard deviations
n_drawsread = 1000 #100 # transform draws 1 to n_drawsread

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
include(readDir *"GiniFracBelowCutoff_Procedures.jl")

#-------------------------------------------------------------
# choose specification files
#-------------------------------------------------------------
nfVARSpec = "10tc"
nModSpec  = "1"
nMCMCSpec = "1"
modName   = "SS"  # VAR or SS
cutoff    = 0.5;
nShockspec= "DistrSh"
#nShockspec= "DistrSh" "AggSh1"
# Aggregate shocks 1: TFP shock, 2: GDP shock,  3: Employment shock
# Distr shock DistrSh is Gini on impact

specDir   = "$(pwd())/CB-fVAR/OVERALL/SpecFiles"
include(specDir * "/fVARspec" * nfVARSpec * ".jl")
include(specDir * "/" * modName * "spec" * nModSpec * ".jl")
include(specDir * "/" * modName * "MCMCspec" * nMCMCSpec * ".jl")

# Percentiles of interest
ngrid     = 20
grid_temp = range(0.2, stop=1.8, length=ngrid);
grid_temp = [0.0; grid_temp]
vec_percs = [0.1; 0.2; 0.5; 0.8; 0.9]

#-------------------------------------------------------------
# load aggregate data
#-------------------------------------------------------------

agg_data, period_agg, mean_unrate = loadaggdata(SampleStart,SampleEnd,juliaversion)
n_agg = size(agg_data)[2]

#-------------------------------------------------------------
# Load coefficients from density estimation
#-------------------------------------------------------------
sNameLoadDir = "fVAR" * nfVARSpec
loaddir  = "$(pwd())/CB-fVAR/OVERALL/results/" * sNameLoadDir *"/";

knots_all = CSV.read(loaddir * sNameLoadDir * "_knots_all.csv", DataFrame, header = true);
knots_all = Array(knots_all)'

ii = getindex.(findall(K_vec.-K.==0),[1 2])[1] # find index ii where K==K_vec
knot = knots_all[quant_sel[ii,:].==1]

PhatDensCoef_factor, MDD_term1, VinvLam_all, period_Dens, PhatDensCoef_lambda, PhatDensCoef_mean, PhatDensCoef_mean_allt = loaddensdata(SampleStart,SampleEnd,K,nfVARSpec,juliaversion)
n_cross = size(PhatDensCoef_factor)[2]

#-------------------------------------------------------------
# Load posterior means
#-------------------------------------------------------------

sName    = "fVAR" * nfVARSpec * "_" * modName * nModSpec * "_" * "MCMC" * nMCMCSpec;
loadDir  = "$(pwd())/CB-fVAR/OVERALL/results/" * sName *"/";

PHIpdraw     = load(loadDir * sName * "_PostDraws.jld", "PHIpdraw")
#n_subseq                 = floor(Int,size(PHIpdraw)[1]/n_every)

YY_IRF = CSV.read(loadDir * sName * "_IRF_YY_" *nShockspec * "_pmean.csv", DataFrame, header = true);
PhatDens_IRF = CSV.read(loadDir * sName * "_IRF_PhatDens_" *nShockspec*"_pmean.csv", DataFrame, header = true);
YY_IRF = Array(YY_IRF)
PhatDens_IRF = Array(PhatDens_IRF)

#-------------------------------------------------------------
# Generate transformed IRFs at posterior mean of (Phi,Sigmatr)
#-------------------------------------------------------------

println("")
println("Generating transformed IRFs at posterior mean... ")
time_init_loop = time_ns();

emp_percs_all    = zeros(H+1,length(vec_percs))
densvalues_ss    = PhatDens_IRF[1,:]; # steady state density
PhatMassDiff_IRF = zeros(H+1,1);
Gini_IRF         = zeros(H+1,1);

for hh = 1:(H+1)
    # read estimated states in period t, 1:n_agg are agg data, n_agg+1:end are dens coeff
    statepmean_hh = YY_IRF[hh,:]
    # statepmean contains dev of unemployment from its mean -> add mean to obtain levels
    unrate_hh = statepmean_hh[3] + mean_unrate

    PhatDensCoef_hh    = coefRecover(statepmean_hh[n_agg+1:end]', PhatDensCoef_lambda, PhatDensCoef_mean)
    emp_percs_all[hh,:] = DensPercentiles(PhatDensCoef_hh, knots, unrate_hh, xgrid, grid_temp, vec_percs)
#    emp_percs_all[hh,:] = DensPercentiles(statepmean_hh[n_agg+1:end]', knots, unrate_hh, xgrid, grid_temp, vec_percs)
    #PhatMassDiff_IRF[hh,1] = ProbMassDiff_Cutoff(PhatDens_IRF[hh,:], densvalues_ss, xgrid);
    #Gini_IRF[hh] = GiniCoef(PhatDens_IRF[hh,:], xgrid)

end

savedir = "$(pwd())/CB-fVAR/OVERALL/results/" * sName *"/";
CSV.write(savedir * sName * "_IRF_Pctl_" *nShockspec*"_pmean.csv", DataFrame(emp_percs_all,:auto))
#CSV.write(savedir * sName * "_IRF_BelowCutoff_" *nShockspec*"_pmean.csv", DataFrame(PhatMassDiff_IRF,:auto))
#CSV.write(savedir * sName * "_IRF_Gini_" *nShockspec*"_pmean.csv", DataFrame(Gini_IRF,:auto))

time_loop=signed(time_ns()-time_init_loop)/1000000000
println("Elapsed time = $(time_loop)")
println("")


#-------------------------------------------------------------
# Generate IRFs for a subset of posterior draws
#-------------------------------------------------------------

println("")
println("Generating transformed IRFs for posterior draws... ")
println("")

YY_IRF_uncertainty        = zeros(H+1,n_agg+n_cross,n_drawsread)
PhatDens_IRF_uncertainty  = zeros(H+1,length(xgrid),n_drawsread)

if juliaversion == 13
    for pp = 1:n_drawsread
    YY_IRF_pp = CSV.read(loadDir * sName * "_IRF_YY_" *nShockspec * "_" * string(pp) * ".csv", header = true);
    PhatDens_IRF_pp = CSV.read(loadDir * sName * "_IRF_PhatDens_" *nShockspec*  "_" * string(pp) * ".csv", header = true);

    YY_IRF_uncertainty[:,:,pp] = Matrix(YY_IRF_pp)
    PhatDens_IRF_uncertainty[:,:,pp] = Matrix(PhatDens_IRF_pp)
    end
else
    for pp = 1:n_drawsread
    YY_IRF_pp = CSV.read(loadDir * sName * "_IRF_YY_" *nShockspec * "_"*string(pp) *".csv", DataFrame, header = true);
    PhatDens_IRF_pp = CSV.read(loadDir * sName * "_IRF_PhatDens_" *nShockspec* "_" * string(pp) * ".csv", DataFrame, header = true);

    YY_IRF_uncertainty[:,:,pp] = Matrix(YY_IRF_pp)
    PhatDens_IRF_uncertainty[:,:,pp] = Matrix(PhatDens_IRF_pp)
    end
end

emp_percs_uncertainty        = zeros(H+1,length(vec_percs),n_drawsread)
PhatMassDiff_IRF_uncertainty = zeros(H+1,n_drawsread);
Gini_IRF_uncertainty         = zeros(H+1,n_drawsread);

errort = 0
for pp = 1:n_drawsread

    time_init_loop = time_ns();
    println("Draw number:  $pp")
    println("Remaining draws:  $(n_drawsread-pp)")

    for hh = 1:(H+1)
        # read estimated states in period t, 1:n_agg are agg data, n_agg+1:end are dens coeff
        statepmean_hh = YY_IRF_uncertainty[hh,:,pp]
        # statepmean contains dev of unemployment from its mean -> add mean to obtain levels
        #unrate_hh = statepmean_hh[3] + mean_unrate
        
        
        PhatDensCoef_hh    = coefRecover(statepmean_hh[n_agg+1:end]', PhatDensCoef_lambda, PhatDensCoef_mean)
        emp_percs_uncertainty[hh,:,pp] = DensPercentiles(PhatDensCoef_hh, knots, xgrid, grid_temp, vec_percs)
        #emp_percs_uncertainty[hh,:,pp] = DensPercentiles(statepmean_hh[n_agg+1:end]', knots, unrate_hh, xgrid, grid_temp, vec_percs)
        #PhatMassDiff_IRF_uncertainty[hh,pp] = ProbMassDiff_Cutoff(PhatDens_IRF_uncertainty[hh,:,pp], densvalues_ss, cutoff, xgrid);
        #Gini_IRF_uncertainty[hh,pp] = GiniCoef(PhatDens_IRF_uncertainty[hh,:,pp], xgrid)
    end

    time_loop=signed(time_ns()-time_init_loop)/1000000000
    println("Elapsed time = $(time_loop)")
    println("")
    CSV.write(savedir * sName * "_IRF_Pctl_" *nShockspec * "_" * string(pp)* ".csv", DataFrame(emp_percs_uncertainty[:,:,pp],:auto))

end



for pp = 1:n_drawsread

    time_init_loop = time_ns();
    println("Draw number:  $pp")
    println("Remaining draws:  $(n_drawsread-pp)")

    for hh = 1:(H+1)
        # read estimated states in period t, 1:n_agg are agg data, n_agg+1:end are dens coeff
        
        # statepmean contains dev of unemployment from its mean -> add mean to obtain levels
#         unrate_hh = statepmean_hh[3] + mean_unrate

        
        try
            # read estimated states in period t, 1:n_agg are agg data, n_agg+1:end are dens coeff
            statepmean_hh = YY_IRF_uncertainty[hh,:,pp]
            # statepmean contains dev of unemployment from its mean -> add mean to obtain levels
            unrate_hh = statepmean_hh[3] + mean_unrate

            PhatDensCoef_hh    = coefRecover(statepmean_hh[n_agg+1:end]', PhatDensCoef_lambda, PhatDensCoef_mean)
            emp_percs_uncertainty[hh,:,pp] = DensPercentiles(PhatDensCoef_hh, knots, unrate_hh, xgrid, grid_temp, vec_percs)

        catch
            statepmean_hh = YY_IRF_uncertainty[hh,:,(pp-1)]
            PhatDensCoef_hh    = coefRecover(statepmean_hh[n_agg+1:end]', PhatDensCoef_lambda, PhatDensCoef_mean)
            
            emp_percs_uncertainty[hh,:,pp] = DensPercentiles(PhatDensCoef_hh, knots, xgrid, grid_temp, vec_percs)
            errort = errort + 1
        end
        #emp_percs_uncertainty[hh,:,pp] = DensPercentiles(statepmean_hh[n_agg+1:end]', knots, unrate_hh, xgrid, grid_temp, vec_percs)
        #PhatMassDiff_IRF_uncertainty[hh,pp] = ProbMassDiff_Cutoff(PhatDens_IRF_uncertainty[hh,:,pp], densvalues_ss, xgrid);
        #Gini_IRF_uncertainty[hh,pp] = GiniCoef(PhatDens_IRF_uncertainty[hh,:,pp], xgrid)
    end

    time_loop=signed(time_ns()-time_init_loop)/1000000000
    println("Elapsed time = $(time_loop)")
    println("")
    CSV.write(savedir * sName * "_IRF_Pctl_" *nShockspec * "_" * string(pp)* ".csv", DataFrame(emp_percs_uncertainty[:,:,pp],:auto))

end

#CSV.write(savedir * sName * "_IRF_BelowCutoff_" *nShockspec*"_uncertainty.csv", DataFrame(PhatMassDiff_IRF_uncertainty,:auto))
#CSV.write(savedir * sName * "_IRF_Gini_" *nShockspec*"_uncertainty.csv", DataFrame(Gini_IRF_uncertainty,:auto))
