using DataFrames
using QuadGK
using Distributions
using JLD
using Optim
using ForwardDiff
using CSV
using LinearAlgebra

#clearconsole()

#-------------------------------------------------------------
# include Functions
#-------------------------------------------------------------
#cd("$(pwd())/Dropbox/Heterogeneity/Software/Empirical_Analysis/")
readDir = "$(pwd())/CB-fVAR/OVERALL/Functions/"
include(readDir *"logSpline_Procedures.jl");

#-------------------------------------------------------------
# load data
#-------------------------------------------------------------
dataDir = "$(pwd())/CB-fVAR/OVERALL/Data/"

unrate_data     = CSV.read(dataDir * "No_superrate.csv", DataFrame, header = true);
earnings_data   = CSV.read(dataDir * "cross_all.csv", DataFrame, header = true);

unrate     = Matrix(unrate_data)[:,2]
#period_UNR = Matrix(unrate_data)[:,1]
earnings_detrended = Matrix(earnings_data)[:,3]
earnings_t = Matrix(earnings_data)[:,2]
Tend       = length(unrate)

#-------------------------------------------------------------
# choose specification file
#-------------------------------------------------------------
nfVARSpec =  "10tc"
specDir   = "$(pwd())/CB-fVAR/SpecFiles/"
include(specDir * "/fVARspec" * nfVARSpec * ".jl")

# subsequently use the same knots regardless of sample size N,T
knots_all = quantile(earnings_detrended, quant_vec)

#-------------------------------------------------------------
# log spline density estimation over K and t
#-------------------------------------------------------------
MDD_GoF_sum  = zeros(K_vec_n, 3);

for ii = 1:K_vec_n
    K             = K_vec[ii]
    knots         = knots_all[quant_sel[ii,:].==1]
    PhatDensValue = zeros(Tend, length(xgrid));
    PhatDensCoef  = zeros(Tend, K);
    PhatDensNorm  = zeros(Tend, 1);
    PhatlogLike   = zeros(Tend, 1);
    Vinv_all      = zeros(K*Tend, K);
    N_all         = zeros(Tend, 1);
    Period_all    = zeros(Tend, 1);
    N_details     = zeros(Tend, 4+K);

    timeidx = 1998

    for tt = 1:Tend

        time_init_loop = time_ns()
        Period_all[tt] = timeidx

        # time t data
        selecteddraws_t = earnings_detrended[earnings_t.==timeidx]
        timeidx = timeidx + 1
        N_all[tt]       = length(selecteddraws_t)

        # count observations with knot restriction
        # recall that there are K-1 knots
        N_knots      = zeros(1,K);
        N_knots[1,1] = sum(selecteddraws_t.<=knots[1])
        for kk = 2:(K-1)
            N_knots[1,kk] = sum((selecteddraws_t.>knots[kk-1]) .& (selecteddraws_t.<=knots[kk]))
        end
        N_knots[1,K] = sum(selecteddraws_t.>knots[end])
        println("Number of obs in knot brackets: $(N_knots)")
        println("Max knot:        $(knots[K-1])")

        # compute MLE
        if tt == 1
            alpha_initial = zeros(K)
        else
            alpha_initial = PhatDensCoef[tt-1,:]
        end

        if TopCodeFlag == 0
            f1(x)  = logspline_obj(selecteddraws_t, x, knots, minimum(xgrid), maximum(xgrid)) # without top-coding
            td     = TwiceDifferentiable(f1, alpha_initial; autodiff  = :forward )
            results_t = try
                optimize(td, alpha_initial, Newton())
            catch
                optimize(f1, alpha_initial, BFGS())
            end
        else
            f2(x)  = logspline_obj_topcode(selecteddraws_t, x, knots, minimum(xgrid), maximum(xgrid)) # with top-coding
            td     = TwiceDifferentiable(f2, alpha_initial; autodiff  = :forward )
            results_t = try
                optimize(td, alpha_initial, Newton())
            catch
                optimize(f2, alpha_initial, BFGS())
            end
        end

        coef_t = results_t.minimizer
        C_topcode = maximum(selecteddraws_t)
        N_max     = sum(selecteddraws_t.==C_topcode)

        if (TopCodeFlag == 0) | (N_max == 1)
            # likelihood w/o top coding
            PhatlogLike[tt] = - N_all[tt]*results_t.minimum
            pi_hat = 0
            println("No top coding")
        else
            # likelihood w top coding
            println("Top coded value is $(C_topcode), number of obs is $(N_max)")
            pi_hat = N_max/N_all[tt]
            PhatlogLike[tt] = - N_all[tt]*results_t.minimum + (N_all[tt]-N_max)*log(1-pi_hat) + N_max*log(pi_hat)
        end

        # results
        PhatDensCoef[tt,:]  = coef_t
        PhatDensNorm[tt,:]  = lnpdfNormalize_unrate(coef_t',knots, unrate, minimum(xgrid), maximum(xgrid))
        PhatDensValue[tt,:] = pdfEval(xgrid,coef_t,knots,[PhatDensNorm[tt,1]])';
        N_details[tt,:]     = [N_all[tt] N_max pi_hat C_topcode N_knots]

        # compute negative inverse hessian
        # changed type of hessian_sqrtK output to "Symmetric"
        # note that V_t type is also Symmetric

        if (TopCodeFlag == 0) | (N_max == 1)
            # Hessian w/o top coding
            Hess_t = hessian_loglh(PhatDensCoef[tt,:], knots, minimum(xgrid), maximum(xgrid))
        else
            # Hessian w top coding
            Hess_t = (N_all[tt]-N_max)/N_all[tt]*hessian_loglh(PhatDensCoef[tt,:], knots, minimum(xgrid), C_topcode)
        end

        Vinv_t = - Hess_t
        if isposdef(Vinv_t) == false
            Vinv_eig = eigen(Vinv_t)
            Vinv_t = Symmetric(Vinv_eig.vectors*Diagonal(abs.(Vinv_eig.values))*Vinv_eig.vectors')
            println("Flipped neg eigenvalues")
        end
        Vinv_all[K*(tt-1)+1:K*tt,:] = Vinv_t

        println("K = $(K), Period $(tt)")
        println("Vinv_t is positive definite: $(isposdef(Vinv_t))")
        time_loop=signed(time_ns()-time_init_loop)/1000000000
        println("Elapsed Time $(time_loop) seconds")
        println("")

    end

    # seasonality adjustment
    nq = 1# quarterly data
    start_season = 1 # starting from 1989Q1 (first quarter)
    (PhatDensCoef_adj, PhatDensCoef_mean, PhatDensCoef_mean_allt) = seasonality_adj(PhatDensCoef, nq, start_season, Tend)

    # compress coefficients, ~ = PhatDensCoef_mean
    (PhatDensCoef_factor, PhatDensCoef_lambda, ~ ) = coefCompress(PhatDensCoef_adj)
    Ktilde = size(PhatDensCoef_factor)[2]
    println("----------------")
    println("Compression Step")
    println("K = $(K), K-tilde = $(Ktilde)")
    println("----------------")
    println("")

    # Goodness of Fit (GoF) is log likelihood
    MDD_GoF   = zeros(Tend, 1);
    for tt = 1:Tend
        MDD_GoF[tt] = PhatlogLike[tt]
    end

    MDD_GoF_sum[ii,:] = [K Ktilde sum(MDD_GoF)]

    # save results
    sNameDir  = "fVAR" * nfVARSpec
    sNameFile = "K" * string(K) * "_fVAR" * nfVARSpec
    savedir = "$(pwd())/CB-fVAR/results/" * sNameDir *"/";
    try mkdir(savedir) catch; end
    CSV.write(savedir * sNameFile * "_DensityPeriod.csv", DataFrame(Period_all,:auto))
    CSV.write(savedir * sNameFile * "_PhatDensValue.csv", DataFrame(PhatDensValue,:auto))
    CSV.write(savedir * sNameFile * "_PhatDensCoef.csv", DataFrame(PhatDensCoef,:auto))
    CSV.write(savedir * sNameFile * "_PhatDensCoef_factor.csv", DataFrame(PhatDensCoef_factor,:auto))
    CSV.write(savedir * sNameFile * "_PhatDensCoef_lambda.csv", DataFrame(PhatDensCoef_lambda,:auto))
    CSV.write(savedir * sNameFile * "_PhatDensCoef_mean.csv", DataFrame(PhatDensCoef_mean,:auto))
    CSV.write(savedir * sNameFile * "_PhatDensCoef_mean_allt.csv", DataFrame(PhatDensCoef_mean_allt,:auto))
    CSV.write(savedir * sNameFile * "_Vinv_all.csv", DataFrame(Vinv_all,:auto))
    CSV.write(savedir * sNameFile * "_N_all.csv", DataFrame(N_all,:auto))
    CSV.write(savedir * sNameFile * "_N_details.csv", DataFrame(N_details,:auto))
    CSV.write(savedir * sNameFile * "_MDD_GoF.csv", DataFrame(MDD_GoF,:auto))

end

sNameDir  = "fVAR" * nfVARSpec;
sNameFile = "fVAR" * nfVARSpec;
savedir = "$(pwd())/CB-fVAR/results/" * sNameDir *"/";
try mkdir(savedir) catch; end
CSV.write(savedir * sNameFile * "_MDD_GoF_sum.csv", DataFrame(MDD_GoF_sum,:auto))
CSV.write(savedir * sNameFile * "_knots_all.csv", DataFrame(knots_all',:auto))
