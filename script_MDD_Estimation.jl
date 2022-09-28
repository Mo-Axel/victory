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
readDir = "$(pwd())/Functions/"
include(readDir *"vech.jl");
include(readDir *"VAR_Procedures.jl");
include(readDir *"Loaddata.jl");

#-------------------------------------------------------------
# choose specification files
#-------------------------------------------------------------
nfVARSpec = "10tc"
nMDDSpec  = "1"
nMCMCSpec = "1"

specDir   = "$(pwd())/SpecFiles/"
include(specDir * "/fVARspec" * nfVARSpec * ".jl")
include(specDir * "/MDDspec" * nMDDSpec * ".jl")
include(specDir * "/MDDMCMCspec" * nMCMCSpec * ".jl")

#-------------------------------------------------------------
# load aggregate data
#-------------------------------------------------------------
juliaversion = 15 #use 13 for 1.3 and 15 for 1.5
agg_data, period_agg, ~ = loadaggdata(SampleStart,SampleEnd,juliaversion)
n_agg = size(agg_data)[2]

#-------------------------------------------------------------
# Iterate over number of knots in specification file
#-------------------------------------------------------------
sNameLoadDir = "fVAR" * nfVARSpec
sNameSaveDir = "fVAR" * nfVARSpec * "_MDD" * nMDDSpec*"_MCMC" * nMCMCSpec
loaddir  = "$(pwd())/results/" * sNameLoadDir *"/";
savedir  = "$(pwd())/results/" * sNameSaveDir *"/";

try mkdir(savedir) catch; end

MDD_Laplace_sum = zeros(hyper_n, K_vec_n+size(hyper_grid)[2])
MDD_Bayes_sum   = zeros(hyper_n, K_vec_n+size(hyper_grid)[2])
Pen_Laplace_sum = zeros(hyper_n, K_vec_n+size(hyper_grid)[2])
Pen_Bayes_sum   = zeros(hyper_n, K_vec_n+size(hyper_grid)[2])

lambda_hat_K_BayesMax   = zeros(K_vec_n,6)
lambda_hat_K_LaplaceMax = zeros(K_vec_n,6)

for ii = 1:K_vec_n

    K = K_vec[ii]

    #-------------------------------------------------------------
    # Load coefficients from density estimation
    #-------------------------------------------------------------

    PhatDensCoef_factor, MDD_GoF, VinvLam_all, period_Dens_ind, ~, ~, ~ = loaddensdata(SampleStart,SampleEnd,K,nfVARSpec,juliaversion)
    n_cross = size(PhatDensCoef_factor)[2]

    #-------------------------------------------------------------
    # Combine agg data and density coefficients
    #-------------------------------------------------------------

    data    = [ agg_data PhatDensCoef_factor ]
    Ktilde = size(PhatDensCoef_factor)[2]

    #-------------------------------------------------------------
    # Compute MDD term2 over hyperparameter grid
    #-------------------------------------------------------------
    silent  = 1
    sNameFile   = "K" * string(K) * "_fVAR" * nfVARSpec * "_MDD" * nMDDSpec*"_MCMC" * nMCMCSpec
    LogFile    = open(savedir * sNameFile * "_MDDLogFile.txt", "w")

    Pen_Bayes   = zeros(size(MDD_GoF)[1],hyper_n)
    Pen_Laplace = zeros(size(MDD_GoF)[1],hyper_n)

    for ll=1:hyper_n

        time_init_loop = time_ns()

        lambda1 = hyper_grid[ll,1]
        lambda2 = hyper_grid[ll,2]
        lambda3 = hyper_grid[ll,3]

        if ii == 1
            MDD_Bayes_sum[ll,1:3]   = hyper_grid[ll,:]'
            MDD_Laplace_sum[ll,1:3] = hyper_grid[ll,:]'
            Pen_Bayes_sum[ll,1:3]   = hyper_grid[ll,:]'
            Pen_Laplace_sum[ll,1:3] = hyper_grid[ll,:]'
        end

        PHIpdraw, SIGMAtrpdraw, YYpred, lnpYY = bayesianVAR_wokron(data, Tdrop, nlags, nsim, nburn,
                                                n_agg, lambda1, lambda2, lambda3, tau, silent, seedoffset, iwishdf);
        # Estimation sample: SampleStart+Tdrop+1 to SampleEnd


        SIGMApmean = zeros(size(SIGMAtrpdraw)[2],size(SIGMAtrpdraw)[2])
            for nn = 1:size(SIGMAtrpdraw)[1]
                SIGMApmean = SIGMApmean + SIGMAtrpdraw[nn,:,:]*SIGMAtrpdraw[nn,:,:]'
            end
        SIGMApmean = 1/size(SIGMAtrpdraw)[1]*SIGMApmean
        SIGMAaa_inv = inv(SIGMApmean[n_agg+1:end, n_agg+1:end])

        tt_sel = 1
        # VinvLam_all contains Vinv for SamleStart:SampleEnd
        # However, in constructing YY, we drop the first Tdrop observations
        for tt = 1:size(MDD_GoF)[1]
            if tt > Tdrop
                Pen_Bayes[tt,ll]   = -0.5*logdet(VinvLam_all[:,:,tt]+SIGMAaa_inv) + Ktilde/2*log(2*pi)
                Pen_Laplace[tt,ll] = -0.5*logdet(VinvLam_all[:,:,tt]) + Ktilde/2*log(2*pi)
            end
        end

        Pen_Bayes_sum[ll,3+ii]   = sum(Pen_Bayes[:,ll]) + lnpYY
        Pen_Laplace_sum[ll,3+ii] = sum(Pen_Laplace[:,ll]) + lnpYY
        MDD_Bayes_sum[ll,3+ii]   = sum(MDD_GoF[Tdrop+1:end]) + Pen_Bayes_sum[ll,3+ii]
        MDD_Laplace_sum[ll,3+ii] = sum(MDD_GoF[Tdrop+1:end]) + Pen_Laplace_sum[ll,3+ii]

        print("Approximation order K $K \n")
        print("Compression Ktilde $n_cross \n")
        print("Step $ll out of $hyper_n \n")
        print("Hyperparameter lambda1 = $lambda1 \n")
        print("Hyperparameter lambda2 = $lambda2 \n")
        print("Hyperparameter lambda3 = $lambda3 \n")
        print("lnpYY = $(lnpYY) \n")
        print("Rel. Marginal data density  = $(MDD_Bayes_sum[ll,3+ii]-MDD_Bayes_sum[1,3+ii]) \n")
        time_loop=signed(time_ns()-time_init_loop)/1000000000
        print("Elapsed time           = $(time_loop) \n")
        print("================================= \n")

        write(LogFile,"Hyperparameter lambda1 = $lambda1 \n")
        write(LogFile,"Hyperparameter lambda2 = $lambda2 \n")
        write(LogFile,"Hyperparameter lambda3 = $lambda3 \n")
        write(LogFile,"lnpYY = $(lnpYY) \n")
        write(LogFile,"Rel. Marginal data density  = $(MDD_Bayes_sum[ll,3+ii]-MDD_Bayes_sum[1,3+ii]) \n")
        write(LogFile,"Elapsed time           = $(time_loop) \n")
        write(LogFile,"================================= \n")

    end

    close(LogFile)

    lambda_hat_K_BayesMax[ii,:]   = [ K n_cross MDD_Bayes_sum[argmax(MDD_Bayes_sum[:,3+ii]),[1 2 3 (3+ii)]] ]
    lambda_hat_K_LaplaceMax[ii,:] = [ K n_cross MDD_Laplace_sum[argmax(MDD_Laplace_sum[:,3+ii]),[1 2 3 (3+ii)]] ]

end


sNameFile = "fVAR" * nfVARSpec * "_MDD" * nMDDSpec*"_MCMC" * nMCMCSpec
CSV.write(savedir * sNameFile * "_MDD_Bayes_sum.csv", DataFrame(MDD_Bayes_sum))
CSV.write(savedir * sNameFile * "_MDD_Laplace_sum.csv", DataFrame(MDD_Laplace_sum))
CSV.write(savedir * sNameFile * "_Pen_Bayes_sum.csv", DataFrame(Pen_Bayes_sum))
CSV.write(savedir * sNameFile * "_Pen_Laplace_sum.csv", DataFrame(Pen_Laplace_sum))

sNameCols = ["K", "Ktilde", "Lam1", "Lam2", "Lam3", "MDD"]
CSV.write(savedir * sNameFile * "_lambda_hat_K_BayesMax.csv", DataFrame(permutedims(sNameCols)), header=false)
CSV.write(savedir * sNameFile * "_lambda_hat_K_BayesMax.csv", DataFrame(lambda_hat_K_BayesMax),append=true)
CSV.write(savedir * sNameFile * "_lambda_hat_K_LaplaceMax.csv", DataFrame(permutedims(sNameCols)), header=false)
CSV.write(savedir * sNameFile * "_lambda_hat_K_LaplaceMax.csv", DataFrame(lambda_hat_K_LaplaceMax),append=true)
