# Specify parameters for prior
mutable struct priors
    mu_phi::Array{Float64,2}
    P_phi::Array{Float64,2}
    nu_iwish::Int64
    sig_iwish::Array{Float64,2}
end

##

function data_to_YY(data, T0, nlags)

    # Data handling
    Tn    = size(data)[1]
    YY    = data[T0+1:Tn,:]
    T     = size(YY)[1] # sample size

    XX    = data[T0:Tn-1,:]
    if nlags > 1
        for  ii = 2:nlags
            XX = [ XX data[T0+1-ii:Tn-ii,:]]
        end
    end

    # OLS estimator
    PHIhat   = \((XX'XX),(XX'YY))
    Shat     = (YY'*YY) - (YY'*XX)*PHIhat
    SIGMAhat = Shat/T

return YY, XX, PHIhat, Shat, SIGMAhat
end

##

function generate_prior(YY, PHIhat, SIGMAhat, n_agg, iwishdf, lambda1, lambda2, lambda3)

    T = size(YY)[1]
    n_cross = size(YY)[2]-n_agg

    mu_phi         = zeros(size(PHIhat))
    SIGMAhatinv    = inv(SIGMAhat)
    SIGMAhatinv_11 = SIGMAhatinv[1:n_agg, 1:n_agg]
    SIGMAhatinv_12 = SIGMAhatinv[1:n_agg, n_agg+1:end]
    SIGMAhatinv_22 = SIGMAhatinv[n_agg+1:end, n_agg+1:end]

    std_data = sum((YY-repeat(mean(YY,dims=1),T)).^2,dims=1)/T

    # changed diagm to Diagonal
    P_phi          = lambda1*[kron(SIGMAhatinv_11,Diagonal([ones(n_agg).*std_data[1:n_agg]; lambda2*ones(n_cross).*std_data[n_agg+1:end]])) kron(SIGMAhatinv_12,Diagonal([sqrt(lambda3)*ones(n_agg).*std_data[1:n_agg]; sqrt(lambda2)*ones(n_cross).*std_data[n_agg+1:end]])); kron(SIGMAhatinv_12',Diagonal([sqrt(lambda3)*ones(n_agg).*std_data[1:n_agg]; sqrt(lambda2)*ones(n_cross).*std_data[n_agg+1:end]])) kron(SIGMAhatinv_22,Diagonal([lambda3*ones(n_agg).*std_data[1:n_agg]; ones(n_cross).*std_data[n_agg+1:end]]))]

    # Prior for Sigma
    nu_iwish    = size(YY)[2]+iwishdf

    SIGMAdiag = zeros(size(SIGMAhat))
    for nn = 1:size(SIGMAdiag)[2]
        Ydiag = YY[2:end,nn]
        Xdiag = YY[1:end-1,nn]
        PHIdiag = \((Xdiag'Xdiag),(Xdiag'Ydiag))
        Sdiag     = (Ydiag'*Ydiag) - (Ydiag'*Xdiag)*PHIdiag
        SIGMAdiag[nn,nn] = Sdiag/length(Ydiag)
    end

    #sig_iwish   = nu_iwish*SIGMAhat
    sig_iwish   = nu_iwish*SIGMAdiag # use diagonal matrix
    prior       = priors(mu_phi, P_phi, nu_iwish, sig_iwish);

return prior
end

##

function bayesianVAR_j_wokron(data, T0, nlags, n_agg, prior, SIGMAj, seedoffset, seedindx)

YY, XX, PHIhat, Shat, SIGMAhat = data_to_YY(data, T0, nlags)
T  = size(YY)[1]
ny = size(YY)[2]

# Draw from Phi|(Sigma,YY)
kronSIGXX = kron(inv(SIGMAj),XX'*XX)
P_phiup   = Symmetric(prior.P_phi + kronSIGXX)
var_phiup = inv(P_phiup)
mu_phiup  = var_phiup*(prior.P_phi*vec(prior.mu_phi) + kronSIGXX*vec(PHIhat))
Random.seed!(10*seedindx+1+seedoffset)
PHIj      = rand(MvNormal(vec(mu_phiup), var_phiup))
PHIj      = reshape(PHIj, ny*nlags, ny)

# Draw from Sigma|(Phi,YY)
Sup       = Symmetric(prior.sig_iwish + (YY-XX*PHIj)'*(YY-XX*PHIj))
nup       = prior.nu_iwish + T
Random.seed!(10*seedindx+2+seedoffset)
SIGMAj    = rand(InverseWishart(nup,cholesky(Sup)))

return PHIj, SIGMAj

end

##

function bayesianVAR_wokron(data, T0, nlags, nsim, nburn, n_agg, lambda1, lambda2, lambda3, tau, silent, seedoffset, iwishdf)

YY, XX, PHIhat, Shat, SIGMAhat = data_to_YY(data, T0, nlags)
prior = generate_prior(YY, PHIhat, SIGMAhat, n_agg, iwishdf, lambda1, lambda2, lambda3)
T  = size(YY)[1]
ny = size(YY)[2]

# Normalization constant
ii = 1
consgam_low = 0
while ii<=ny
   #consgam_low = consgam_low + (logabsgamma(0.5*(prior.nu_iwish+1-ii)))[1]
   consgam_low = consgam_low + (lgamma(0.5*(prior.nu_iwish+1-ii)))[1]
   ii = ii+1
end
normalize_const = 0.5*prior.nu_iwish*logdet(prior.sig_iwish)-0.5*(ny*prior.nu_iwish)*log(2) - consgam_low

# Initialize SIGMAjj for Gibbs sampler
SIGMAjj = SIGMAhat

# Initialize matrices for collecting draws from posterior
PHIpdraw     = zeros(nsim, ny*nlags, ny);
SIGMAtrpdraw = zeros(nsim, ny, ny);
lnpYSIGMA    = zeros(nsim);

# Gibbs Sampling
counter = 0
if silent == 0
   println("")
   println(" Bayesian estimation of VAR: Direct Sampling... ")
   println("")
end

for jj = 1:(nsim+nburn)

    # Draw from Phi|(Sigma,YY)
    kronSIGXX = kron(inv(SIGMAjj),XX'*XX)
    P_phiup   = Symmetric(prior.P_phi + kronSIGXX)
    var_phiup = inv(P_phiup)
    mu_phiup  = var_phiup*(prior.P_phi*vec(prior.mu_phi) + kronSIGXX*vec(PHIhat))
    #mu_phiup  = inv(P_phiup)*(prior.P_phi*vec(prior.mu_phi) + kronSIGXX*vec(PHIhat))
    #mu_phiup  = \(P_phiup,(prior.P_phi*vec(prior.mu_phi) + kronSIGXX*vec(PHIhat)))
    #var_phiup = inv(P_phiup)
    Random.seed!(10*jj+1+seedoffset)
    PHIjj     = rand(MvNormal(vec(mu_phiup), var_phiup))
    PHIjj     = reshape(PHIjj, ny*nlags, ny)

    # Draw from Sigma|(Phi,YY)
    Sup     = Symmetric(prior.sig_iwish + (YY-XX*PHIjj)'*(YY-XX*PHIjj))
    nup     = prior.nu_iwish + T
    Random.seed!(10*jj+2+seedoffset)
    SIGMAjj = rand(InverseWishart(nup,cholesky(Sup)))

    if jj > nburn

       PHIpdraw[jj-nburn,:,:]     = PHIjj
       SIGMAtrpdraw[jj-nburn,:,:] = Matrix(sparse(cholesky(Symmetric(SIGMAjj)).L))

       lnpSIGMA            = normalize_const - 0.5*(prior.nu_iwish + ny + 1)*logdet(SIGMAjj)-0.5*tr(inv(SIGMAjj)*prior.sig_iwish)
       lnpYSIGMA[jj-nburn] = lnpSIGMA - 0.5*(ny*T)*log(2*pi) -0.5*T*logdet(SIGMAjj) + 0.5*logdet(prior.P_phi) - 0.5*logdet(P_phiup) - 0.5*tr(inv(SIGMAjj)*Shat) - 0.5*(
                             vec(PHIhat)'*kronSIGXX*vec(PHIhat) + (vec(prior.mu_phi)'*prior.P_phi*vec(prior.mu_phi)) - (mu_phiup'*P_phiup*mu_phiup))

    end

    # output
    counter = counter + 1
    if counter == 10000 & silent == 0
       println("")
       println(" Draw number:  $jj")
       println(" Remaining draws:  $(nsim+nburn-jj)")
       println(" lnpYSIGMA: $(lnpYSIGMA[jj-nburn])")
       println("")
       counter = 0
    end

end

PHIpmean     = dropdims(mean(PHIpdraw,dims=1),dims=1)
SIGMAtrpmean = dropdims(mean(SIGMAtrpdraw,dims=1),dims=1)
YYpred       = XX*PHIpmean;

# compute MDD
SIGMA_vech = zeros(nsim, Int(ny*(ny+1)/2), 1)
for jj = 1:nsim
    SIGMA_vech[jj,:,1]=vech(SIGMAtrpdraw[jj,:,:])
end
lnpYY = GewekeMDD(SIGMA_vech,lnpYSIGMA,tau)

return PHIpdraw, SIGMAtrpdraw, YYpred, lnpYY

end

##

function GewekeMDD(Bpdraw,lnpYB,tau)

ndraws = size(lnpYB)[1]
Bpdraw = reshape(Bpdraw, (size(Bpdraw)[1], size(Bpdraw)[2]*size(Bpdraw)[3], 1))
B_vec  = dropdims(Bpdraw,dims=3) # make each Bpdraw into vector
B_bar  = mean(B_vec,dims=1)'
npara  = length(B_bar)
VB_bar = (B_vec - repeat(mean(B_vec,dims=1),nsim))'*(B_vec-repeat(mean(B_vec,dims=1),nsim))/ndraws
invVB_bar = inv(VB_bar)
logdetVB_bar = logdet(VB_bar)
lnfB   = zeros(ndraws)

for jj = 1:ndraws
   temp_quad = ((B_vec[jj,:]-B_bar)'*invVB_bar*(B_vec[jj,:]-B_bar))[1]
   lnfB[jj]  = ( (-log(tau)-(npara/2)*log(2*pi)-0.5*logdetVB_bar
                -0.5*temp_quad)*( temp_quad <= cquantile(Chisq(npara),tau)) )
end

lnfpYB    = lnfB - lnpYB
maxlnfpYB = maximum(lnfpYB)

fB_YB = exp.(lnfpYB .- ones(ndraws)*maxlnfpYB)
lnpYY =  - ( log( sum(fB_YB)/ndraws ) + maxlnfpYB)

return lnpYY

end
