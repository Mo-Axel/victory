# function basis_logspline(x,knots)
#     # this procedure evaluates the spline basis functions
#     # at the n*1 vector x, given the (K-1)*1 vector of knots
#     # output is  n*K
#
#     x_dim  = length(x)
#     if x_dim == 1
#         x=[x]
#     end
#
#     basis_fcn = copy(x)
#
#     # cubic part
#     for i=1:length(knots)
#         basis_fcn = [basis_fcn  max.((x .- ones(x_dim)*knots[i]),zeros(x_dim)).^3]
#     end
#     return basis_fcn
# end

function pdfEval_noNorm(x,coef,knots)
    # this procedure evaluates the unnormalized exp(log spline density)
    # x is n*1
    # knots is (K-1)*1
    # basis is n*(K+1)
    # coef is (k+1)*1
    pdf_logspline = exp.(basis_logspline(x,knots)*coef)
    return pdf_logspline
end

function pdfEval(x,coef,knots,lnnorm)
    # this procedure evaluates the exp(1/sqrt(K)*log spline density)
    pdf_logspline = exp.(basis_logspline(x,knots)*coef - lnnorm.*ones(Float64, size(x)[1]) )
    return pdf_logspline
end

function logspline_obj(data, coef, knots, lb, ub)
    (IntPdf,error) = quadgk(x -> pdfEval_noNorm(x,coef,knots),lb,ub)
    out = (mean(basis_logspline(data,knots),dims=1)*coef)[1]-log(IntPdf[1])
return -out
end

function logspline_obj_topcode(data, coef, knots, lb, ub)
    C_topcode = maximum(data)
    N = length(data)
    N_max = sum(data.==C_topcode)
    if N_max == 1
        # no top coding in data because only one obs = max
        (IntPdf,error) = quadgk(x -> pdfEval_noNorm(x,coef,knots),lb,ub)
        N_max_factor = 1
    else
        # top coding is relevant because multiple obs = max
        # note that default order is 7
        (IntPdf,error) = quadgk(x -> pdfEval_noNorm(x,coef,knots),lb,C_topcode,order=10)
        N_max_factor = (N-N_max)/N
    end
    # out = (mean(basis_logspline(data[data.<C_topcode],knots),dims=1)*coef)[1]-N_max_factor*log(IntPdf[1])
    out = (sum(basis_logspline(data[data.<C_topcode],knots),dims=1)*coef)[1]/N-N_max_factor*log(IntPdf[1])
return -out
end


function lnpdfNormalize(coef,knots,lb=0,ub=10)
    # this procedure generates the normalization constant for the density
    # coef is T*K
    # knots is (K-1)*1
    lnnorm_out = zeros(size(coef)[1])
    for i = 1:size(coef,1)
       (IntPdf,error) = quadgk(x -> pdfEval_noNorm(x,coef[i,:],knots),lb,ub)
       lnnorm_out[i]  = log(IntPdf[1,1])
    end
    return lnnorm_out
end


function lnpdfNormalize_unrate(coef,knots,unrate,lb=0,ub=10)
    # this procedure generates the normalization constant for the density
    # coef is T*K
    # knots is (K-1)*1
    lnnorm_out = zeros(size(coef)[1])
    for i = 1:size(coef,1)
       (IntPdf,error) = quadgk(x -> pdfEval_noNorm(x,coef[i,:],knots),lb,ub)
       lnnorm_out[i]  = log(IntPdf[1,1]/(1-unrate[i]))
    end
    return lnnorm_out
end

function lnpdfNormalize_unrate_Rsum(coef,knots,unrate,lb=0,ub=10)
    # this procedure generates the normalization constant for the density
    # coef is T*K
    # knots is (K-1)*1

    xs = range(lb, stop=ub, length=1000)
    delta_x = xs[2]-xs[1]

    lnnorm_out = zeros(size(coef)[1])
    for i = 1:size(coef,1)
        IntPdf = 0.0
        for ii = 2:length(xs)
            IntPdf = IntPdf + (pdfEval_noNorm(xs[ii],coef[i,:],knots)*delta_x)[1]
        end
        lnnorm_out[i]  = log(IntPdf[1,1]/(1-unrate[i]))
    end

return lnnorm_out

end


function hessian_loglh(coef, knots, lb, ub)

    K = length(knots)+1
    (IntPdf,error) = quadgk(x -> pdfEval_noNorm(x,coef,knots),lb,ub)

    hess_l_jk = zeros(K,K)

    for j = 1:K
        for k = j:K
            (Int_j,error) = quadgk(x -> basis_logspline(x,knots)[j]*pdfEval_noNorm(x,coef,knots),lb,ub)
            Int_j_norm    = Int_j[1]/IntPdf[1]
            (Int_k,error) = quadgk(x -> basis_logspline(x,knots)[k]*pdfEval_noNorm(x,coef,knots),lb,ub)
            Int_k_norm    = Int_k[1]/IntPdf[1]
            (Int_jk,error) = quadgk(x -> (basis_logspline(x,knots)[j] - Int_j_norm[1])*(basis_logspline(x,knots)[k] - Int_k_norm[1])*pdfEval_noNorm(x,coef,knots),lb,ub)
            hess_l_jk[j,k] = -Int_jk[1]/IntPdf[1]
        end
    end

    hess_sym = Symmetric(hess_l_jk)

return hess_sym

end


function coefCompress(coef)
    # this procedure demeans and reduces the dimension of
    # the coefficients
    # coef is T*(k+2)
    # demean coefficients
    coef_mean   = mean(coef,dims=1)
    coef_normal = coef - ones(size(coef,1))*coef_mean
    coef_normal_cov = coef_normal'coef_normal/size(coef,1)

    # compute eigenvalues and eigenvectors
#    (D, V) = eig(coef_normal_cov) # previous verion of LinearAlgebra
    D = eigvals(coef_normal_cov) # version 1.5.
    V = eigvecs(coef_normal_cov)
    # Eliminate eigenvectors associated with small eigenvalues
    cutoff      = 1E-10
    #cutoff      = 0.1

    select      = D .> cutoff
    Dselect     = D[select]
    Vselect     = V[:,select]
    coef_factor = coef_normal*Vselect
    # standardization is new
    coef_factor_std = std(coef_factor,corrected=false, dims=1)
    coef_factor = coef_factor ./ coef_factor_std
    # end standardization

    coef_lambda = \((coef_factor'*coef_factor),(coef_factor'*coef_normal))
    # return results
    return (coef_factor, coef_lambda, coef_mean)
end

function coefRecover(coef_factor, coef_lambda, coef_mean)
    # this procedure recovers the actual from the compressed coefficients
    # coef_factor is m*n_f
    coef = coef_factor*coef_lambda+ones(size(coef_factor,1))*coef_mean
    return coef
end

function true_pdf(a,measureCoef_2_1t,measureCoef_2_2t,measureCoef_2_3t,moment_2_1t,moment_2_2t,moment_2_3t, aggEmployment)
    employed_g = exp(measureCoef_2_1t*(a - moment_2_1t) + measureCoef_2_2t*((a-moment_2_1t)^2 - moment_2_2t) + measureCoef_2_3t*((a-moment_2_1t)^3 - moment_2_3t))
    integrand = employed_g
    return integrand
end

function pdfWinberry(PDensCoef, moment_2_1, moment_2_2, moment_2_3, xgrid, aggEmployment)
    Tn = size(PDensCoef)[1]
    normalization_cons = zeros(Tn,1)
    for tt = 1:Tn
        (IntPdf,error) = quadgk(a -> true_pdf(a,PDensCoef[tt,1],PDensCoef[tt,2],PDensCoef[tt,3],moment_2_1[tt],moment_2_2[tt],moment_2_3[tt], aggEmployment), minimum(xgrid), maximum(xgrid))
        normalization_cons[tt] = IntPdf
    end

    PDens = zeros(Tn, length(xgrid))
    for tt = 1:Tn
        PDens[tt,:] = exp.(PDensCoef[tt,1]*(xgrid .- moment_2_1[tt]) + PDensCoef[tt,2]*((xgrid.-moment_2_1[tt]).^2 .- moment_2_2[tt]) + PDensCoef[tt,3]*((xgrid.-moment_2_1[tt]).^3 .- moment_2_3[tt]))
        PDens[tt,:] = PDens[tt,:]/normalization_cons[tt]
    end
    return PDens
end

function seasonality_adj(data_cross, nq, start_season, T)
# take the original coefficients and do seasonality adjustment

    # Generate seasonal dummies
    #start_season = 1 # start from 1989Q1 (first quarter)
    DD_seas = repeat(Matrix(I,nq,nq),Int(ceil(T/nq)),1)[start_season:end,:]  # seasonal dummies

    if size(DD_seas)[1]>=T
        #trim at T?
        DD_seas = DD_seas[1:T,:]
    else
        #otherwise continue with another identity matrix,
        #trimmed to restrict sample size to T
        temp = Matrix(I,nq,nq)[1:(T-size(DD_seas)[1]),:]
        DD_seas = [DD_seas; temp]
    end

    # seasonality adjustment
    # nq = 4 # quarterly data
    n_cross = size(data_cross)[2]

    #BIC for seansonal model
    C_ols = \((DD_seas'*DD_seas),(DD_seas'*data_cross))
    CSIGMAhat       = (data_cross'*data_cross - data_cross'*DD_seas*C_ols)/T
    data_cross_seas = data_cross-DD_seas*C_ols

    BIC_seas = zeros(n_cross)
    for nc = 1:n_cross
        RSS_seas = sum((data_cross_seas[:,nc]).^2)
        BIC_seas[nc] = T + T*log(2*pi) + T*log(RSS_seas/T) + log(T)*(nq+1)
    end

    # BIC for model with intercept only
    DD_int = ones(T,1) # intercept, having the mean only
    C_int_ols = \((DD_int'*DD_int),(DD_int'*data_cross))
    data_cross_int = data_cross-DD_int*C_int_ols

    BIC_int = zeros(n_cross)
    for nc = 1:n_cross
        RSS_int = sum((data_cross_int[:,nc]).^2)
        BIC_int[nc] = T + T*log(2*pi) + T*log(RSS_int/T) + log(T)*(1+1)
    end

    seasonal_mean = zeros(size(data_cross))

    # pick the specification with lower BIC (mean only vs. seasonal dummies)
    for nc = 1:n_cross
        if BIC_seas[nc] < BIC_int[nc]
            seasonal_mean[:,nc] = (DD_seas*C_ols)[:,nc]
        else
            seasonal_mean[:,nc] = (DD_int*C_int_ols)[:,nc]
        end
    end

    data_cross = data_cross - seasonal_mean

    avg_seasonal_mean = mean(seasonal_mean[1:nq,:],dims=1)

return (data_cross, avg_seasonal_mean, seasonal_mean)

end
