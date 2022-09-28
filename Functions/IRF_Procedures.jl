function IRF_qSh(PHI, SIGMAtr, qstar, sh_size, Hmax, xgrid)

    YY_IRF         = zeros(Hmax+1,n_agg+n_cross)
    PhatDens_IRF   = zeros(Hmax+1,length(xgrid));
    Gini_IRF       = zeros(Hmax+1,1)
    # hh=1 is steady state
    # hh=2 is period of impact
    YY_IRF[1+1,:]  = SIGMAtr*qstar*sh_size
    if Hmax > 1
        for hh = 1+2:H+1
            YY_IRF[hh,:]  = YY_IRF[hh-1,:]'*PHI[:,:]
        end
    end
    # compute density IRFs
    for hh = 1:Hmax+1
        PhatDensCoef_hh    = coefRecover(YY_IRF[hh,n_agg+1:end]', PhatDensCoef_lambda, PhatDensCoef_mean)
        PhatDensNorm_hh    = lnpdfNormalize(PhatDensCoef_hh,knots, minimum(xgrid), maximum(xgrid))
        #PhatDensNorm_hh    = lnpdfNormalize_unrate_Rsum(PhatDensCoef_hh,knots,YY_IRF[hh,3]+mean_unrate, minimum(xgrid), maximum(xgrid))
        PhatDens_IRF[hh,:] = pdfEval(xgrid,PhatDensCoef_hh',knots,[PhatDensNorm_hh[1]]);
        Gini_IRF[hh]       = GiniCoef(PhatDensNorm_hh, xgrid)
    end
    
    return YY_IRF, PhatDens_IRF, Gini_IRF
end
    