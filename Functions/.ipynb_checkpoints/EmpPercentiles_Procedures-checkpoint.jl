function phatEval(x,PhatDensCoef, knots, xgrid)
    # compute normalization constant of cross-sectional density, density normalizes to 1-unempl.
    #PhatDensCoef_tt    = coefRecover(statepmean_tt',PhatDensCoef_lambda, PhatDensCoef_mean)
#    PhatDensCoef_tt    = lnpdfNormalize_unrate(PhatDensCoef_tt, knots, (statepmean_t[3] + mean(UNR)),0,3)
    PhatDensNorm    = lnpdfNormalize(PhatDensCoef, knots, minimum(xgrid),maximum(xgrid))
    #PhatDensNorm    = lnpdfNormalize_unrate_Rsum(PhatDensCoef, knots, unrate, minimum(xgrid),maximum(xgrid))
    out             = exp.((basis_logspline(x,knots)*PhatDensCoef')[1] - PhatDensNorm[1])
    return out
end

function phatIntegrate(PhatDensCoef, knots, xgrid, lb, ub)
    # integrate normalized density from lb to ub. Here normalized density integrates to 1-u
    (Intout,error) = quadgk(x -> phatEval(x,PhatDensCoef, knots, xgrid), lb, ub)
    return Intout
end

function phatIntegrate_lmass(PhatDensCoef, knots, xgrid, lb, ub)
    # integrate normalized density from lb to ub. Add unemployment rate to integral so that overall density would normalize to one
    (Intout,error) = quadgk(x -> phatEval(x,PhatDensCoef, knots, xgrid),lb, ub)
    Intout = Intout
    return Intout
end

function DensPercentiles(PhatDensCoef, knots, xgrid, grid_temp, vec_percs)

    valuesIntegrate = zeros(length(grid_temp)-1,1)

    # iterate over grid
    for nn = 2:length(grid_temp)
        # nn=2 corresponds to the first non-zero grid value
        if nn  == 2
            # add unemployment rate to cdf (mass of individuals earning zero)
            valuesIntegrate[nn-1] = phatIntegrate_lmass(PhatDensCoef, knots, xgrid, grid_temp[nn-1], grid_temp[nn])[1]
        else
            valuesIntegrate[nn-1] = phatIntegrate(PhatDensCoef, knots, xgrid, grid_temp[nn-1], grid_temp[nn])[1]
        end
    end

    emp_cdf = cumsum(valuesIntegrate,dims=1)
    emp_percs = zeros(length(vec_percs))

    itp = LinearInterpolation(vec(emp_cdf), vec(grid_temp[2:end]), extrapolation_bc = Line())
    for pp = 1:length(vec_percs)
        emp_percs[pp] = itp(vec_percs[pp])
    end

    #theta = 1.0
    #emp_percs = 1/(2*theta)*(exp.(theta*emp_percs) - exp.(-theta*emp_percs))

    return emp_percs

end
