function ProbMassDiff_Cutoff(PhatDens, densvalues_ss, xgrid)

    ngrid = length(xgrid);
    theta_sinh = 1;
    ygrid = 1/(2*theta_sinh)*(exp.(theta_sinh*xgrid) - exp.(-theta_sinh*xgrid));
    Jacobian = 1/2*(exp.(theta_sinh*xgrid) + exp.(-theta_sinh*xgrid));

    densvalues_diff = (PhatDens - densvalues_ss)./Jacobian;
    ygrid_diff      = ygrid[2:xn] - ygrid[1:xn-1];
    probmass_diff   = densvalues_diff[2:xn].*ygrid_diff;
    out = sum(probmass_diff[ygrid[2:xn]]);

    return out

end


function GiniCoef(PhatDens, xgrid)

    ngrid = length(xgrid);
    #theta_sinh = 1;
    #ygrid = 1/(2*theta_sinh)*(exp.(theta_sinh*xgrid) - exp.(-theta_sinh*xgrid));
    ygrid = 1*xgrid
    Jacobian = ones(size(xgrid));

    densvalues_pmean   = PhatDens/Jacobian;
    # compute probability mass function
    ygrid_diff        = ygrid;
    probmassfcn       = densvalues_pmean[1:ngrid].*ygrid_diff;
    probmassfcn[1]    = probmassfcn[1]+densvalues_pmean[1]; # for mixed distribution

    # compute Gini coefficient
    Silag = 0;
    G     = 0;
    Si    = 1;
    for ii = 1:(ngrid-1)
        Si = sum(probmassfcn[1:ii].*ygrid[2:(1+ii)]);
        G  = G+probmassfcn[ii]*(Silag+Si);
        Silag = Si;
    end

    out = 1-G/Si;

    return out

end

function GiniCoe(PhatDens, xgrid)

    ngrid = length(xgrid);
    theta_sinh = 1;
    ygrid = 1/(2*theta_sinh)*(exp.(theta_sinh*xgrid) - exp.(-theta_sinh*xgrid));
    Jacobian = 1/2*(exp.(theta_sinh*xgrid) + exp.(-theta_sinh*xgrid));

    densvalues_pmean   = PhatDens./Jacobian;
    # compute probability mass function
    ygrid_diff        = ygrid[2:ngrid] - ygrid[1:ngrid-1];
    probmassfcn       = densvalues_pmean[2:ngrid].*ygrid_diff;
    probmassfcn[1]    = probmassfcn[1]+densvalues_pmean[1]; # for mixed distribution

    # compute Gini coefficient
    Silag = 0;
    G     = 0;
    Si    = 1;
    for ii = 1:(ngrid-1)
        Si = sum(probmassfcn[1:ii].*ygrid[2:(1+ii)]);
        G  = G+probmassfcn[ii]*(Silag+Si);
        Silag = Si;
    end

    out = 1-G/Si;

    return out

end







function Gini(PhatDens, xgrid)
    
    
    densvalues_pmean   = PhatDens
    probmassfcn       = densvalues_pmean[1:ngrid].*xgrid;
    probmassfcn[1]    = probmassfcn[1]+densvalues_pmean[1]; # for mixed distribution
    
    ngrid     = 1000
    grid_temp = range(0.001, stop=1, length=ngrid);
    grid_temp = [0.0; grid_temp]
    vec_percs = range(0.001, stop=1, length=ngrid);
    emp_percs = DensPercentiles(probmassfcn[1], knots_all, xgrid, grid_temp, grid_temp)
    GINI = var(emp_percs)
    return GINI
end




var(range(0.001, stop=1, length=1000))







