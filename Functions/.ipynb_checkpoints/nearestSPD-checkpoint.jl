function nearestSPD(A)

    B = 0.5*(A+A')
    
    # compute the symmetric polar factor of B. Call it H. Clearly H is itself SPD.
    U,Sigma,V = svd(B)
    H = V*diagm(Sigma)*V'
    
    # get Ahat in the above formula
    Ahat = 0.5*(B+H)
    
    # ensure symmetry
    Ahat = 0.5*(Ahat+Ahat')
    
    # test that Ahat is in fact PD. If it is not so, then tweak it just a bit
    
    p = 1
    k = 0
    
    while p == 1
    if (isposdef(Ahat)) && (det(Ahat)>0.0)
        p = 0
    end
    
    k = k+1
    
    if  p == 1
        # Ahat fail the chol test. It must have been just a hair off,
        # due to floating point trash, so it is simplest now just to
        # tweak by adding a tinay multiple of an identity matrix
    
        #mineig = minimum(eig(Ahat)[1])
        #maxeig = maximum(eig(Ahat)[1])
        mineig = minimum(eigvals(Ahat))
        maxeig = maximum(eigvals(Ahat))
    
    #    Ahat = Ahat + (-mineig*k.^2 + eps(maxeig)*eye(size(A)[1]))
        Ahat = Ahat + (-mineig*k.^2 .+ eps(maxeig)*LinearAlgebra.I(size(A)[1]))
    end
    
    end
    
        return Ahat
    
    end
    