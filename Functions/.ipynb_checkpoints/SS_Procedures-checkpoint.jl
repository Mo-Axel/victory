function FFBS(A,B,H_t,R,Se,Phi1,y,n_agg,seedoffset,seedindx)

(T,l) = size(y)
(n,n) = size(Phi1);
s     = zeros(T+1,n);
P     = zeros(T+1,n,n);

# Initialization
s[1,:]   = zeros(n);
a        = inv(Matrix(I,n^2,n^2) - kron(Phi1,Phi1))*reshape(R*Se*R',n*n,1);
P[1,:,:] = reshape(a,n,n);

# Kalman Filter Recursion
sprime   = zeros(n,1);
liki     = ones(T,1);

# Forward Iterations
for t=1:T

    # Updating Step
    sprime = Phi1*s[t,:];
    Pprime = Phi1*P[t,:,:]*Phi1' + R*Se*R';

    # prediction step
    yprediction = A + B*sprime;
    v = y[t,:] - yprediction;
    F = B*Pprime*B' + H_t[:,:,t];
    if det(F) <= 0.0
        U,Sig,V = svd(F)
        detF = prod(Sig[Sig.>1e-12])
    else
        detF = det(F)
    end
    liki[t]    = -0.5*l*log(2*pi) - 0.5*log(detF) - 0.5*v'*inv(F)*v;

    # updating step
    kgain      = Pprime*B'*inv(F);
    s[t+1,:]   = (sprime + kgain*v);
    P[t+1,:,:] = Pprime - kgain*B*Pprime;

end

s = s[2:end,:]
P = P[2:end,:,:]

# initialize backward simulations for t=T
U,Sig,V       = svd(P[T,:,:])
Psqrt         = U*sqrt.(Diagonal(Sig))
s_smooth      = zeros(size(s));
s_smooth[T,:] = s[T,:]+(Psqrt*randn(size(s)[2],1))

# backward iterations
for t = 1:(T-1)

    st   = s[T-t,:]
    Pt   = P[T-t,:,:]
    Phat = Phi1*Pt*Phi1'+ + R*Se*R'
    Up,Sigp,Vp = svd(Phat);
    kp         = sum(Sigp.>1e-12);
    inv_Phat   = Symmetric(Up[:,1:kp]*Diagonal(1 ./ Sigp[1:kp])*Vp[:,1:kp]')

    vt    = s_smooth[T-t+1,:]- Phi1*st
    smean = st + Pt*Phi1'*inv_Phat*vt
    Pmean = Pt - Pt*Phi1'*inv_Phat*Phi1*Pt;
    Ufin, Sigfin, ~ = svd(Pmean);
    Pfinsqrt        = Ufin*sqrt.(Diagonal(Sigfin));

    Random.seed!(100000*seedindx+10*t+3+seedoffset)
    s_smooth[T-t,:] = smean+Pfinsqrt*randn(length(smean),1);

end

return s_smooth, liki

end
