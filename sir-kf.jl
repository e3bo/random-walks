using Distributions
using LinearAlgebra
using Optim

function obj(pvar::Vector; γ::Float64 = 365.25 / 9, dt::Float64 = 1 / 365.25, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 20e6, ρ::Float64 = 0.4, τ::Float64 = 0.01, beta_sd::Float64 = 0.1)
    # prior for time 0
    x0 = [19e6; pvar[1]; pvar[2]; 0]
    p0 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]

    dstate = size(x0, 1)

    # observation
    h = [0 0 0 ρ]
    dobs = size(h, 1)
    r = Matrix(undef, dobs, dobs)

    # data
    z = [[754], [665], [955], [908], [1005], [866], [834], [1189], [1000], [1382], [1598], [1731], [1222], [933], [1393], [1360], [1836], [1592], [1447], [1143]]

    # filter (assuming first observation at time 1)
    nobs = length(z)

    Σ = Array{eltype(pvar)}(undef, dobs, dobs, nobs)
    ytkkmo = Array{eltype(pvar)}(undef, dobs, nobs)
    k = Array{eltype(pvar)}(undef, dstate, nobs)
    xkk = Array{eltype(pvar)}(undef, dstate, nobs)
    xkkmo = Array{eltype(pvar)}(undef, dstate, nobs)
    pkk = Array{eltype(pvar)}(undef, dstate, dstate, nobs)
    pkkmo = Array{eltype(pvar)}(undef, dstate, dstate, nobs)
    
    bvec = pvar[3:22]
    
    for i in 1:nobs
        β = bvec[i]
        if (i == 1)
            xlast =  x0
            plast = p0
        else
            xlast = xkk[:,i - 1]
            plast = pkk[:,:,i-1]
        end
        
        x = xlast[1]
        l = xlast[2]
        y = xlast[3]
        xlast[4] = 0
        vf = [-β*x*y/N - ι*x, β*x*y/N + ι*x - η*l, η*l - γ*y, γ*y]
        xnext = xlast + dt * vf
        
        for j in 1:dstate
            if xnext[j] < 0
                xnext[j] = 0
            end
        end
        xkkmo[:,i] = xnext
        
        f = [0, β*x/N*y/N + ι*x/N, η*l/N, γ*y/N]
        
        q = [  f[1]+f[2]     -f[2]            0      0
                   -f[2] f[2]+f[3]        -f[3]      0
                       0     -f[3]    f[3]+f[4]  -f[4]
                       0         0        -f[4]   f[4]]
                       
        jac= [-β*y/N    0    -β*x/N     0
              β*y/N   -η     β*x/N     0
                  0    η        -γ     0
                  0    0         γ     0]
        
        dp = jac * plast + plast * jac' + q
        pkkmo[:,:,i] = plast + dp * dt
        
        r[1,1] = z[i][1] * τ 
        Σ[:,:,i] = h * pkkmo[:,:,i] * h' + r
        k[:,i] = pkkmo[:,:,i] * h' / Σ[:,:,i]
        ytkkmo[:,i] = z[i] - h * reshape(xkkmo[:,i], dstate, 1)
        xkk[:,i] = reshape(xkkmo[:,i], dstate, 1) + reshape(k[:,i], dstate, 1) * ytkkmo[:,i]
        pkk[:,:,i] = (I - reshape(k[:,i], dstate, 1) * h) * pkkmo[:,:,i]
    end
    
    jumpdensity = Normal(0, beta_sd)
    dbeta = [bvec[i] - bvec[i - 1] for i in 2:length(bvec)]
    rwlik = 0
    for diff in dbeta
        rwlik += logpdf(jumpdensity, diff)
    end
    
    nll = 0.5 * (sum(ytkkmo[1,:] .^2 ./ Σ[1,1,:] + map(log, Σ[1,1,:])) + nobs * log(2 * pi)) - rwlik
    #nll, ytkkmo, Σ, xkkmo
    nll
end


init = [1e3; 1e3; [42. for i in 1:20]]
lower = [1e2; 1e2; [0. for i in 1:20]]
upper = [1e5; 1e5; [420. for i in 1:20]]

ans0 = optimize(obj, lower, upper, init, Fminbox(BFGS())) # works, but no autodiff
#ans1 = optimize(obj, [100.0], LBFGS()) # works if init is close to minimizer
#ans2 = optimize(obj, [1.0], LBFGS(); autodiff = :forward) # reduced iterations
ans3 = optimize(obj, lower, upper, init, Fminbox(LBFGS()); autodiff = :forward) # can start far away 