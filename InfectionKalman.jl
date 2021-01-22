module InfectionKalman

using DataFrames
using Distributions
using ForwardDiff
using LinearAlgebra
using Optim

export fit

function obj(pvar::Vector, z, w; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 20e6, ρ1::Float64 = 0.4, just_nll::Bool = true, betasd::Float64 = 1.)
    # prior for time 0
    x0 = [19e6; pvar[1]; pvar[2]; 0]
    p0 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]
    
    τ = pvar[3]
    ρ2 = pvar[4]
    
    dstate = size(x0, 1)

    # cyclic observation matrix
    hm = [0 0 0 ρ2
         0 0 0 ρ2
         0 0 0 ρ2
         0 0 0 ρ1
         0 0 0 ρ1
         0 0 0 ρ1
         0 0 0 ρ1]
    
    dobs = 1
    r = Matrix(undef, dobs, dobs)

    # filter (assuming first observation at time 1)
    nobs = length(z)

    Σ = Array{eltype(pvar)}(undef, dobs, dobs, nobs)
    ytkkmo = Array{eltype(pvar)}(undef, dobs, nobs)
    k = Array{eltype(pvar)}(undef, dstate, nobs)
    xkk = Array{eltype(pvar)}(undef, dstate, nobs)
    xkkmo = Array{eltype(pvar)}(undef, dstate, nobs)
    pkk = Array{eltype(pvar)}(undef, dstate, dstate, nobs)
    pkkmo = Array{eltype(pvar)}(undef, dstate, dstate, nobs)
    
    bvec = pvar[5:end]
    @assert length(bvec) == nobs "length of bvec should equal number of observations"
    
    for i in 1:nobs
        h = reshape(hm[w[i],:], dobs, dstate)
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
        plast[4,:] .= 0
        plast[:,4] .= 0
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
        
        dp = jac * plast + plast * jac' + q * N
        pkkmo[:,:,i] = plast + dp * dt
        
        r[1,1] = z[i][1] * τ 
        Σ[:,:,i] = h * pkkmo[:,:,i] * h' + r
        k[:,i] = pkkmo[:,:,i] * h' / Σ[:,:,i]
        ytkkmo[:,i] = z[i] - h * reshape(xkkmo[:,i], dstate, 1)
        xkk[:,i] = reshape(xkkmo[:,i], dstate, 1) + reshape(k[:,i], dstate, 1) * ytkkmo[:,i]
        pkk[:,:,i] = (I - reshape(k[:,i], dstate, 1) * h) * pkkmo[:,:,i]
    end
    
    jumpdensity = Normal(0, betasd)
    dbeta = [bvec[i] - bvec[i - 1] for i in 2:length(bvec)]
    rwlik = 0
    for diff in dbeta
        rwlik += logpdf(jumpdensity, diff)
    end
    
    nll = 0.5 * (sum(ytkkmo[1,:] .^2 ./ Σ[1,1,:] + map(log, Σ[1,1,:])) + nobs * log(2 * pi)) - rwlik
    if just_nll
       return nll
    else
       return nll, ytkkmo, Σ, xkkmo, pkkmo, pkk
    end
end

function hess(par, z, w)
    @time h = ForwardDiff.hessian(pvar -> obj(pvar, z, w), par)
    h
end

function fit(cdata, pdata; detailed_results::Bool = false, hessian::Bool = false, time_limit = 600, show_trace::Bool = false) 

    wsize = size(pdata)[1] - 4
    z = [[el] for el in cdata.reports[end-wsize+1:end]]
    w = [el for el in cdata.wday[end-wsize+1:end]]

    res = optimize(pvar -> obj(pvar, z, w), pdata.lower, pdata.upper, pdata.init, Fminbox(LBFGS()), Optim.Options(show_trace = show_trace, time_limit = time_limit); autodiff = :forward)

    if hessian 
        h = hess(res.minimizer, z, w)
        res = [res, h]
    end
    
    if detailed_results
        n, r, s, x, pk, pkk = obj(res.minimizer, z, w; just_nll = false)
        res = [res, n, r, s, x, pk, pkk]
    end
    res
end

end