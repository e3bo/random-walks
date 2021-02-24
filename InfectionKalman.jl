module InfectionKalman

using DataFrames
using Distributions
using ForwardDiff
using LinearAlgebra
using Optim

export hess
export obj
export grad

function obj(pvar::Vector, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, ρ::Float64 = 0.4, a::Float64 = 1., betasd::Float64 = 1., just_nll::Bool = true)

    # prior for time 0
    l0 = exp(pvar[1])
    τ = exp(pvar[2])
    bvec = pvar[3:end]
    
    y0 = l0 * η / γ
    x0 = [N - l0 - y0; l0; y0; 0]
    p0 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]

    #println(pvar)
    dstate = size(x0, 1)

    # observation matrix
    h = [0 0 0 ρ]
    
    dobs = 1
    r = Matrix(undef, dobs, dobs)

    # filter (assuming first observation at time 1)
    nobs = length(z)

    Σ = Array{eltype(bvec)}(undef, dobs, dobs, nobs)
    ytkkmo = Array{eltype(bvec)}(undef, dobs, nobs)
    k = Array{eltype(bvec)}(undef, dstate, nobs)
    xkk = Array{eltype(bvec)}(undef, dstate, nobs)
    xkkmo = Array{eltype(bvec)}(undef, dstate, nobs)
    pkk = Array{eltype(bvec)}(undef, dstate, dstate, nobs)
    pkkmo = Array{eltype(bvec)}(undef, dstate, dstate, nobs)
    logβ = Array{eltype(bvec)}(undef, nobs)

    @assert length(bvec) == nobs "length of bvec should equal number of observations"

    logβ[nobs] = bvec[nobs]
    for i in (nobs - 1):-1:1
        logβ[i] = (logβ[i + 1] - log(γ) - bvec[i]) / a + log(γ)
    end

    for i in 1:nobs
        if i == 1
            xlast =  x0
            plast = p0
        else
            xlast = xkk[:,i - 1]
            plast = pkk[:,:,i-1]
        end
        β = exp(logβ[i])
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
        pnext = plast + dp * dt
        for j in 1:dstate
            if pnext[j, j] < 0
                pnext[j, :] .= 0
                pnext[:, j] .= 0
            end
        end
        pkkmo[:,:,i] = pnext 
        
        
        r[1,1] = τ
        Σ[:,:,i] = h * pkkmo[:,:,i] * h' + r
        k[:,i] = pkkmo[:,:,i] * h' / Σ[:,:,i]
        ytkkmo[:,i] = z[i] - h * reshape(xkkmo[:,i], dstate, 1)
        xkk[:,i] = reshape(xkkmo[:,i], dstate, 1) + reshape(k[:,i], dstate, 1) * ytkkmo[:,i]
        pkk[:,:,i] = (I - reshape(k[:,i], dstate, 1) * h) * pkkmo[:,:,i]
    end
    
    stepdensity = Normal(0, betasd)
    rwlik = 0
    for bval in bvec[1:(end-1)]
        rwlik += logpdf(stepdensity, bval)
    end
    nll = 0.5 * (sum(ytkkmo[1,:] .^2 ./ Σ[1,1,:] + map(log, Σ[1,1,:])) + nobs * log(2 * pi)) - rwlik
    if just_nll
       return nll
    else
       return nll, ytkkmo, Σ, xkkmo, pkkmo, pkk
    end
end

function grad(pvar, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, ρ::Float64 = 0.4, a::Float64 = 1., betasd::Float64 = 1.)
    g = ForwardDiff.gradient(par -> obj(par, z; ρ = ρ, N = N, η = η, γ = γ, a = a, betasd = betasd), pvar)
    g
end

function hess(logE0::Float64, logtau::Float64, bpars, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, ρ::Float64 = 0.4, a::Float64 = 1., betasd::Float64 = 1.)
    pvar = [logtau bpars[end]]
    h = ForwardDiff.hessian(par -> obj(vec([logE0 par[1] bpars[1:(end -1)]' par[2]]), z; ρ = ρ, N = N, η = η, γ = γ, a = a, betasd = betasd), pvar)
    h
end

function fit(cdata, pdata; detailed_results::Bool = false, hessian::Bool = false, time_limit = 600, show_trace::Bool = false, betasd::Float64 = 1., N::Float64 = 1e7, a::Float64 = 1.) 

    wsize = size(pdata)[1] - 2
    z = [[el] for el in cdata.smooth[end-wsize+1:end]]
    w = [el for el in cdata.wday[end-wsize+1:end]]

    res = optimize(pvar -> obj(pvar, z, w; betasd = betasd, N = N, a = a), pdata.lower, pdata.upper, pdata.init, Fminbox(LBFGS()), Optim.Options(show_trace = show_trace, time_limit = time_limit); autodiff = :forward)

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