module InfectionKalman

using Distributions
using ForwardDiff
using LinearAlgebra

export hess
export obj
export grad

function obj(pvar::Vector, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, ρ::Float64 = 0.4, a::Float64 = 1., betasd::Float64 = 1., just_nll::Bool = true)

    # prior for time 0
    l0 = exp(pvar[1])
    τc = exp(pvar[2])
    τh = exp(pvar[3])
    chr = exp(pvar[4])
    bvec = pvar[5:end]
    
    y0 = l0 * η / γ
    x0 = [N - l0 - y0; l0; y0; 0; 0]
    p0 = [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 0]

    #println(pvar)
    dstate = size(x0, 1)

    # observation matrix
    h = [0 0 0 ρ 1; 0 0 0 0 1]
    
    dobs = size(z, 2)
    r = Matrix(undef, dobs, dobs)
    r[1,1] = τc
    r[1,2] = 0
    r[2,1] = 0
    r[2,2] = τh
    
    zmiss = [ismissing(x) for x in z]
    zz = Array{eltype(z)}(undef, dobs)

    # filter (assuming first observation at time 1)
    nobs = size(z, 1)

    Σ = Array{eltype(bvec)}(undef, dobs, dobs, nobs)
    ytkkmo = Array{eltype(bvec)}(undef, dobs, nobs)
    k = Array{eltype(bvec)}(undef, dstate, dobs, nobs)
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
        xlast[5] = 0
        plast[4,:] .= 0
        plast[:,4] .= 0
        plast[5,:] .= 0
        plast[:,5] .= 0
        vf = [-β*x*y/N - ι*x, β*x*y/N + ι*x - η*l, η*l - γ*y - chr*γ*y, γ*y, chr*γ*y]
        xnext = xlast + dt * vf
        for j in 1:dstate
            if xnext[j] < 0
                xnext[j] = 0
            end
        end
        xkkmo[:,i] = xnext
        
        f = [0, β*x/N*y/N + ι*x/N, η*l/N, γ*y/N, chr*γ*y/N]
        
        q = [  f[1]+f[2]     -f[2]               0      0        0
                   -f[2] f[2]+f[3]           -f[3]      0        0
                       0     -f[3]  f[3]+f[4]+f[5]  -f[4]    -f[5]
                       0         0           -f[4]   f[4]        0
                       0         0           -f[5]      0     f[5] ]
                       
        jac= [-β*y/N    0         -β*x/N     0     0
               β*y/N   -η          β*x/N     0     0
                   0    η   -γ*(1 + chr)     0     0
                   0    0              γ     0     0
                   0    0       chr *  γ     0     0]
        
        dp = jac * plast + plast * jac' + q * N
        pnext = plast + dp * dt
        for j in 1:dstate
            if pnext[j, j] < 0
                pnext[j, :] .= 0
                pnext[:, j] .= 0
            end
        end
        pkkmo[:,:,i] = pnext 
        
        Σ[:,:,i] = h * pkkmo[:,:,i] * h' + r
        k[:,:,i] = pkkmo[:,:,i] * h' / Σ[:,:,i]
        for j in 1:dobs
            if zmiss[i,j]
                k[:, j ,i] .= 0
                zz[j] = 0
            else 
                zz[j] = z[i,j]
            end
        end
        ytkkmo[:,i] = zz - h * reshape(xkkmo[:,i], dstate, 1)
        xkk[:,i] = reshape(xkkmo[:,i], dstate) + reshape(k[:,:,i], dstate, dobs) * ytkkmo[:,i]
        pkk[:,:,i] = (I - reshape(k[:,:,i], dstate, dobs) * h) * pkkmo[:,:,i]
    end
    
    stepdensity = Normal(0, betasd)
    rwlik = 0
    for bval in bvec[1:(end-1)]
        rwlik += logpdf(stepdensity, bval)
    end
    kflik = 0
    for i in 1:nobs
         sel = [!x for x in zmiss[i,:]]
         kflik -= 0.5 * (ytkkmo[sel,i]' / Σ[sel, sel, i] * ytkkmo[sel,i] +  log(det(Σ[sel, sel, i]))  + dobs * log(2 * pi))
    end
    nll = -kflik - rwlik
    
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

end