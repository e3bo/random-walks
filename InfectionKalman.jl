module InfectionKalman

using Distributions
using ForwardDiff
using LinearAlgebra

export hess
export obj
export grad

function obj(pvar::Vector, z; γ::Float64 = 365.25 / 9,  γd::Float64 = 365.25 / 10,  γh::Float64 = 365.25 / 10,   dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, ρ::Float64 = 0.4, a::Float64 = 1., betasd::Float64 = 1., just_nll::Bool = true,
maxlogRt::Float64 = 1.6)

    # prior for time 0
    l0 = exp(pvar[1])
    h0 = exp(pvar[2])
    d0 = h0 * γh / γd
    
    τc = exp(pvar[3])
    τh = exp(pvar[4])
    τd = exp(pvar[5])
    
    chp = exp(pvar[6])
    hfp = exp(pvar[7])
    bvec = pvar[8:end]
    
    y0 = l0 * η / γ
    
    x0 = [N - l0 - y0 - h0 - d0; l0; y0; 0; 0; h0; d0; 0]
    p0 = convert(Array{eltype(bvec), 2}, Diagonal([1, 1, 1, 0, 0, 1, 1, 0]))
    
    dstate = size(x0, 1)

    # observation matrix
    hmat = [0 0 0 ρ 1 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1]
    
    dobs = size(z, 2)
    r = Diagonal([τc, τh, τd]) 

    zmiss = [ismissing(x) for x in z]
    zz = Array{eltype(z)}(undef, dobs)
    maxzscore = Inf

    # filter (assuming first observation at time 1)
    nobs = size(z, 1)
    rdiagadj = Array{eltype(bvec)}(undef, dobs, nobs)
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
        logβ[i] = min((logβ[i + 1] - log(γ) - bvec[i]) / a + log(γ), maxlogRt + log(γ))
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
        h = xlast[6]
        d = xlast[7]
        
        xlast[4] = 0
        xlast[5] = 0
        xlast[8] = 0
        plast[4,:] .= 0
        plast[:,4] .= 0
        plast[5,:] .= 0
        plast[:,5] .= 0
        plast[8,:] .= 0
        plast[:,8] .= 0
        
        vf = [
        -β*x*y/N - ι*x, 
        β*x*y/N + ι*x - η*l, 
        η*l - γ*y, 
        (1 - chp) * γ*y, 
        chp*γ*y,
        chp*γ*y - γh*h,
        hfp*γh*h - γd*d,
        γd*d
        ]
        xnext = xlast + dt * vf
        for j in 1:dstate
            if xnext[j] < 0
                xnext[j] = 0
            end
        end
        xkkmo[:,i] = xnext
        
        f = [0, 
             β*x/N*y/N + ι*x/N, 
             η*l/N, 
             (1 - chp)*γ*y/N,
             chp*γ*y/N,
             (1 - hfp) * γh * h / N,
             hfp * γh * h / N,
             γd * d / N
        ]
        
        q = [  f[1]+f[2]     -f[2]               0      0        0              0          0      0
                   -f[2] f[2]+f[3]           -f[3]      0        0              0          0      0
                       0     -f[3]  f[3]+f[4]+f[5]  -f[4]    -f[5]           -f[5]         0      0
                       0         0           -f[4]   f[4]        0              0          0      0
                       0         0           -f[5]      0     f[5]              0          0      0
                       0         0           -f[5]      0        0  f[7]+f[5]+f[6]      -f[7]     0
                       0         0              0       0        0           -f[7]  f[8]+f[7]  -f[8]
                       0         0              0       0        0              0       -f[8]   f[8] ]
                       
        jac= [-β*y/N    0         -β*x/N     0     0       0    0 0
               β*y/N   -η          β*x/N     0     0       0    0 0
                   0    η             -γ     0     0       0    0 0
                   0    0    (1 - chp)*γ     0     0       0    0 0
                   0    0          chp*γ     0     0       0    0 0
                   0    0          chp*γ     0     0     -γh    0 0  
                   0    0              0     0     0  hfp*γh  -γd 0
                   0    0              0     0     0       0   γd 0]
        
        dp = jac * plast + plast * jac' + q * N
        pnext = plast + dp * dt
        for j in 1:dstate
            if pnext[j, j] < 0
                pnext[j, :] .= 0
                pnext[:, j] .= 0
            end
        end
        pkkmo[:,:,i] = pnext 
        Σ[:,:,i] = hmat * pkkmo[:,:,i] * hmat' + r

        for j in 1:dobs
            if zmiss[i,j]
                zz[j] = 0
            else 
                zz[j] = z[i,j]
            end
        end       
        ytkkmo[:,i] = zz - hmat * reshape(xkkmo[:,i], dstate, 1)
        for j in 1:dobs
            if zmiss[i, j]
                zscore = 0
                rdiagadj[j,i] = 0
            else
                sd = sqrt(Σ[j,j,i])
                zscore = ytkkmo[j,i] / sd 
                if abs(zscore) > maxzscore
                    adjzscore = maxzscore / (1 + abs(zscore) - maxzscore)
                    newsd = abs(ytkkmo[j,i]) / adjzscore
                    rdiagadj[j,i] = newsd ^ 2 - sd ^ 2
                else 
                    rdiagadj[j,i] = 0
                end
            end
            Σ[j,j,i] += rdiagadj[j,i]
        end
        k[:,:,i] = pkkmo[:,:,i] * hmat' / Σ[:,:,i]
        for j in 1:dobs
            if zmiss[i,j]
                k[:, j ,i] .= 0
            end
        end
        xkk[:,i] = reshape(xkkmo[:,i], dstate) + reshape(k[:,:,i], dstate, dobs) * ytkkmo[:,i]
        pkk[:,:,i] = (I - reshape(k[:,:,i], dstate, dobs) * hmat) * pkkmo[:,:,i]
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
       return nll, ytkkmo, Σ, xkkmo, pkkmo, pkk, rdiagadj
    end
end

function grad(pvar, z; γ::Float64 = 365.25 / 9,  γd::Float64 = 365.25 / 10,  γh::Float64 = 365.25 / 10, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, ρ::Float64 = 0.4, a::Float64 = 1., betasd::Float64 = 1.)
    g = ForwardDiff.gradient(par -> obj(par, z; ρ = ρ, N = N, η = η, γ = γ,  γd = γd,  γh = γh, a = a, betasd = betasd), pvar)
    g
end

function hess(logE0::Float64, logtau::Float64, bpars, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, ρ::Float64 = 0.4, a::Float64 = 1., betasd::Float64 = 1.)
    pvar = [logtau bpars[end]]
    h = ForwardDiff.hessian(par -> obj(vec([logE0 par[1] bpars[1:(end -1)]' par[2]]), z; ρ = ρ, N = N, η = η, γ = γ, a = a, betasd = betasd), pvar)
    h
end

end