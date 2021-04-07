module InfectionKalman

using DataFrames
using DiffEqSensitivity
using Distributions
using ForwardDiff
using LinearAlgebra
using OrdinaryDiffEq

export hess
export obj
export grad
export vectorfield

function detect_frac(t; max_detect_par::Float64 = 0.4, detect_inc_rate::Float64 = 1.1, half_detect::Float64 = 30., base_detect_frac::Float64 = 0.1)
    rho =  max_detect_par * (t ^ detect_inc_rate)  / ((half_detect ^ detect_inc_rate) + (t ^ detect_inc_rate)) + base_detect_frac
    return rho
end

function hmat(t; t0 = 2020.164)
    day = (t - t0) * 365.25
    H =  [0 0 0 detect_frac(day) 1 0 0 0; 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1] 
    return H
end

function vectorfield(du, u, par, t)
  x = u[1,1]
  l = u[2,1]
  y = u[3,1]
  h = u[6,1]
  d = u[7,1]
  
  β = par[1]
  N = par[2]
  ι = par[3]
  η = par[4]
  γ = par[5]
  γd = par[6]
  γh = par[7]
  chp = par[8]
  hfp = par[9]
          
  du[1,1] = dx = -β*x*y/N - ι*x
  du[2,1] = dl = β*x*y/N + ι*x - η*l
  du[3,1] = dy = η*l - γ*y 
  du[4,1] = dc = (1 - chp) * γ*y
  du[5,1] = dhnew = chp*γ*y
  du[6,1] = dh = chp*γ*y - γh*h
  du[7,1] = dd = hfp*γh*h - γd*d
  du[8,1] = ddrep = γd*d

  f = [0, 
       β*x/N*y/N + ι*x/N, 
       η*l/N, 
       (1 - chp)*γ*y/N,
       chp*γ*y/N,
       (1 - hfp) * γh * h / N,
       hfp * γh * h / N,
       γd * d / N]

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
   
   p = u[:,2:end]
   du[:,2:end] .= jac * p + p * jac' + q * N
end

function obj(pvar::Vector, cov, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, a::Float64 = 1., betasd::Float64 = 1., just_nll::Bool = true, maxlogRt::Float64 = 1.6, γd::Float64 = 365.25 / 1, γh::Float64 = 365.25 / 1, h0::Float64 = 10., τh::Float64 = 10., τd::Float64 = 10., chp::Float64 = 0.01, hfp::Float64 = 0.01)

    zloc = deepcopy(z)
    if size(z, 2) == 3

    
        l0 = exp(pvar[1])
        h0 = exp(pvar[2])
        d0 = h0 * γh / γd
    
        τc = exp(pvar[3])
        τh = exp(pvar[4])
        τd = exp(pvar[5])
    
        chp = exp(pvar[6])
        hfp = exp(pvar[7])
        
        γh = exp(pvar[8])
        γd = exp(pvar[8])
        
        doseeffect = exp(pvar[9])
        prophomeeffect = exp(pvar[10])

        bvec = pvar[11:end]
    elseif size(z, 2) == 1 && "cases" in names(z)
        l0 = exp(pvar[1])
        τc = exp(pvar[2])
        doseeffect = exp(pvar[3])
        
        zloc[!, "hospitalizations"] .= missing
        zloc[!, "deaths"] .= missing
    
        bvec = pvar[4:end]
    end
    zloc = Matrix(select(zloc, :cases, :hospitalizations, :deaths)) # ensure assumed column order
    
    d0 = h0 * γh / γd
    y0 = l0 * η / γ
    
    x0 = [N - l0 - y0 - h0 - d0; l0; y0; 0; 0; h0; d0; 0]
    p0 = convert(Array{eltype(bvec), 2}, Diagonal([1, 1, 1, 0, 0, 1, 1, 0]))
    
    dstate = size(x0, 1)
    dobs = size(zloc, 2)
    r = Diagonal([τc, τh, τd]) 

    zmiss = [ismissing(x) for x in zloc]
    zz = Array{eltype(zloc)}(undef, dobs)
    maxzscore = Inf

    # filter (assuming first observation at time 1)
    nobs = size(zloc, 1)
    rdiagadj = ones(eltype(bvec), dobs, nobs) # make non-zero to ensure Σ is not singular
    Σ = Array{eltype(bvec)}(undef, dobs, dobs, nobs)
    ytkkmo = Array{eltype(bvec)}(undef, dobs, nobs)
    k = Array{eltype(bvec)}(undef, dstate, dobs, nobs)
    xkk = Array{eltype(bvec)}(undef, dstate, nobs)
    xkkmo = Array{eltype(bvec)}(undef, dstate, nobs)
    pkk = Array{eltype(bvec)}(undef, dstate, dstate, nobs)
    pkkmo = Array{eltype(bvec)}(undef, dstate, dstate, nobs)
    logβ = bvec

    @assert length(bvec) == nobs "length of bvec should equal number of observations"

    for i in 1:nobs
        if i == 1
            xlast = x0
            plast = p0
        else
            xlast = xkk[:,i - 1]
            plast = pkk[:,:,i-1]
        end
        β = exp(logβ[i] - cov.doses[i]*doseeffect - cov.prophome[i] * prophomeeffect)
        xlast[4] = 0
        xlast[5] = 0
        xlast[8] = 0
        plast[4,:] .= 0
        plast[:,4] .= 0
        plast[5,:] .= 0
        plast[:,5] .= 0
        plast[8,:] .= 0
        plast[:,8] .= 0
        par = [β, N,  ι, η, γ, γd, γh, chp, hfp]
        xplast = hcat(xlast, plast)
        prob = ODEProblem(vectorfield, xplast, (0.0, dt), par)
        xpnext = solve(prob, Tsit5(), saveat = dt).u[2]
        xnext = xpnext[:,1]
        pnext = xpnext[:,2:end]
        for j in 1:dstate
            if xnext[j] < 0
                xnext[j] = 0
            end
        end
        xkkmo[:,i] .= xnext
        for j in 1:dstate
            if pnext[j, j] < 0
                pnext[j, :] .= 0
                pnext[:, j] .= 0
            end
        end
        pkkmo[:,:,i] .= pnext
        
        Σ[:,:,i] = hmat(cov.time[i]) * pkkmo[:,:,i] * hmat(cov.time[i])' + r

        for j in 1:dobs
            if zmiss[i,j]
                zz[j] = 0
            else 
                zz[j] = zloc[i,j]
            end
        end       
        ytkkmo[:,i] = zz - hmat(cov.time[i]) * reshape(xkkmo[:,i], dstate, 1)
        for j in 1:dobs
            if zmiss[i, j]
                zscore = 0
                rdiagadj[j,i] += 0
            else
                sd = sqrt(Σ[j,j,i])
                zscore = ytkkmo[j,i] / sd 
                if abs(zscore) > maxzscore
                    adjzscore = maxzscore / (1 + abs(zscore) - maxzscore)
                    newsd = abs(ytkkmo[j,i]) / adjzscore
                    rdiagadj[j,i] += newsd ^ 2 - sd ^ 2
                else 
                    rdiagadj[j,i] += 0
                end
            end
            Σ[j,j,i] += rdiagadj[j,i]
        end
        k[:,:,i] = pkkmo[:,:,i] * hmat(cov.time[i])' / Σ[:,:,i]
        for j in 1:dobs
            if zmiss[i,j]
                k[:, j ,i] .= 0
            end
        end
        xkk[:,i] = reshape(xkkmo[:,i], dstate) + reshape(k[:,:,i], dstate, dobs) * ytkkmo[:,i]
        pkk[:,:,i] = (I - reshape(k[:,:,i], dstate, dobs) * hmat(cov.time[i])) * pkkmo[:,:,i]
    end
    
    stepdensity = Normal(0, betasd)
    rwlik = 0
    for i in 1:(nobs - 1)
        step = (logβ[i + 1] - log(γ)) - a * (logβ[i] - log(γ))
        rwlik += logpdf(stepdensity, step)
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

function grad(pvar, x, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, a::Float64 = 1., betasd::Float64 = 1., just_nll::Bool = true, γd::Float64 = 365.25 / 1, γh::Float64 = 365.25 / 1, h0::Float64 = 10., τh::Float64 = 10., τd::Float64 = 10., chp::Float64 = 0.01, hfp::Float64 = 0.01)
    g = ForwardDiff.gradient(par -> obj(par, x, z; N = N, η = η, γ = γ, a = a, betasd = betasd, dt = dt, ι = ι, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τd, chp = chp, hfp = hfp), pvar)
    g
end

function hess(pvar, x, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, a::Float64 = 1., betasd::Float64 = 1., just_nll::Bool = true, γd::Float64 = 365.25 / 1, γh::Float64 = 365.25 / 1, h0::Float64 = 10., τh::Float64 = 10., τd::Float64 = 10., chp::Float64 = 0.01, hfp::Float64 = 0.01)
    h = ForwardDiff.hessian(par -> obj(par, x, z; N = N, η = η, γ = γ, a = a, betasd = betasd, dt = dt, ι = ι, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τd, chp = chp, hfp = hfp), pvar)
    h
end

end