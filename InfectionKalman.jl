module InfectionKalman

using Catalyst
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

function genvectorfield()

  rn = @reaction_network begin
     β, S + Y --> L + Y
     η, L --> Y
     ρ * (1 - p_h) * γ, Y --> R
     (1 - ρ) * (1 - p_h) * γ, Y --> 0
     γ_r, R --> C
     p_h * γ, Y --> H + A
     γ_a, A --> C
     γ_r * p_d, H --> D
     γ_r * (1 - p_d), H --> 0
     γ_d, D --> D_r
  end β η γ γ_d γ_r γ_a γ_h p_h p_d ρ

  odesys = convert(ODESystem, rn)
  odefun = ODEFunction(odesys, jac=true)
  
  F = (u, p, t) -> odefun(u, p, t)
  J = (u, p, t) -> odefun.jac(u,p, t)
  S = netstoichmat(rn)

  #u0 = [1e7, 10, 1, 10, 0, 0, 0, 0, 0] 
  #pars = [400, 365/4, 91.24, 365 / 10, 365 / 1, 365 / 1, 365 / 10, .1, .1, .4]
  
  jls = [jumpratelaw(r) for r in reactions(rn)]
  ffuns =  [build_function(jl, states(rn), parameters(rn), independent_variable(rn), expression=Val{false}) for jl in jls] # there may be a nicer way to do this soon: https://github.com/SciML/Catalyst.jl/issues/306
  
  function vf(du, u, par, t)
    du[:,1] .= F(u[:,1], par, t)
    f = [ffun(u, par, t) for ffun in ffuns]
    q = S' * Diagonal(f) * S
    jac = J(u[:,1], par, t)
    du[:,2:end] .= jac * p + p * jac' + q
  end
end


function obj(pvar::Vector, cov, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, η::Float64 = 365.25 / 4, N::Float64 = 7e6, β_0sd::Float64 = 1., τ_csd::Float64 = 0.1, just_nll::Bool = true, γ_d::Float64 = 365.25 / 1, γ_h::Float64 = 365.25 / 1, γ_a::Float64 = 365.25 / 1, γ_r::Float64 = 365.25 / 1, H0::Float64 = 10., τ_h::Float64 = 10., τ_d::Float64 = 10., p_h::Float64 = 0.01, p_d::Float64 = 0.01)

    zloc = deepcopy(z)
    nτ_c = cov.τ_cmap[end]
    if size(z, 2) == 3
        L0 = exp(pvar[1])
        H0 = exp(pvar[2])
        τ_h = exp(pvar[3])
        τ_d = exp(pvar[4])
        p_h = 1 / (1 + exp(-pvar[5]))
        prophomeeffect = pvar[6]
        p_d = exp(pvar[7])
        γ_d12 = exp(pvar[8])
        γ_d34 = exp(pvar[9])
        τ_c = [exp(p) for p in pvar[10:(10 + nτ_c - 1)]]
        β_0 = pvar[(10 + nτ_c):end]
    elseif size(z, 2) == 2 && !("hospitalizations" in names(z))
        L0 = exp(pvar[1])
        H0 = exp(pvar[2])
        τ_d = exp(pvar[3])
        prophomeeffect = pvar[4]
        p_d = exp(pvar[5])
        γ_d12 = exp(pvar[6])
        γ_d34 = exp(pvar[7])
        τ_c = [exp(p) for p in pvar[8:(8 + nτ_c - 1)]]
        β_0 = pvar[(8 + nτ_c):end]
        zloc[!, "hospitalizations"] .= missing
    end
    zloc = Matrix(select(zloc, :cases, :hospitalizations, :deaths)) # ensure assumed column order
    
    D0 = H0 * γ_h / γ_d * p_d
    Y0 = L0 * η / γ
    
    x0 = [         max(N - L0 - Y0 - H0 - D0, N * 0.1) #S
                                            min(Y0, N) #Y
                                            min(L0, N) #L
           min(cov.ρ[1] * Y0 * γ * (1 - p_h) / γ_r, N) #R
                                                    0  #C
                                            min(H0, N) #H
                            min(Y0 * γ * p_h / γ_a, N) #A
                                            min(D0, N) #D
                                                    0] #D_r
                                                     
    p0 = convert(Array{eltype(β_0), 2}, Diagonal([1, 1, 1, 1, 0, 1, 1, 1, 0]))
    
    if (eltype(β_0) == Float64)
        println(pvar)
    end
    
    dstate = size(x0, 1)
    dobs = size(zloc, 2)

    zmiss = [ismissing(x) for x in zloc]
    zz = Array{eltype(zloc)}(undef, dobs)
    maxzscore = Inf

    nobs = size(zloc, 1)
    rdiagadj = ones(eltype(β_0), dobs, nobs) # make non-zero to ensure Σ is not singular
    Σ = Array{eltype(β_0)}(undef, dobs, dobs, nobs)
    ytkkmo = Array{eltype(β_0)}(undef, dobs, nobs)
    k = Array{eltype(β_0)}(undef, dstate, dobs, nobs)
    xkk = Array{eltype(β_0)}(undef, dstate, nobs)
    xkkmo = Array{eltype(β_0)}(undef, dstate, nobs)
    pkk = Array{eltype(β_0)}(undef, dstate, dstate, nobs)
    pkkmo = Array{eltype(β_0)}(undef, dstate, dstate, nobs)
    
    @assert length(β_0) == cov.β_0map[end] "length of β_0 should be equal to the final value in cov.β_0map"
    
    hmat =  [0 0 0 0 0 1 0 0 0; 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0 1] 
    vectorfield = genvectorfield()

    for i in 1:nobs
        if i == 1
            xlast = x0
            plast = p0
        else
            xlast = xkk[:,i - 1]
            plast = pkk[:,:,i-1]
        end
        β = exp(β_0[cov.β_0map[i]] + cov.prophomeiqr[i] * prophomeeffect)
        xlast[6] = 0
        xlast[9] = 0
        plast[6,:] .= 0
        plast[:,6] .= 0
        plast[9,:] .= 0
        plast[:,9] .= 0
        hfp = hfpvec[1]
      
        par = [β, N,  ι, η, γ, γd, γh,  γr, γhnew, chp, hfp, cov.rhot[i]]
        if cov.wday[i] == 1 || cov.wday[i] == 2
           par[6] = gammad12
        elseif cov.wday[i] == 3 || cov.wday[i] == 4
           par[6] = gammad34
        end
        xplast = hcat(xlast, plast)
        println(i)
        println(xlast)
        println(par)
        println("\n")
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
        r = Diagonal([τcvec[cov.τcvecmap[i]] * xlast[3], τh, τd]) 
        
        Σ[:,:,i] = hmat * pkkmo[:,:,i] * hmat' + r

        for j in 1:dobs
            if zmiss[i,j]
                zz[j] = 0
            else 
                zz[j] = zloc[i,j]
            end
        end       
        ytkkmo[:,i] = zz - hmat * reshape(xkkmo[:,i], dstate, 1)
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
    for i in 1:(length(β_0) - 1)
        step = β_0[i + 1] - β_0[i]
        rwlik += logpdf(stepdensity, step)
    end
    
    tcstepdensity = Normal(0, tcsd)
    for i in 1:(length(τcvec) - 1)
        tcstep =  log(τcvec[i + 1]) - log(τcvec[i])
        rwlik += logpdf(tcstepdensity, tcstep)
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

function grad(pvar, x, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, a::Float64 = 1., betasd::Float64 = 1., tcsd::Float64 = 0.1, just_nll::Bool = true, γd::Float64 = 365.25 / 1, γh::Float64 = 365.25 / 1, h0::Float64 = 10., τh::Float64 = 10., τd::Float64 = 10., chp::Float64 = 0.01, hfp::Float64 = 0.01)
    g = ForwardDiff.gradient(par -> obj(par, x, z; N = N, η = η, γ = γ, a = a, betasd = betasd, tcsd = tcsd, dt = dt, ι = ι, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τd, chp = chp, hfp = hfp), pvar)
    g
end

function hess(pvar, x, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, a::Float64 = 1., betasd::Float64 = 1., tcsd::Float64 = 0.1, just_nll::Bool = true, γd::Float64 = 365.25 / 1, γh::Float64 = 365.25 / 1, h0::Float64 = 10., τh::Float64 = 10., τd::Float64 = 10., chp::Float64 = 0.01, hfp::Float64 = 0.01)
    h = ForwardDiff.hessian(par -> obj(par, x, z; N = N, η = η, γ = γ, a = a, betasd = betasd, tcsd = tcsd, dt = dt, ι = ι, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τd, chp = chp, hfp = hfp), pvar)
    h
end

end