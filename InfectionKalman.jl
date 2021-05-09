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
export genvectorfield

function genvectorfield()

  rn = @reaction_network begin
     β / N, X + Y --> L + Y
     η, L --> Y
     ρ * (1 - p_h) * γ, Y --> Z
     (1 - ρ) * (1 - p_h) * γ, Y --> 0
     γ_z, Z --> Z_r
     p_h * γ, Y --> H + A + Z
     γ_h * p_d, H --> D
     γ_h * (1 - p_d), H --> 0
     γ_d, D --> D_r
  end β N η γ γ_d γ_z γ_h p_h p_d ρ

  odesys = convert(ODESystem, rn)
  odefun = ODEFunction(odesys, jac=true)
  
  s = netstoichmat(rn)

  #u0 = [1e7, 10, 1, 10, 0, 0, 0, 0, 0] 
  #pars = [400, 365/4, 91.24, 365 / 10, 365 / 1, 365 / 10, .1, .1, .4]
  
  jls = [jumpratelaw(r) for r in reactions(rn)]
  ffuns =  [build_function(jl, states(rn), parameters(rn), independent_variable(rn), expression=Val{false}) for jl in jls] # there may be a nicer way to do this soon: https://github.com/SciML/Catalyst.jl/issues/306
  
  function vf(du, u, par, t)
    du[:,1] .= odefun(u[:,1], par, t)
    f = [ffun(u, par, t) for ffun in ffuns]
    q = s' * Diagonal(f) * s
    jac = odefun.jac(u[:,1], par, t)
    p = u[:,2:end]
    du[:,2:end] .= jac * p + p * jac' + q
  end
end

function obj(w::Vector, cov, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, η::Float64 = 365.25 / 4, N::Float64 = 7e6, β_0sd::Float64 = 1., τ_csd::Float64 = 0.1, p_hsd::Float64 = 0.5, just_nll::Bool = true, γ_d::Float64 = 365.25 / 1, γ_h::Float64 = 365.25 / 1, γ_z::Float64 = 365.25 / 1, H0::Float64 = 10., τ_h::Float64 = 10., τ_d::Float64 = 10., p_h::Float64 = 0.01, p_d::Float64 = 0.01)

    zloc = deepcopy(z)
    
    nτ_c = cov.τ_cmap[end]
    np_h = cov.p_hmap[end]
    
    if size(z, 2) == 3
        τ_h = exp(w[3])
        if "doses_scaled" in names(cov)
            doseeffect = w[6]
            w2 = vcat(w[1:2], w[4:5], w[7:end])
        else 
            w2 = vcat(w[1:2], w[4:end])
        end
        p_h = [1 / (1 + exp(-p)) for p in w2[(9 + nτ_c):(9 + nτ_c + np_h - 1)]]
    elseif size(z, 2) == 2 && !("hospitalizations" in names(z))
        np_h = 0
        w2 = w
        zloc[!, "hospitalizations"] .= missing
    end
    
    L0 = exp(w2[1])
    H0 = exp(w2[2])
    τ_d = exp(w2[3])
    residentialeffect = w2[4]
    p_d = 1 / (1 + exp(-w2[5]))
    γ_d12 = exp(w2[6])
    γ_d34 = exp(w2[7])
    γ_z17 = exp(w2[8])
    
    τ_c = [exp(p) for p in w2[9:(9 + nτ_c - 1)]]
    β_0 = w2[(9 + nτ_c + np_h):end]

    zloc = Matrix(select(zloc, :cases, :hospitalizations, :deaths)) # ensure assumed column order
    
    D0 = H0 * γ_h / γ_d * p_d
    Y0 = L0 * η / γ
    x0 = [         max(N - L0 - Y0 - H0 - D0, N * 0.1) #X
                                            min(Y0, N) #Y
                                            min(L0, N) #L
           min(cov.ρ[1] * Y0 * γ * (1 - p_h[cov.p_hmap[1]]) / γ_z, N) #Z
                                                    0  #Z_r
                                            min(H0, N) #H
                                                    0  #A
                                            min(D0, N) #D
                                                    0] #D_r
                                                     
    p0 = convert(Array{eltype(β_0), 2}, Diagonal([1, 1, 1, 1, 0, 1, 0, 1, 0]))
    
    if (eltype(β_0) == Float64)
        println(w)
    end
    
    dstate = size(x0, 1)
    dobs = size(zloc, 2)

    zmiss = [ismissing(x) for x in zloc]
    zz = Array{eltype(zloc)}(undef, dobs)

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
    
    hmat =  [0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 0 1] 
    vectorfield = genvectorfield()
    zerovars = [5, 7, 9]

    for i in 1:nobs
        if i == 1
            xlast = x0
            plast = p0
        else
            xlast = xkk[:,i - 1]
            plast = pkk[:,:,i-1]
        end
        
        if "doses_scaled" in names(cov)
            β = exp(β_0[cov.β_0map[i]] + cov.residential[i] * residentialeffect + cov.dose_scaled[i] * doseeffect)
        else
            β = exp(β_0[cov.β_0map[i]] + cov.residential[i] * residentialeffect)
        end
        if β > 1000
            β = 1000
        end
        for zv in zerovars
            xlast[zv] = 0
            plast[zv,:] .= 0
            plast[:,zv] .= 0
        end
        par = [β, N, η, γ, γ_d, γ_z, γ_h, p_h[cov.p_hmap[i]], p_d, cov.ρ[i]]
        if cov.wday[i] == 1
           par[5] = γ_d12
           par[6] = γ_z17
        elseif cov.wday[i] == 2
           par[5] = γ_d12
        elseif cov.wday[i] == 3 || cov.wday[i] == 4
           par[5] = γ_d34
        elseif cov.wday[i] == 7
           par[6] = γ_z17
        end
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
        r = Diagonal([τ_c[cov.τ_cmap[i]] * xlast[3], τ_h, τ_d]) 
        
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
    
    stepdensity = Normal(0, β_0sd)
    rwlik = 0
    for i in 1:(length(β_0) - 1)
        step = β_0[i + 1] - β_0[i]
        rwlik += logpdf(stepdensity, step)
    end
    
    τ_cstepdensity = Normal(0, τ_csd)
    for i in 1:(length(τ_c) - 1)
        τ_cstep =  log(τ_c[i + 1]) - log(τ_c[i])
        rwlik += logpdf(τ_cstepdensity, τ_cstep)
    end
    
    p_hstepdensity = Normal(0, p_hsd)
    for i in 1:(length(p_h) - 1)
        p_hstep =  log(p_h[i + 1]) - log(1 - p_h[i + 1]) - (log(p_h[i]) - log(1 - p_h[i]))
        rwlik += logpdf(p_hstepdensity, p_hstep)
    end
    
    kflik = 0
    for i in 1:nobs
         sel = [!x for x in zmiss[i,:]]
         kflik -= 0.5 * (ytkkmo[sel,i]' / Σ[sel, sel, i] * ytkkmo[sel,i] + log(det(Σ[sel, sel, i]))  + dobs * log(2 * pi))
    end
    nll = -kflik - rwlik
    
    if just_nll
       return nll
    else
       return nll, ytkkmo, Σ, xkkmo, pkkmo, pkk, rdiagadj
    end
end

function grad(w, x, z; args...)
    g = ForwardDiff.gradient(par -> obj(par, x, z; args...), w)
    g
end

function hess(w, x, z; args...)
    h = ForwardDiff.hessian(par -> obj(par, x, z; args...), w)
    h
end

end