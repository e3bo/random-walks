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


function genvectorfield

  rn = @reaction_network begin
     β, S + Y --> L + Y
     η, L --> Y
     rhot * (1 - chp) * γ, Y --> R
     γr, R --> C
     chp * γ, Y --> H
     chp * γ, Y --> Hnew
     Hnew * γhnew, Hnew --> C
     γr * hfp, H --> D
     γr * (1 - hfp), H --> 0
     γd, D --> Drep
  end β η γ γr γhnew γh chp hfp rhot



end



function vectorfield(du, u, par, t)
  x, l, y, removed, hnew, crep, h, d, drep = u[:,1]

  β, N,  ι, η, γ, γd, γh,  γr, γhnew, chp, hfp, rhot = par

  du[1,1] = dx = -β*x*y/N - ι*x
  du[2,1] = dl = β*x*y/N + ι*x - η*l
  du[3,1] = dy = η*l - γ*y
  du[4,1] = dremoved = rhot * (1 - chp) * γ * y - removed * γr
  du[5,1] = dhnew = chp*γ*y - hnew * γhnew
  du[6,1] = dcrep = hnew * γhnew + removed * γr
  du[7,1] = dh = chp*γ*y - γh*h
  du[8,1] = dd = hfp*γh*h - γd*d
  du[9,1] = ddrep = γd*d

  f = [0,                             # x; 
       β*x/N*y/N + ι*x/N,             # x-> l; 1, 2
       η*l/N,                         # l -> y;, 2, 3
       rhot * (1 - chp)*γ*y/N,        # y -> removed; 3, 4
       chp*γ*y/N,                     # y -> hnew, y -> h; 3, 5  (3 , 7)
       (1 - hfp) * γh * h / N,        # h -> not death; 7
       hfp * γh * h / N,              # h -> d; 7, 8
       γd * d / N,                    # d -> drep; 8, 9
       γhnew * hnew,                   # hnew -> crep; 5, 6
       γr * removed]                  # removed -> crep; 4, 6

  q = [  f[1]+f[2]     -f[2]               0              0          0              0                0             0       0
             -f[2] f[2]+f[3]           -f[3]              0          0              0                0             0       0
                 0     -f[3]  f[3]+f[4]+f[5]          -f[4]      -f[5]              0            -f[5]             0       0
                 0         0           -f[4]   f[4] + f[10]          0         -f[10]                0             0       0
                 0         0           -f[5]              0  f[5]+f[9]          -f[9]                0             0       0
                 0         0              0          -f[10]      -f[9]   f[10] + f[9]                0             0       0
                 0         0          -f[5]               0          0              0   f[7]+f[5]+f[6]         -f[7]       0
                 0         0              0               0          0              0            -f[7]   f[7] + f[8]   -f[8]
                 0         0              0               0          0              0                0         -f[8]    f[8]]

  jac= [-β*y/N - ι    0               -β*x/N     0       0       0       0    0 0
         β*y/N + ι   -η                β*x/N     0       0       0       0    0 0
                 0    η                   -γ     0       0       0       0    0 0
                 0    0    rhot *(1 - chp)*γ   -γr       0       0       0    0 0
                 0    0                chp*γ     0  -γhnew       0       0    0 0
                 0    0                    0    γr      γh       0       0    0 0 
                 0    0                chp*γ     0       0       0     -γh    0 0
                 0    0                    0     0       0       0  γh*hfp  -γd 0
                 0    0                    0     0       0       0       0   γd 0]
   
   p = u[:,2:end]
   du[:,2:end] .= jac * p + p * jac' + q * N
end

function obj(pvar::Vector, cov, z; γ::Float64 = 365.25 / 9, dt::Float64 = 0.00273224, ι::Float64 = 0., η::Float64 = 365.25 / 4, N::Float64 = 7e6, a::Float64 = 1., betasd::Float64 = 1., tcsd::Float64 = 0.1, just_nll::Bool = true, maxlogRt::Float64 = 1.6, γd::Float64 = 365.25 / 1, γh::Float64 = 365.25 / 1, γhnew::Float64 = 365.25 / 1, γr::Float64 = 365.25 / 1, h0::Float64 = 10., τh::Float64 = 10., τd::Float64 = 10., chp::Float64 = 0.01, hfp::Float64 = 0.01)

    zloc = deepcopy(z)
    ntauc = cov.τcvecmap[end]
    if size(z, 2) == 3
        l0 = exp(pvar[1])
        h0 = exp(pvar[2])
        τh = exp(pvar[3])
        τd = exp(pvar[4])
        chp = 1 / (1 + exp(-pvar[5]))
        prophomeeffect = pvar[6]
        hfpvec = exp(pvar[7])
        gammad12 = exp(pvar[8])
        gammad34 = exp(pvar[9])
        τcvec = [exp(p) for p in pvar[10:(10 + ntauc - 1)]]
        bvec = pvar[(10 + ntauc):end]
    elseif size(z, 2) == 2 && !("hospitalizations" in names(z))
        l0 = exp(pvar[1])
        h0 = exp(pvar[2])
        τd = exp(pvar[3])
        prophomeeffect = pvar[4]
        hfpvec = exp(pvar[5])
        gammad12 = exp(pvar[6])
        gammad34 = exp(pvar[7])
        τcvec = [exp(p) for p in pvar[8:(8 + ntauc - 1)]]
        bvec = pvar[(8 + ntauc):end]
        zloc[!, "hospitalizations"] .= missing
    end
    zloc = Matrix(select(zloc, :cases, :hospitalizations, :deaths)) # ensure assumed column order
    
    d0 = h0 * γh / γd * hfpvec[1]
    
    y0 = l0 * η / γ
    # x, l, y, removed, hnew, crep, h, d, drep
    
    x0 = [             max(N - l0 - y0 - h0 - d0, 100) 
                                            min(l0, N) 
                                            min(y0, N)  
          min(cov.rhot[1] * y0 * γ * (1 - chp)/ γr, N)
                          min(chp * y0 * γ / γhnew, N)
                                                     0 
                                            min(h0, N)
                                            min(d0, N)    
                                                     0]
                                                     
    p0 = convert(Array{eltype(bvec), 2}, Diagonal([1, 1, 1, 1, 1, 0, 1, 1, 0]))
    
    if (eltype(bvec) == Float64)
        println(pvar)
    end
    
    dstate = size(x0, 1)
    dobs = size(zloc, 2)

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
    
    @assert length(bvec) == cov.bvecmap[end] "length of bvec should equal the final value in cov.bvecmap"
    
    hmat =  [0 0 0 0 0 1 0 0 0; 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0 1] 
    
    for i in 1:nobs
        if i == 1
            xlast = x0
            plast = p0
        else
            xlast = xkk[:,i - 1]
            plast = pkk[:,:,i-1]
        end
        β = exp(bvec[cov.bvecmap[i]] + cov.prophomeiqr[i] * prophomeeffect)
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
    for i in 1:(length(bvec) - 1)
        step = bvec[i + 1] - bvec[i]
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