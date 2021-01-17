using LinearAlgebra
using Optim

function obj(x::Vector; γ::Float64 = 365 / 9)
    # prior for time 0
    x0 = [19e6; 1e4; 1e4; 0]
    p0 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]

    dstate = size(x0, 1)
    # dynamics
    Φ = Matrix(x[1] * I, dstate, dstate)
    q = Matrix(γ*I, dstate, dstate)

    # observation
    h = [0 0 0 1]
    dobs = size(h, 1)
    r = Matrix(0.3I, dobs, dobs)

    # data
    z = [[1393], [1360], [1836], [1592], [1447], [1143]]

    # filter (assuming first observation at time 1)
    nobs = length(z)

    s = Array{eltype(x)}(undef, dobs, dobs, nobs)
    ytkkmo = Array{eltype(x)}(undef, dobs, nobs)
    k = Array{eltype(x)}(undef, dstate, nobs)
    xkk = Array{eltype(x)}(undef, dstate, nobs)
    xkkmo = Array{eltype(x)}(undef, dstate, nobs)
    pkk = Array{eltype(x)}(undef, dstate, dstate, nobs)
    pkkmo = Array{eltype(x)}(undef, dstate, dstate, nobs)
    
    for i in 1:nobs
        if (i == 1)
            xkkmo[:,i] = Φ*x0
            pkkmo[:,:,i] = Φ*p0*Φ' + q
        else
            xkkmo[:,i] = Φ*xkk[:,i - 1]
            pkkmo[:,:,i] = Φ*pkk[:,:,i-1]*Φ' + q
        end
        s[:,:,i] = h * pkkmo[:,:,i] * h' + r
        k[:,i] = pkkmo[:,:,i] * h' / s[:,:,i]
        ytkkmo[:,i] = z[i] + h * reshape(xkkmo[:,i], dstate, 1)
        xkk[:,i] = reshape(xkkmo[:,i], dstate, 1) + reshape(k[:,i], dstate, 1) * ytkkmo[:,i]
        pkk[:,:,i] = (I - reshape(k[:,i], dstate, 1) * h) * pkkmo[:,:,i]
    end        
    nll = 0.5 * (sum(ytkkmo[1,:] .^2 ./ s[1,1,:] + map(log, s[1,1,:])) + nobs * log(2 * pi))
    nll
end

ans0 = optimize(obj, [-10], [10], [0.3], Fminbox(BFGS())) # works, but no autodiff
ans1 = optimize(obj, [-1.0], LBFGS()) # works if init is close to minimizer
ans2 = optimize(obj, [-1.0], LBFGS(); autodiff = :forward) # reduced iterations
ans3 = optimize(obj, [-10], [10], [8.], Fminbox(LBFGS()); autodiff = :forward) # can start far away 