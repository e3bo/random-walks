using LinearAlgebra
using Optim

function obj(x::AbstractArray{T}) where {T}
    # prior for time 0
    x0 = [-1., 1.]
    p0 = Matrix(1.0I, 2, 2)

    # dynamics
    Φ = [0.8 0.2; -0.1 x[1]]
    q = [0.2 0.0; 0.0 0.5]

    # observation
    h = [1.0 0.0]
    r = Matrix(0.3I, 1, 1)

    
    # (mock) data
    z = [[-1.77], [-0.78], [-1.28], [-1.06], [-3.65], [-2.47], [-0.06], [-0.91], [-0.80], [1.48]]

    # filter (assuming first observation at time 1)
    n = length(z)

    s = Array{T}(undef, 1, 1, n)
    ytkkmo = Array{T}(undef, 1, n)
    k = Array{T}(undef, 2, n)
    xkk = Array{T}(undef, 2, n)
    xkkmo = Array{T}(undef, 2, n)
    pkk = Array{T}(undef, 2, 2, n)
    pkkmo = Array{T}(undef, 2, 2, n)
    
    for i in 1:n
        if (i == 1)
            xkkmo[:,i] = Φ*x0
            pkkmo[:,:,i] = Φ*p0*Φ' + q
        else
            xkkmo[:,i] = Φ*xkk[:,i - 1]
            pkkmo[:,:,i] = Φ*pkk[:,:,i-1]*Φ' + q
        end
        s[:,:,i] = h * pkkmo[:,:,i] * h' + r
        k[:,i] = pkkmo[:,:,i] * h' / s[:,:,i]
        ytkkmo[:,i] = z[i] + h * reshape(xkkmo[:,i], 2, 1)
        xkk[:,i] = reshape(xkkmo[:,i], 2, 1) + reshape(k[:,i], 2, 1) * ytkkmo[:,i]
        pkk[:,:,i] = (I - reshape(k[:,i], 2, 1) * h) * pkkmo[:,:,i]
    end        
    nll = 0.5 * (sum(ytkkmo[1,:] .^2 ./ s[1,1,:] + map(log, s[1,1,:])) + n * log(2 * pi))
    nll
end

ans0 = optimize(obj, [-10], [10], [0.3], Fminbox(BFGS())) # works, but no autodiff
ans1 = optimize(obj, [-1.0], LBFGS()) # works if init is close to minimizer
ans2 = optimize(obj, [-1.0], LBFGS(); autodiff = :forward) # reduced iterations
ans3 = optimize(obj, [-10], [10], [8.], Fminbox(LBFGS()); autodiff = :forward) # can start far away 