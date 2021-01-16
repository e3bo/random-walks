

z = [1393, 1360, 1592, 1143]
times = [2020.768, 2020.770, 2020.773, 2020.776]

ytkk = similar(times)
ytk = similar(times)
s = similar(times)
n = length(z)

using LinearAlgebra

function obj(a)
    # prior for time 0
    x0 = [-1., 1.]
    p0 = Matrix(1.0I, 2, 2)

    # dynamics
    Φ = [0.8 0.2; -0.1 a]
    q = [0.2 0.0; 0.0 0.5]

    # observation
    h = [1.0 0.0]
    r = Matrix(0.3I, 1, 1)

    # (mock) data
    z = [[-1.77], [-0.78], [-1.28], [-1.06], [-3.65], [-2.47], [-0.06], [-0.91], [-0.80], [1.48]]

    # filter (assuming first observation at time 1)
    n = length(z)

    s = zeros(Float64, 1, 1, n)
    ytkkmo = zeros(Float64, 1, n)
    k = zeros(Float64, 2, n)
    xkk = zeros(Float64, 2, n)
    xkkmo = zeros(Float64, 2, n)
    pkk = zeros(Float64, 2, 2, n)
    pkkmo = zeros(Float64, 2, 2, n)

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
    (nll, ytkkmo, xkkmo)
end