

z = [1393, 1360, 1592, 1143]
times = [2020.768, 2020.770, 2020.773, 2020.776]

ytkk = similar(times)
ytk = similar(times)
s = similar(times)
n = length(z)


## hacked example from Kalman.jl, by Moritz Schauer (GitHub user mschauer)
using Kalman, GaussianDistributions, LinearAlgebra
using GaussianDistributions: ⊕ # independent sum of Gaussian r.v.
using Statistics

# prior for time 0
x0 = [-1., 1.]
P0 = Matrix(1.0I, 2, 2)

# dynamics
Φ = [0.8 0.2; -0.1 0.8]
b = zeros(2)
Q = [0.2 0.0; 0.0 0.5]

# observation
H = [1.0 0.0]
R = Matrix(0.3I, 1, 1)

# (mock) data
ys = [[-1.77], [-0.78], [-1.28], [-1.06], [-3.65], [-2.47], [-0.06], [-0.91], [-0.80], [1.48]]


# filter (assuming first observation at time 1)
N = length(ys)

s = zeros(Float64, 1, 1, N)
yres = zeros(Float64, 1, N)
p = Gaussian(x0, P0)
ps = [p] # vector of filtered Gaussians
for i in 1:N
    global p
    # predict
    p = Φ*p ⊕ Gaussian(zero(x0), Q) #same as Gaussian(Φ*p.μ, Φ*p.Σ*Φ' + Q)
    # correct
    p, yres[:,i], s[:,:,i] = Kalman.correct(Kalman.JosephForm(), p, (Gaussian(ys[i], R), H))
    push!(ps, p) # save filtered density
end

nll = 0.5 * (sum(yres[1,:] .^2 ./ s[1,1,:] + map(log, s[1,1,:])) + N * log(2 * pi))
