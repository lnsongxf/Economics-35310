# Crews_Econ35310_Q2
## Author: Levi Crews
### I acknowledge help from Younghun Shim

###############################################################################

# Load packages
using Distributions, StatsBase, Random, Statistics
using LinearAlgebra, SparseArrays, NLsolve, IterativeSolvers

# Load user-defined functions
# Change path as needed
include("Crews_Econ35310_FindW.jl")
include("Crews_Econ35310_FindA.jl")
include("Crews_Econ35310_Q2Solver.jl");

# Set seed
Random.seed!(4321)

# Initialize network
num = 100000 # note: I was running out of memory with 100,000 firms, so I did 10,000
maxcard = 300
Nmu = 3.5
Nsigma = 0.8

N = zeros(num, 2)
N[:, 1] = collect(1:num)
for i = 1:maxcard
    temp = min(rand(LogNormal(Nmu, Nsigma)), i-1)
    N[i, 2] = ceil(Int64, temp)
end
temp = min.(rand(LogNormal(Nmu, Nsigma), num-maxcard), maxcard*ones(num-maxcard,1))
N[maxcard+1:end, 2] = ceil.(Int64, temp)
N = Array{Int64}(N);

# Initialize parameters
Phi_ADscale = -4.42
Phi_AFscale = -2.22
Phi_BFscale = -2.01
Phi_ABdisp = 2.34
f = 2.0

# Draw firm characteristics
phi = rand(LogNormal(-1.52, 0.85), num) # \phi_j
X = rand(Beta(0.1, 0.9), num) # X_j
Y = rand(LogNormal(Phi_ADscale, Phi_ABdisp), num) # Y_{kj} for all k
alpha_H = X .* Y
alpha_F = rand(LogNormal(Phi_AFscale, Phi_ABdisp), num) # \alpha_{Fj}
beta_F = rand(LogNormal(Phi_BFscale, Phi_ABdisp), num) # \beta_{jF}
fixed = f*ones(num); # fixed cost

# Calibrate elasticities and labor supply
sigma = 4.0
rho = 2.0
tau = 1.05
mu = sigma / (sigma - 1)
L = 1.0;

# Arbitrarily choose foreign variable price vector
p_F = (rand(num,1) .+ 1.25) # p_F = 1.75 * ones(num,1)
P_F = 0.05 # P_F = (sum(p_F.^(1-sigma)))^(1/(sigma-1))
E_F = 0.7;

###############################################################################

w, totFshare = Crews_FindW(L, N, tau, p_F, P_F, E_F, sigma, rho, phi, alpha_H, alpha_F, beta_F, fixed)
if w >= 0
    println("w = ", w)
    println("Hooray! We found an equilibrium.")
else
    println("w = ", w)
    println("We failed to find an equilibrium.")
end

###################################################################################

#=
# 10% increase in foreign price
p_F_hat = 1.1

# guess w_hat
x_hat = ones(num,1)
x_hat_previous = 1.5*ones(num,1)
w_hat = 1    # guess
w_hat_previous = 0.9

while abs(w_hat_previous-w_hat) > 10e-9
    w_hat_previous = w_hat
    c_hat = ((ones(num,1)-totFshare)*(w_hat^(1-rho)) + totFshare*(p_F_hat^(1-rho))).^(1/(1-rho))
    x_hat_F = c_hat.^(1-sigma) .* IndI
    P_hat = (c_hat'.^(1-sigma)*share_H)^(1/(1-sigma))
    share_L_hat = w_hat^(1-rho) * c_hat.^(rho-1)
    share_hat = c_hat'.^(1-rho).*SparseNet.*c_hat.^(1-rho)
    while norm(x_hat-x_hat_previous) > 10e-30
        x_hat_previous = x_hat
        E_hat = (w*L/E*w_hat) + (total_profit[:,1]'*x_hat_previous[:,1]/E)
        x_hat = x_F./total_sale.*x_hat_F + x_H./total_sale.*c_hat.^(1-sigma)*P_hat^(sigma-1)*E_hat + share_hat.*x_D*x_hat_previous./total_sale
        w_hat = (share_L .* share_L_hat)'*(total_cost .* x_hat/(w*L))
        w_hat = w_hat[1]
    end
end

real_wage_change = w_hat/P_hat
=#
