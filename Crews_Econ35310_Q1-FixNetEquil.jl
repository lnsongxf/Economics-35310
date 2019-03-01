## Crews_Econ35310_FixNetEquil

function Crews_FixNetEquil(E_EF::Array{Float64,2},N::Array{Float64,2},
    S::SparseMatrixCSC{Float64,Int64},
    P::Float64,p_H::Array{Float64,2},P_F::Float64,p_F::Array{Float64,2},
    w::Float64,tau::Float64,sigma::Float64,rho::Float64)
    #=
    The function takes eleven arguments:
    + an initial guess for aggregate expenditures E_EF = [E, E_F]
    + a fixed network N,
    + a sparse matrix of shares S,
    + foreign price index P,
    + foreign price vector p_H,
    + domestic price index P_F,
    + domestic price vector p_F
    + a proposed wage w,
    + iceberg costs tau,
    + consumption elasticity of substitution sigma,
    + production elasticity of substitution rho,

    Given N, P, p_H, P_F, p_F, w, tau, sigma, rho, and initial guesses E and E_F,
    the function solves for
    + the labor endowment L,
    + aggregate expenditure E,
    + aggregate foreign expenditure E_F
    such that equations (13)-(17) from TKMD (2018) hold.
    =#

    # Retrieve number of firms in simulation
    num = size(N)[1]
    cards = N[:,2]
    maxcard = Int64(maximum(cards))

    # Retrieve firm characteristic vectors
    IndI = N[:,maxcard+3]
    IndE = N[:,maxcard+4]
    phi = N[:,maxcard+5]
    X = N[:,maxcard+6]
    Y = N[:,maxcard+7]
    alpha_F = N[:,maxcard+8]
    beta_F = N[:,maxcard+9]
    mu = sigma / (sigma-1)

    # Demand and sales
    E = E_EF[1,1]
    E_F = E_EF[1,2]
    q_F = beta_F.^(sigma-1) .* ((tau*p_H).^(-sigma) / (P_F^(1-sigma))) .* IndE
    x_F = tau * p_H .* q_F * E_F # Foreign sales to final consumers


    q = p_H.^(-sigma)/(P^(1-sigma))
    x_H = p_H .* q * E

    domestic_sale_sum = gmres(sparse(I,num,num) - S[:,1:num], S[:,1:num]*(x_H + x_F)/mu)
    total_cost = (x_H + x_F)/mu + domestic_sale_sum
    x_D = total_cost' .* S[:,1:num]
    total_sale = x_H + x_F + domestic_sale_sum          # Nx1 vector
    total_profit = total_sale - total_cost              # Nx1 vector

    #=
    # Total sales and profits
    DomHH = zeros(num,1)
    x = zeros(num,1)
    pi = zeros(num,1)
    for j=1:num
        Z = N[j,3:maxcard+2]
        filter!(e -> e>0, Z)
        Z = Array{Int64}(Z)
        domHH = (phi[j,1]*P / mu)^(sigma-1) * (Theta[j,1])^((sigma-1)/(rho-1)) * E
        profHH = domHH / sigma
        exports = IndE[j,1] * (phi[j,1]*beta_F[j,1]*P / (tau*mu))^(sigma-1) * (Theta[j,1])^((sigma-1)/(rho-1)) * E_F
        profE = exports / sigma
        domFirms = 0
        for i = Z
            domFirms += (phi[j,1]*X[j,1]*Y[i,1])^(rho-1) * (x[i,1]*c[i,1]) * (Theta[j,1] / Theta[i,1])
        end
        DomHH[j,1] = domHH
        x[j,1] = domHH + exports + domFirms      # TKMD eq (13)
        pi[j,1] = profHH + profE      # TKMD eq (14)
    end

    total_cost = x[:,1] .- pi[:,1]
    =#

    share_L = Array{Float64}(S[:, num+1])
    share_F = Array{Float64}(S[:, num+2])

    L = share_L' * total_cost[:,1] / w     # TKMD eq (17) labor market clearing

    TB = 10e10*(sum(x_F[:,1]) - share_F' * total_cost[:,1])    # TKMD eq (16) trade balance

    MC = 10e10*(E - w*L - sum(total_profit[:,1]) + TB)
    #MC = 10e10*(E - w*L - sum(pi[:,1]) + TB)     # TKMD eq (15) goods market clearing

    MC_TB = [MC, TB]

    return MC_TB, L, total_cost, total_sale, total_profit, q, x_H, x_D, x_F
    #return MC_TB, L, pi, x, total_cost
end
