# Crews_Econ35310_FindA

function Crews_FindA(L::Float64,N::Array{Int64,2},w::Float64,tau::Float64,
    p_F::Array{Float64,2},P_F::Float64,E_F::Float64,sigma::Float64,rho::Float64,
    phi::Array{Float64,1},alpha_H::Array{Float64,1},alpha_F::Array{Float64,1},
    beta_F::Array{Float64,1},fixed::Array{Float64,2})
    #=
    The function takes fourteen arguments:
    + labor endowment L
    + fixed network N,
    + proposed wage w,
    + iceberg costs tau,
    + foreign price vector p_F,
    + foreign price index P_F,
    + foreign expenditure E_F,
    + consumption elasticity of substitution sigma,
    + production elasticity of substitution rho,
    + vector: productivities phi,
    + vector: domestic production salience alpha_H,
    + vector: foreign production salience alpha_F,
    + vector: foreign demand salience beta_F,
    + vector: fixed costs =fixed=

    Given these inputs, the function solves for
    + aggregate demand A,
    + total cost =total_cost=,
    + labor share share_L,
    + total fixed costs (without =w= weight) total_fixed_noW,
    + total foreign input share totFshare
    such that equations (15) and Appx G.1(2b) from TKMD (2018) hold.
    =#

    A_old = 1.1
    A_new = 1.0

    while abs(A_new - A_old) > 10e-7
        A_old = A_new
        P, total_profit, total_cost, share_L, total_fixed_noW, totFshare = Crews_Q2Solver(N, A_old, w, tau, p_F, P_F, E_F, sigma, rho, phi, alpha_H, alpha_F, beta_F, fixed)
        A_new = (w*L + sum(total_profit)) * P^(sigma-1)   # TKMD eq. 15 & G.1 step 2b
        A_new = 0.05*A_new + 0.95*A_old
        if abs(A_new - A_old) > 10e50
            println("Aggregate demand blew up.")
            break
        end
    end

    return A_new, total_cost, share_L, total_fixed_noW, totFshare
end
