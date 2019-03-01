# Crews_Econ35310_FindW

function Crews_FindW(L::Float64,N::Array{Int64,2},tau::Float64,
    p_F::Array{Float64,2},P_F::Float64,E_F::Float64,sigma::Float64,rho::Float64,
    phi::Array{Float64,1},alpha_H::Array{Float64,1},alpha_F::Array{Float64,1},
    beta_F::Array{Float64,1},fixed::Array{Float64,2})
    #=
    The function takes eight arguments:
    + labor endowment L
    + fixed network N,
    + proposed aggregate demand A,
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
    + wage =w=,
    + total foreign input share totFshare
    such that equations (26) of TKMD (2018) holds.
    =#

    w_old = 1.1
    w_new = 1.0

    while abs(w_new - w_old) > 10e-7
        w_old = w_new
        A, total_cost, share_L, total_fixed_noW, totFshare = Crews_FindA(L, N, w_old, tau, p_F, P_F, E_F, sigma, rho, phi, alpha_H, alpha_F, beta_F, fixed)
        w_new = (share_L' * total_cost) / (L - sum(total_fixed_noW))   # TKMD eq. 26
        w_new = 0.05*w_new + 0.95*w_old
        if abs(w_new - w_old) > 10e15
            println("The nominal wage blew up.")
            break
        end
    end

    return w, totFshare
end
