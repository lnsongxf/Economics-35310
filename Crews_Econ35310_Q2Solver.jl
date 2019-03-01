# Crews_Econ35310_Q2Solver

function Crews_Q2Solver(N::Array{Int64,2},A::Float64,w::Float64,tau::Float64,
    p_F::Array{Float64,2},P_F::Float64,E_F::Float64,sigma::Float64,rho::Float64,
    phi::Array{Float64,1},alpha_H::Array{Float64,1},alpha_F::Array{Float64,1},
    beta_F::Array{Float64,1},fixed::Array{Float64,2})
    #=
    The function takes fourteen arguments:
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
    + CPI =P=,
    + total profit =total_profit=,
    + total cost =total_cost=,
    + labor share share_L,
    + total fixed costs (without =w= weight) total_fixed_noW,
    + total foreign input share totFshare
    =#

    ################################ Parameters ###################################

    mu = sigma / (sigma - 1)

    ############################## Allocate memory ################################

    c = zeros(num,1)
    share = zeros(maxcard,num)
    share_F = zeros(num,1)
    share_L = zeros(num,1)
    totFshare = zeros(num,1)
    trade = zeros(Int64,num,1)
    indegree = zeros(Int64,num,1)
    buy = zeros(Int64,maxcard,num)
    IndI = zeros(Int64,num,1)
    #FirmID = zeros(Int64,maxcard,num)
    buyer_list = [];

    ############################# Problem of firm 1 ###############################
    Theta = zeros(1,2)
    profit = zeros(1,4)
    ## sourcing capability
    Theta[1,1] = alpha_F[1,1]^(rho-1)*p_F[1,1]^(1-rho) + w^(1-rho) #import
    Theta[1,2] = w.^(1-rho); #don't import

    #=
    Four cases:
    1. Import and export
    2. Import, no export
    3. Export, no import
    4. No export, no import
    =#
    profit[1,1] = (1/sigma) * mu^(1-sigma) * phi[1,1]^(sigma-1) * Theta[1,1]^((sigma-1)/(rho-1)) * (beta_F[1,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w * (2*fixed[1,1])
    profit[1,2] = (1/sigma) * mu^(1-sigma) * phi[1,1]^(sigma-1) * Theta[1,1]^((sigma-1)/(rho-1)) - w * (fixed[1,1])
    profit[1,3] = (1/sigma) * mu^(1-sigma) * phi[1,1]^(sigma-1) * Theta[1,2]^((sigma-1)/(rho-1)) * (beta_F[1,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w * (fixed[1,1])
    profit[1,4] = (1/sigma) * mu^(1-sigma) * phi[1,1]^(sigma-1) * Theta[1,2]^((sigma-1)/(rho-1))
    indegree[1,1]=0
    trade[1,1] = findmax(profit[1,:])[2][1] # findmax() returns (max, CartesianIndex(row, column))
    IndI[1,1] = trade[1,1] < 3 # indicator for importing
    c[1,1] = Theta[1,Int64(2-IndI[1,1])]^(1/(1-rho)) / phi[1,1]; # cost function
    share_L[1,1] = w^(1-rho) / Theta[1,Int64(2-IndI[1,1])]
    share_F[1,1] = IndI[1,1] * ((p_F[1,1]/alpha_F[1,1])^(1-rho)) / Theta[1,Int64(2-IndI[1,1])]
    totFshare[1,1] = share_F[1,1];

    ############################# Problem of firm 2 ###############################
    Theta = zeros(2,2)
    profit = zeros(2,4)

    suppliers = sample(N[1:2, 1], N[2,2]; replace=false, ordered=false) # random draw of length N[j,2] from preceding firms
    suppliers_sorted = suppliers[sortperm(phi[suppliers,1], rev=true)] # because fixed costs are identical, can go in order of productivty (=lowest cost) through eligible supplier set

    Theta[1,1] = alpha_F[2,1]^(rho-1) * p_F[2,1]^(1-rho) + w^(1-rho) #import
    Theta[1,2] = w.^(1-rho) #no import
    Theta_cumsum = alpha_H[2,1]^(rho-1) * c[1,1]^(1-rho)

    Theta[2,1] = Theta_cumsum + alpha_F[2,1]^(rho-1) * p_F[2,1]^(1-rho) + w^(1-rho) # import
    Theta[2,2] = Theta_cumsum + w^(1-rho) #no import

    # profit[1,:] = profit without any indegree
    profit[1,1] = (1/sigma) * mu^(1-sigma) * phi[2,1]^(sigma-1) * Theta[1,1]^((sigma-1)/(rho-1)) * (beta_F[2,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w * (2*fixed[2,1])
    profit[1,2] = (1/sigma) * mu^(1-sigma) * phi[2,1]^(sigma-1) * Theta[1,1]^((sigma-1)/(rho-1)) - w*fixed[2,1]
    profit[1,3] = (1/sigma) * mu^(1-sigma) * phi[2,1]^(sigma-1) * Theta[1,2]^((sigma-1)/(rho-1)) * (beta_F[2,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w * fixed[2,1]
    profit[1,4] = (1/sigma) * mu^(1-sigma) * phi[2,1]^(sigma-1) * Theta[1,2]^((sigma-1)/(rho-1))

    # profit[2,:] = profit for all possible lengths of in degree
    profit[2,1] = (1/sigma) * mu^(1-sigma) * phi[2,1]^(sigma-1) * Theta[2,1].^((sigma-1)/(rho-1)) * (beta_F[2,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w * (fixed[2,1] .+ 2*fixed[2,1])
    profit[2,2] = (1/sigma) * mu^(1-sigma) * phi[2,1]^(sigma-1) * Theta[2,1].^((sigma-1)/(rho-1)) - (w * fixed[2,1] .+ fixed[2,1])
    profit[2,3] = (1/sigma) * mu^(1-sigma) * phi[2,1]^(sigma-1) * Theta[2,2].^((sigma-1)/(rho-1)) * (beta_F[2,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w*(fixed[2,1] .+ fixed[2,1])
    profit[2,4] = (1/sigma) * mu^(1-sigma) * phi[2,1]^(sigma-1) * Theta[2,2].^((sigma-1)/(rho-1)) - w * fixed[2,1]

    # get profit-maximizing sourcing strategy
    PiMaxInx = findmax(profit)[2]
    indegree[2,1] = PiMaxInx[1] - 1
    trade[2,1] = PiMaxInx[2]

    # fill in buyer list and shares
    if indegree[2,1] == 0
        IndI[2,1] = trade[2,1] < 3
        c[2,1] = Theta[1, Int64(2-IndI[2,1])]^(1/(1-rho)) / phi[2,1]
        share_L[2,1] = w^(1-rho) / Theta[1, Int64(2-IndI[2,1])]
        share_F[2,1] = IndI[2,1] * ((p_F[2,1]/alpha_F[2,1])^(1-rho)) / Theta[1, Int64(2-IndI[2,1])]
        share[1,2] = 0
        totFshare[2,1] = share_F[2,1];
    else
        buy[1:indegree[2,1],2] = suppliers_sorted[1:indegree[2,1]]
        buyer_list = [buyer_list; 2*ones(indegree[2,1],1)]
        #FirmID[1:indegree[2,1],2] = 2*ones(indegree[2,1],1)
        IndI[2,1] = trade[2,1] < 3
        c[2,1] = Theta[indegree[2,1]+1, Int64(2-IndI[2,1])]^(1/(1-rho)) / phi[2,1]
        share_L[2,1] = w^(1-rho) / Theta[indegree[2,1]+1, Int64(2-IndI[2,1])]
        share_F[2,1] = IndI[2,1] * ((p_F[2,1]/alpha_F[2,1])^(1-rho)) / Theta[indegree[2,1]+1, Int64(2-IndI[2,1])]
        share[1:indegree[2,1],2] = alpha_H[2,1]^(rho-1) * c[buy[1:indegree[2,1],2]].^(1-rho) ./ Theta[indegree[2,1]+1, Int64(2-IndI[2,1])]
        totFshare[2,1] = share_F[2,1] + totFshare[buy[1:indegree[2,1],2]]' * share[1:indegree[2,1],2];
    end


    ####################### Problem of firms 3 through N ###########################
    for j=3:num  # i is order number
        suppliers = sample(N[1:j-1, 1], N[j,2]; replace=false, ordered=false) # random draw of length N[j,2] from preceding firms
        suppliers_sorted = suppliers[sortperm(phi[suppliers,1], rev=true)]
        # because fixed costs are identical, can go in order of productivty (=lowest cost) through eligible supplier set
        Theta = zeros(length(suppliers_sorted)+1,2)
        profit = zeros(length(suppliers_sorted)+1,4)

        Theta[1,1] = alpha_F[j,1]^(rho-1) * p_F[j,1]^(1-rho) + w^(1-rho) #import
        Theta[1,2] = w.^(1-rho) #no import
        Theta_cumsum = cumsum(alpha_H[j,1]^(rho-1) * c[suppliers_sorted,1]'.^(1-rho), dims=2)

        Theta[2:length(suppliers_sorted)+1,1] = Theta_cumsum .+ (alpha_F[j,1]^(rho-1) * p_F[j,1]^(1-rho)) .+ w^(1-rho) # import
        Theta[2:length(suppliers_sorted)+1,2] = Theta_cumsum .+ w^(1-rho) #no import

        # profit[1,:] = profit without any indegree
        profit[1,1] = (1/sigma) * mu^(1-sigma) * phi[j,1]^(sigma-1) * Theta[1,1]^((sigma-1)/(rho-1)) * (beta_F[j,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w * (2*fixed[j,1])
        profit[1,2] = (1/sigma) * mu^(1-sigma) * phi[j,1]^(sigma-1) * Theta[1,1]^((sigma-1)/(rho-1)) - w*fixed[j,1]
        profit[1,3] = (1/sigma) * mu^(1-sigma) * phi[j,1]^(sigma-1) * Theta[1,2]^((sigma-1)/(rho-1)) * (beta_F[j,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w * fixed[j,1]
        profit[1,4] = (1/sigma) * mu^(1-sigma) * phi[j,1]^(sigma-1) * Theta[1,2]^((sigma-1)/(rho-1))

        # profit[2:length(suppliers_sorted)+1,:] = profit for all possible lengths of in degree
        # we can just add one new supplier to the list each time because we already sorted eligible suppliers by productivty and all fixed costs are the same
        ## import, export
        profit[2:length(suppliers_sorted)+1,1] = (1/sigma) * mu^(1-sigma) * phi[j,1]^(sigma-1) * Theta[2:length(suppliers_sorted)+1,1].^((sigma-1)/(rho-1)) * (beta_F[j,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w * (cumsum(ones(length(suppliers_sorted),1), dims=1) * fixed[j,1] .+ 2*fixed[j,1])

        ## import, no export
        profit[2:length(suppliers_sorted)+1,2] = (1/sigma) * mu^(1-sigma) * phi[j,1]^(sigma-1) * Theta[2:length(suppliers_sorted)+1,1].^((sigma-1)/(rho-1)) - w * (cumsum(ones(length(suppliers_sorted),1), dims=1) * fixed[j,1] .+ fixed[j,1])

        ## no import, export
        profit[2:length(suppliers_sorted)+1,3] = (1/sigma) * mu^(1-sigma) * phi[j,1]^(sigma-1) * Theta[2:length(suppliers_sorted)+1,2].^((sigma-1)/(rho-1)) * (beta_F[j,1]^(sigma-1) * tau^(1-sigma) * E_F/(P_F^(1-sigma))) - w*(cumsum(ones(length(suppliers_sorted),1), dims=1) * fixed[j,1] .+ fixed[j,1])

        ## no import, no export
        profit[2:length(suppliers_sorted)+1,4] = (1/sigma) * mu^(1-sigma) * phi[j,1]^(sigma-1) * Theta[2:length(suppliers_sorted)+1,2].^((sigma-1)/(rho-1)) - w * (cumsum(ones(length(suppliers_sorted),1), dims=1) * fixed[j,1])

        # get profit-maximizing sourcing strategy
        PiMaxInx = findmax(profit)[2]
        indegree[j,1] = PiMaxInx[1] - 1
        trade[j,1] = PiMaxInx[2]

        # fill in buyer list and shares
        if indegree[j,1] == 0
            IndI[j,1] = trade[j,1] < 3
            c[j,1] = Theta[1, Int64(2-IndI[j,1])]^(1/(1-rho)) / phi[j,1]
            share_L[j,1] = w^(1-rho) / Theta[1, Int64(2-IndI[j,1])]
            share_F[j,1] = IndI[j,1] .* ((p_F[j,1]./alpha_F[j,1])^(1-rho)) ./ Theta[1, Int64(2-IndI[j,1])]
            share[1,j] = 0
            totFshare[j,1] = share_F[j,1]
        else
            buy[1:indegree[j,1],j] = suppliers_sorted[1:indegree[j,1]]
            #FirmID[1:indegree[j,1],j] = j*ones(indegree[j,1],1)
            buyer_list = [buyer_list; j*ones(indegree[j,1],1)]
            IndI[j,1] = trade[j,1] < 3
            c[j,1] = Theta[indegree[j,1]+1, Int64(2-IndI[j,1])]^(1/(1-rho)) / phi[j,1]
            share_L[j,1] = w^(1-rho) / Theta[indegree[j,1]+1, Int64(2-IndI[j,1])]
            share_F[j,1] = IndI[j,1] .* ((p_F[j,1]./alpha_F[j,1])^(1-rho)) ./ Theta[indegree[j,1]+1, Int64(2-IndI[j,1])]
            share[1:indegree[j,1],j] = alpha_H[j,1]^(rho-1) .* c[buy[1:indegree[j,1],j]].^(1-rho) ./ Theta[indegree[j,1]+1, Int64(2-IndI[j,1])]
            totFshare[j,1] = share_F[j,1] + totFshare[buy[1:indegree[j,1],j]]' * share[1:indegree[j,1],j]
        end
    end

    # Construct sparse network (value 1 for each link)
    buy_vector = reshape(buy,(num*maxcard,1))
    buy_vector = buy_vector[buy_vector .> 0]
    buy_vector = [buy_vector; num]
    buyer_list = dropdims(buyer_list, dims = tuple(findall(size(buyer_list) .== 1)...))
    buyer_list = [buyer_list; num]
    SparseNet = sparse(buy_vector, buyer_list, [ones(length(buy_vector)-1); 0]);

    # Construct sparse share matrix
    share_vector = reshape(share,(num*maxcard,1))
    share_vector = share_vector[share_vector .> 0]
    share_sparse = sparse(buy_vector, buyer_list, [share_vector; 0]);

    # Prices and equilibrium
    p_H = mu * c
    P = (sum(p_H.^(1-sigma)))^(1/(sigma-1)); # TKMD eq (5)

    IndE = trade .% 2 .== 1

    E = A*P^(1-sigma)   # domestic expenditure given A
    x_H = p_H * A    # domestic final demand vector
    x_F = tau * p_H .* beta_F.^(sigma-1) * (E_F/P_F^(1-sigma))   # foreign final demand vector

    c_H = ((sigma-1)/sigma) * x_H    # vector: cost of inputs sold to domestic final demand
    c_F = ((sigma-1)/sigma) * x_F    # vector: cost of inputs sold to foreign final demand

    total_cost = gmres(sparse(I,num,num) - share_sparse, c_H + c_F)
    total_sale = x_H + x_F + total_cost - c_H
    total_profit = total_sale - total_cost

    total_fixed_noW = fixed' * SparseNet + IndE'

    return P, total_profit, total_cost, share_L, total_fixed_noW, totFshare
end
