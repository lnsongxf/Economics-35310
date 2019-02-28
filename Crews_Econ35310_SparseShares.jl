## Crews_Econ35310_SparseShares

function Crews_SparseShares(N::Array{Float64,2},c::Array{Float64,2},Theta::Array{Float64,2},
    p_F::Array{Float64,2},rho::Float64)
    #=
    The function takes five arguments:
    + a fixed network N,
    + a vector of firm costs c,
    + a vector of Theta,
    + a vector of foreign prices p_F,
    + production elasticity of substitution rho

    Given N, c, and Theta, the function solves for
    + the matrix of input shares S
    that satisfies TKMD equations (8)-(10)
    =#

    # Retrieve number of firms in simulation
    num = size(N)[1]
    cards = N[:,2]
    maxcard = Int64(maximum(cards))

    # Retrieve firm characteristic vectors
    IndI = N[:,maxcard+3]
    IndE = N[:,maxcard+4]
    X = N[:,maxcard+6]
    Y = N[:,maxcard+7]
    alpha_F = N[:,maxcard+8]

    # Initialize sparse matrix
    cnx = Int64.(sum(cards[1:num,1])) # total number of connections
    J = Array{Int64}(zeros(cnx+2*num)) # index of firm with space for all connections + labor & foreign input shares for each firm
    K = Array{Int64}(zeros(cnx+2*num)) # index of share: so (j,k) is s_{kj}
    V = zeros(cnx+2*num) # value at that index pair
    J[1:2] = [1, 1]
    K[1] = num+1
    K[2] = num+2
    V[1] = (w^(1-rho)) / (Theta[1,1])
    V[2] = IndI[1,1]*(alpha_F[1,1])^(rho-1)*(p_F[1,1])^(1-rho) / (Theta[1,1]);


    for j=2:num
        index = findfirst(x -> x == 0, J)
        jcard = Int64.(cards[j,1])
        J[index:index+jcard+1] = repeat([j],jcard+2)
        K[index+jcard] = num+1 #column for labor share
        K[index+jcard+1] = num+2 #column for direct import share
        V[index+jcard] = (w^(1-rho)) / (Theta[j,1])
        V[index+jcard+1] = IndI[j,1]*(alpha_F[j,1])^(rho-1)*(p_F[j,1])^(1-rho) / (Theta[j,1])

        # direct firm-to-firm shares
        Z = N[j,3:maxcard+2]
        filter!(e -> e>0, Z)
        Z = Array{Int64}(Z)
        temp = zeros(length(Z))
        for k = 1:length(Z)
            temp[k] = (X[Z[k],1]*Y[j,1])^(rho-1)*(c[Z[k],1])^(1-rho) / (Theta[j,1]) # = S[k,j]
        end
        if jcard > 0
            K[index:index+jcard-1] = Z
            V[index:index+jcard-1] = temp
        end
    end
    filter!(e -> e>0, J)
    K = K[1:length(J)]
    V = V[1:length(J)]
    S = sparse(J,K,V)

    #sparse network
    SparseNet = sparse(J,K,ones(length(J)))
    SparseNet = SparseNet[:,1:num]
    return S, SparseNet
end
