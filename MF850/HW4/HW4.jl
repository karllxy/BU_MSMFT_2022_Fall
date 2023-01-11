# MF850 Homework 4
# Edited by Xuyang Liu
using Statistics
# Global Variables
r = 0.01
sigma = 0.1
T = 1
K = 95
N = 100000

# Auxiliary Function
function mc_sim(S0, r, sigma, T, n, N)
    dt = T / n
    Sd = S0 .* exp.(sigma * dt^0.5 * randn(N) .+ (r - sigma^2 / 2) * dt)
    return Sd
end

function maxim(mat, n)
    mat2 = zeros(N, n)
    for i = 1:length(mat[1:end, 1])
        for j = 1:length(mat[1, 1:end])
            mat2[i, j] = max(K - mat[i, j], 0)
        end
    end
    return mat2
end
function ls_beta(x, y) # might not be used
    beta = (x' * x)^-1 * x' * y
    return beta
end

# Problem (a)
X0 = 100
n1 = 1
function one_step(X0, n1, tp)
    S_mat = zeros(N, n1)
    mc = mc_sim(X0, r, sigma, T, n1, N)
    St = mc
    S_mat[1:end, 1] .= mc

    V_mat = maxim(S_mat, n1)
    #beta1 = ls_beta(St,V_mat)
    X = ones(N, 2)
    X[1:end, 2] = St
    beta1 = X \ V_mat   # this is using backslash to run the linear regression
    V_pre = X * beta1
    if tp == "LS"
        P0 = mean(V_pre * exp(-r * T))
    end
    if tp == "TvR"
        P0 = mean(V_mat * exp(-r * T))
    end
    return P0
end

# in this way, we can have:
P1_LS = one_step(X0, n1, "LS")
P1_TvR = one_step(X0, n1, "TvR")
P1_LS, P1_TvR





# Problem (b)(c)
S = 100
K = 95
r = 0.01
sigma = 0.1
t = 1
N = 5
P = 10000

function LS_MT(N, type)
    dt = t / N
    R = exp(r * dt)
    T = typeof(S * exp(-sigma^2 * dt / 2 + sigma * dt^0.5 * 0.1) / R)
    X = zeros(N + 1, P)

    for p = 1:P
        X[1, p] = x = S
        for n = 1:N
            x *= R * exp(-sigma^2 * dt / 2 + sigma * dt^0.5 * randn())
            X[n+1, p] = x
        end
    end

    V = [max(K - x, 0) / R for x in X[N+1, :]]

    if type == "exclude"
        for n = N-1:-1:1
            I = V .!= 0
            A = [x^d for d = 0:1, x in X[n+1, :]] # coeff 
            β = A[:, I]' \ V[I] #
            cV = A' * β
            for p = 1:P
                ev = max(K - X[n+1, p], 0)
                if I[p] && cV[p] < ev
                    V[p] = ev / R
                else
                    V[p] /= R
                end
            end
        end
    end

    if type == "include"
        for n = N-1:-1:1
            A = [x^d for d = 0:1, x in X[n+1, :]] # coeff 
            β = A' \ V
            cV = A' * β
            for p = 1:P
                ev = max(K - X[n+1, p], 0)
                if cV[p] < ev
                    V[p] = ev / R
                else
                    V[p] /= R
                end
            end
        end
    end

    return max(mean(V), K - S)
end

LS_in0 = LS_MT(15, "include")
LS_ex0 = LS_MT(15, "exclude")

function TvR_MT(N, type)
    dt = t / N
    R = exp(r * dt)
    T = typeof(S * exp(-sigma^2 * dt / 2 + sigma * dt^0.5 * 0.1) / R)
    X = zeros(N + 1, P)

    for p = 1:P
        X[1, p] = x = S
        for n = 1:N
            x *= R * exp(-sigma^2 * dt / 2 + sigma * dt^0.5 * randn())
            X[n+1, p] = x
        end
    end

    V = [max(K - x, 0) / R for x in X[N+1, :]]

    if type == "exclude"
        for n = N-1:-1:1
            I = V .!= 0
            A = [x^d for d = 0:1, x in X[n+1, :]] # coeff 
            β = A[:, I]' \ V[I] #
            cV = A' * β
            for p = 1:P
                ev = max(K - X[n+1, p], 0)
                if I[p] && cV[p] < ev
                    V[p] = ev / R
                else
                    V[p] = cV[p] / R
                end
            end
        end
    end

    if type == "include"
        for n = N-1:-1:1
            A = [x^d for d = 0:1, x in X[n+1, :]] # coeff 
            β = A' \ V
            cV = A' * β
            for p = 1:P
                ev = max(K - X[n+1, p], 0)
                if cV[p] < ev
                    V[p] = ev / R
                else
                    V[p] = cV[p] / R
                end
            end
        end
    end

    return max(mean(V), K - S)
end

TvR_in0 = TvR_MT(15, "include")
TvR_ex0 = TvR_MT(15, "exclude")

LS_in0, LS_ex0, TvR_in0, TvR_ex0
