# MF850 Homework 4
# Edited by Xuyang Liu
using Statistics

# Global Variables
S = 100
K = 95
r = 0.01
sigma = 0.1
t = 1
N = 15
P = 10000

function LS_MT(N, type)
    dt = t / N
    R = exp(r * dt)
    T = typeof(S * exp(-sigma^2 * dt / 2 + sigma * √dt * 0.1) / R)
    X = zeros(N + 1, P)

    for p = 1:P
        X[1, p] = x = S
        for n = 1:N
            x *= R * exp(-sigma^2 * dt / 2 + sigma * √dt * randn())
            X[n+1, p] = x
        end
    end

    V = [max(K - x, 0) / R for x in X[N+1, :]]

    if type == "exclude"
        for n = N-1:-1:1
            I = V .!= 0
            A = [x^d for d = 0:5, x in X[n+1, :]] # coeff 
            β = A[:, I]' \ V[I] #
            cV = A' * β
            for p = 1:P
                ev = max(K - X[n+1, p], 0)
                if I[p] && cV[p] < ev
                    V[p] = ev / R
                else
                    V[p] = V[p] / R
                end
            end
        end
    end

    if type == "include"
        for n = N-1:-1:1
            A = [x^d for d = 0:5, x in X[n+1, :]] # coeff 
            β = A' \ V
            cV = A' * β
            for p = 1:P
                ev = max(K - X[n+1, p], 0)
                if I[p] && cV[p] < ev
                    V[p] = ev / R
                else
                    V[p] = V[p] / R
                end
            end
        end
    end

    return max(mean(V), K - S)
end



function TvR_MT(N)
    dt = t / N
    R = exp(r * dt)
    T = typeof(S * exp(-sigma^2 * dt / 2 + sigma * √dt * 0.1) / R)
    X = zeros(N + 1, P)

    for p = 1:P
        X[1, p] = x = S
        for n = 1:N
            x *= R * exp(-sigma^2 * dt / 2 + sigma * √dt * randn())
            X[n+1, p] = x
        end
    end

    V = [max(K - x, 0) / R for x in X[N+1, :]]

    for n = N-1:-1:1
        A = [x^d for d = 0:5, x in X[n+1, :]] # coeff 
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
    return max(mean(V), K - S)
end

#(a)
Pa_ls = LS_MT(1, "include")
Pa_TvR = TvR_MT(1)
Pa_ls, Pa_TvR

#(b)
Pb_ls = LS_MT(252, "include")
Pb_TvR = TvR_MT(252)
Pb_ls, Pb_TvR

#(c)
Pc1 = LS_MT(252, "include")
Pc2 = LS_MT(252, "exclude")
Pc1, Pc2
