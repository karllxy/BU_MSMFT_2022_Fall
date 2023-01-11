#using Array
using Statistics
S = 100
K = 95
r = 0.01
sigma = 0.1
t = 1
N = 500
P = 10000

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
    #I = V .!= 0
    A = [x^d for d = 1:1, x in X[n+1, :]] # coeff 
    β = A' \ V #
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


max(mean(V), K - S)
