# MF850 Homework 3 
# Edited by XuyangLiu
using Plots
using SparseArrays
using LinearAlgebra

# Global Variables:
sigma = 0.1
k = 0.0001
fi = 10 * k
alpha = 100 * k
S_low = 20
S_high = 20.1

N = 100
delta = (S_high - S_low) / N
h0 = S_low:delta:S_high
hl = zeros(length(h0))
for i = 1:length(h0)
    hl[i] = h0[i]
end

A = zeros(N + 1, N + 1)
A[1, 1] = 1
A[end, end] = 1
for i = 2:length(h0)-1
    A[i, i-1] = 1 / delta^2
    A[i, i] = -2 / delta^2
    A[i, i+1] = 1 / delta^2
end

fil = ones(N + 1) * fi * -1
fil[1] = (fi * k)^0.5
fil[end] = alpha

function evaluate(h)
    hA = zeros(N + 1, N + 1)
    for i = 1:N+1
        hA[i, i] = h[i]
    end

    NA = 0.5 * sigma^2 * A - 1 / k * hA
    NA[1] = 1
    NA[end] = 1

    new_ha = NA^-1 * fil


    return new_ha

end

function iter(n)
    h = hl
    for i = 1:n
        h = evaluate(h)
    end
    return h
end

it = iter(2000)
