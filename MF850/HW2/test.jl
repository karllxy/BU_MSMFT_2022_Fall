using Plots
using SparseArrays
using LinearAlgebra

function analytical(x)
    c2 = -5 / 3 * exp(1)
    c1 = -1
    y = c1 + c2 * exp(-x) + (x^3) / 3 - x^2 + 2 * x
    return y
end

N = 1000
h = 1 / N
lst = 0:h:1

tg = zeros(length(lst)) # x^2 list
analy = zeros(length(lst)) # the analytical result
for i = 1:length(lst)
    analy[i] = analytical(lst[i])
    tg[i] = lst[i]^2
end

# we descretization v'' , v'
A1 = zeros(N + 1, N + 1)
A2 = zeros(N + 1, N + 1)

for i = 2:length(lst)-1
    A1[i, i-1] = -1 / h
    A1[i, i] = 1 / h
    A2[i, i-1] = 1 / (h^2)
    A2[i, i] = -2 / (h^2)
    A2[i, i+1] = 1 / (h^2)
end

a0, b0, an, bn = 1, 1, 2, 1
A1[1, 1] = (1 - 2 * a0 * h / b0) / h
A1[1, 2] = -1 / h
A1[end, end-1] = 1 / h
A1[end, end] = (-1 - 2 * an * h / bn) / h

A2[1, 1] = (-2 + 2 * a0 * h / b0) / (h^2)
A2[1, 2] = 2 / (h^2)
A2[end, end-1] = 2 / (h^2)
A2[end, end] = (-2 - 2 * an * h / bn) / (h^2)

A = A1 + A2

f = zeros(N + 1)
f[1] = 2 * a0 / (b0 * h) - 2 * a0 / b0
f[end] = 2 * 0 + 2 * 0
y = f + tg
y[end] = 0
y[1]
V = A^-1 * y
