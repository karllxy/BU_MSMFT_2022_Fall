# MF850 Homework2
# Edited by XuyangLiu
using Plots
using SparseArrays
using LinearAlgebra

function analytical(x) # analytical function for 2.1(a)
    c2 = 7 * exp(1) / (3 * (exp(1) - 1))
    c1 = -c2
    y = c1 + c2 * exp(-x) + (x^3) / 3 - x^2 + 2 * x
    return y
end

function analytical2(x) # analytical function for 2.1(d)
    c2 = -5 / 3 * exp(1)
    c1 = -1
    y = c1 + c2 * exp(-x) + (x^3) / 3 - x^2 + 2 * x
    return y
end

# Problem2.1 (a)

function descri_ODE(N, condition)
    h = 1 / N
    lst = 0:h:1


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

    if condition == "Dirichlet"

        tg = zeros(length(lst)) # x^2 list
        analy = zeros(length(lst)) # the analytical result
        for i = 1:length(lst)
            analy[i] = analytical(lst[i])
            tg[i] = lst[i]^2
        end

        # v''+v'-x^2
        A = A1 + A2
        A[1, 1] = 1
        A[N+1, N+1] = 1

        # modify x^2
        f = zeros(N + 1)
        f[1] = 0
        f[N+1] = -1
        y = f + tg
        y[N+1] = -1

        # then we can calculate the V
        V = A^-1 * y
        V[1] = 0
    end

    if condition == "Robin"

        tg = zeros(length(lst)) # x^2 list
        analy = zeros(length(lst)) # the analytical result
        for i = 1:length(lst)
            analy[i] = analytical2(lst[i])
            tg[i] = lst[i]^2
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
        f[1] = 2 * 1 / (b0 * h) - 2 * 1 / b0
        f[end] = 2 * 0 + 2 * 0
        y = f + tg

        V = A^-1 * y
    end
    return V, analy, lst
end

# We plot this two curve
N = 1000
Va, analya, lsta = descri_ODE(N, "Dirichlet")

plot(lsta, [Va, analya], label = ["Descretization" "Analytical"])

#Porblem2.1(b) (c)
error = zeros(1000)
for i = 1:1000
    vb, anab, lstb = descri_ODE(i, "Dirichlet")
    er = sum(vb - anab)
    error[i] = abs(er)
end
plot(error[5:end]) # the error plot

# From the plot we can see that its a sublinear
# order = abs(log(abs((c-b)/(b-a))) / log(abs((b-a)/(a-z)))) 
# the order is about 0.7 , close to 1. So its a linear convergence, to be specific, a sublinear


# Problem 2.1 (D)
vd, analyd, lstd = descri_ODE(1000, "Robin")

# For Problem 2.2, please see the another code file problem2.2.jl
include("problem2.2.jl")
