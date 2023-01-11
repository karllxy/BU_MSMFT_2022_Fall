# MF850 Homework1
# Edited by XuyangLiu
using Plots
using Distributions

S0 = 100
K = 95
T = 1
sigma = 0.1
r = 0.02
D = Normal(0,1)
# First we calculate the price obtained via BS model.
function BSmodel()
    d1 = (log(S0/K)+(r + (sigma^2)/2)*T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
    put = K * exp(-r*T) * cdf(D,-d2) - S0 * cdf(D,-d1)

    return put
end

put = BSmodel()


# Problem 1.1 (a)
function put_binomial(N,type)
    u = exp(sigma * sqrt(1/N))
    d = 1/u
    q = (exp(r/N)-d)/(u-d)
    last_price = zeros(N+1)
    last_put_price = zeros(N+1)
    last_price[1] = S0 * u^N
    last_put_price[1] = max(K - last_price[1],0)

    for i in 2:N+1
        last_price[i] = last_price[i-1] * (d^2)
        last_put_price[i] = max(K - last_price[i], 0)
    end

    if type == "European"
        for j in N+1:-1:1
            for i in 1:(j-1)
                last_put_price[i]  = (q * last_put_price[i] + (1-q) * last_put_price[i+1]) * exp(-r*T/N)
            end
        end
    end

    if type == "American"
        for j in N+1:-1:1
            for i in 1:(j-1)
                last_price[i] = last_price[i+1] * u
                last_put_price[i]  = max((q * last_put_price[i] + (1-q) * last_put_price[i+1]) * exp(-r*T/N),K - last_price[i])
            end
        end
    end

    return last_put_price[1]

end

Problem1a = put_binomial(1000,"European")
print("The answer to the Problem1.1 (a) is about: ",Problem1a)
#Problem 1.1 (b)
# We can also use the function we built in (a)

function path_converg(N,type)
    path = zeros(N)
    for i in 1:N
        path[i] = put_binomial(i,type)
    end

    return path

end
pathN = 300
Europath = path_converg(pathN,"European")
print("\n")
print(put)
print("\n")
Euro_delta =  abs.(Europath.-put)
# From the plot we can see that it converge to BS price.
# Therefore, binomial-Bsprice converges to 0
x = 1:1:pathN
y = [ 1/x for x in x ]
plot(x,[Euro_delta y],label=["Difference of Binomial Price and BS Price" "1/n"],title="Convergence analysis for Euro-put")
# From the plot we can consider 1/n as the e(n).


# Problem1.1 (c)
# We repeat what we have done above, but in American rules.
Problem1c = put_binomial(1000,"American")
Ampath = path_converg(pathN,"American")
Ame_delta = abs.(Ampath.-put)
# e(n) = 1/n + 0.06
y2 = [ (1/x + 0.06) for x in x ]
plot(x,[Ame_delta y2],label=["Difference of Binomial Price and BS Price" "1/n + 0.06"],title="Convergence analysis for American-put")

#According to order of convergence, this sequence is linearly to 0.06,and the rate of convergence is still 1.
