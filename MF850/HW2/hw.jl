using Plots, LaTeXStrings


## P1
begin

include("HWP1.jl")
import .HWP1 as P1
using Statistics: mean

# a
v_true(x) = (x^3 - 3 * x^2 + 6 * x) / 3 + 7 * (1 - exp(- x))/(3 * exp(-1) - 3)
#x^2 / 2 - x - (exp(-x) - 1) / (2 * (exp(-1) - 1))
x, v = P1.solve(100, P1.Dirichlet(0), P1.Dirichlet(-1))
plot(x, v, label="Numeric");
plot!(x, v_true, label="Analytic")

# b
xvs_len = 10
# The added +1 to the number of points makes sure h is halved for each i. Array of (x, v) tuples
xvs = [P1.solve(10 * 2^i + 1, P1.Dirichlet(0), P1.Dirichlet(-1)) for i = 0:(xvs_len-1)]
# Below, last takes v from the tuple (x,v), the zipped iterator is equivalend to getting v[i],v[i+1],v[i+2], the indices [1:2:end] and [1:4:end] make sure to select coinciding grid points despite different grid sizes.  Finally, 2:end-1 skips the first and last points, because it causes division by 0.
conv_order = [
    log2.(abs.((v1 - v2[1:2:end]) ./ (v2[1:2:end] - v4[1:4:end]))[2:end-1]) for
    (v1, v2, v4) in (v -> zip(v, v[2:end], v[3:end]))(last.(xvs))
]
# Almost all are 2, but degrade for low h due to the condition number of the matrix
mean.(conv_order)
minimum.(conv_order)

# c
vs_err = [maximum(abs.(v - v_true.(x))) for (x, v) in xvs]
plot(
    step.(first.(xvs)),
    vs_err,
    xaxis = (L"h", :log),
    yaxis = (L"L^\infty", :log),
    legend = :none,
)

# d
v_true(x) = 1/3 * (x^3 - 3x^2 + 6x - 5 * exp(1 - x) - 3)
x, v = P1.solve(100, P1.Robin(1, 1, 1), P1.Robin(2, 1, 0))
plot(x, v);
plot!(x, v_true)

end



## P2
begin
include("HWP2.jl")
import .HWP2 as P2

v = P2.solve()
surface(P2.xs, P2.ys, v', xaxis = L"x", yaxis = L"y", zaxis = L"v")
# This approximates v(-0.5, 0) and converges as the grid size tends to zero
v[findfirst(x -> x > -0.5, P2.xs), findfirst(y -> y > 0.0, P2.ys)]
end
