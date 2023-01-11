module HWP1

export solve, Dirichlet, Neumann, Robin

using LinearAlgebra: Tridiagonal

struct Domain
    x::AbstractRange{Float64}
    N::Int
    h::Float64
end

function Domain(xmin::Float64, xmax::Float64, N::Int)
    x = range(xmin, stop = xmax, length = N)
    h = step(x)
    Domain(x, N, h)
end

abstract type BoundaryCondition end

struct Dirichlet <: BoundaryCondition
    f::Float64
end

struct Neumann <: BoundaryCondition
    f::Float64
end

struct Robin <: BoundaryCondition
    a::Float64
    b::Float64
    f::Float64
end


# Implementation of Dirichlet
A_lower(condition::Dirichlet, dom) = [Float64(i == 1) for i = 1:dom.N]
A_upper(condition::Dirichlet, dom) = [Float64(i == dom.N) for i = 1:dom.N]
fhat_lower(condition::Dirichlet, dom) = condition.f / dom.h^2
fhat_upper(condition::Dirichlet, dom) = condition.f / dom.h^2


# Implementation of Robin
function A_lower(condition::Robin, dom)
    res = zeros(dom.N)
    a, b = condition.a, condition.b
    h = dom.h
    res[1] = -(2 - 2h * a / b + h^2 * a / b)
    res[2] = 2
    res
end

function A_upper(condition::Robin, dom)
    res = zeros(dom.N)
    a, b = condition.a, condition.b
    h = dom.h
    res[end-1] = 2
    res[end] = -(2 + 2h * a / b + h^2 * a / b)
    res
end

fhat_lower(condition::Robin, dom) =
    - condition.f * (1 - 2 / dom.h) / condition.b
fhat_upper(condition::Robin, dom) =
    condition.f * (1 + 2 / dom.h) / condition.b


# Implementation of Neumann, using Robin
A_lower(condition::Neumann, N) = A_lower(Robin(0, 1, condition.f), N)
A_upper(condition::Neumann, N) = A_upper(Robin(0, 1, condition.f), N)

fhat_lower(condition::Neumann, dom) =
    fhat_lower(Robin(0, 1, condition.f), dom)
fhat_upper(condition::Neumann, dom) =
    fhat_upper(Robin(0, 1, condition.f), dom)



function solve(N::Int, boundary_l::BoundaryCondition, boundary_u::BoundaryCondition)
    xmin = 0.0
    xmax = 1.0
    dom = Domain(xmin, xmax, N)
    h = dom.h

    ld = fill(1 - h / 2, N - 1)
    ud = fill(1 + h / 2, N - 1)
    d = fill(-2.0, N)
    A = Tridiagonal(ld, d, ud)
    A[1, :] = A_lower(boundary_l, dom)
    A[end, :] = A_upper(boundary_u, dom)

    fhat = collect(dom.x.^2)
    fhat[1] += fhat_lower(boundary_l, dom)
    fhat[end] += fhat_upper(boundary_u, dom)

    dom.x, (A \ fhat) .* h^2
end

end
