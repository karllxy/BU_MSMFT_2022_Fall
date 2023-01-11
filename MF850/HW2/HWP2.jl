module HWP2

const Nx = 200
const Ny = 200
const xs = range(-1, stop = 1, length = Nx)
const ys = range(-1, stop = 1, length = Ny)
const hx = step(xs)
const hy = step(ys)

ind(i, j) = (j - 1) * Nx + i

boundary_value(x, y) = cos(π * (x + y^2))

function RHS()
    RHS = zeros(Nx * Ny)
    for j = 1:Ny, i in [1, Nx]
        RHS[ind(i, j)] = boundary_value(xs[i], ys[j])
    end
    for i = 1:Nx, j in [1, Ny]
        RHS[ind(i, j)] = boundary_value(xs[i], ys[j])
    end
    RHS
end

# Dense version.  The struct is just to get convenience with custom indexing
struct DenseCoefficients
    arr::Array{Float64, 2}
end
DenseCoefficients(N) = DenseCoefficients(zeros(N, N))
Base.setindex!(A::DenseCoefficients, val, i, j, ii, jj) = A.arr[ind(i, j), ind(ii, jj)] = val
function makeDenseA()
    A = DenseCoefficients(Nx * Ny)
    for j = 2:Ny-1, i = 2:Nx-1
        A[i, j, i,   j  ] = -2 / hx^2 - 2 / hy^2
        A[i, j, i+1, j  ] = 1 / hx^2
        A[i, j, i,   j+1] = 1 / hy^2
        A[i, j, i-1, j  ] = 1 / hx^2
        A[i, j, i,   j-1] = 1 / hy^2
    end
    for j = 1:Ny, i in [1, Nx]
        A[i, j, i, j] = 1
    end
    for i = 1:Nx, j in [1, Ny]
        A[i, j, i, j] = 1
    end
    A.arr
end

# Sparse version
import SparseArrays
mutable struct SparseComponents
    current_index::Int
    I::Array{Int,1}
    J::Array{Int,1}
    V::Array{Float64,1}
end
SparseComponents(len) = SparseComponents(0, zeros(Int, len), zeros(Int, len), zeros(len))
function Base.setindex!(A::SparseComponents, val, i, j, ii, jj)
    A.current_index += 1
    A.I[A.current_index] = ind(i, j)
    A.J[A.current_index] = ind(ii, jj)
    A.V[A.current_index] = val
end
sparse(A::SparseComponents; m, n, combine = (x, y) -> y) = @views SparseArrays.sparse(
    A.I[1:A.current_index],
    A.J[1:A.current_index],
    A.V[1:A.current_index],
    m,
    n,
    combine,
)
function makeSparseA()
    A = SparseComponents(5 * Nx * Ny)
    for j = 2:Ny-1, i = 2:Nx-1
        A[i, j, i,   j  ] = -2 / hx^2 - 2 / hy^2
        A[i, j, i+1, j  ] = 1 / hx^2
        A[i, j, i,   j+1] = 1 / hy^2
        A[i, j, i-1, j  ] = 1 / hx^2
        A[i, j, i,   j-1] = 1 / hy^2
    end
    for j = 1:Ny, i in [1, Nx]
        A[i, j, i, j] = 1
    end
    for i = 1:Nx, j in [1, Ny]
        A[i, j, i, j] = 1
    end
    sparse(A, m = Nx * Ny, n = Nx * Ny)
end

# Note that the dense version is never really better, just kept for illustration
makeA() = (Nx * Ny <= 10^3) ? makeDenseA() : makeSparseA()

solve() = reshape(makeA() \ RHS(), Ny, Nx)

end
