using Plots
using SparseArrays
using LinearAlgebra

N = 100
h = (1 + 1) / N
lst = -1:h:1


#T = Tridiagonal(ones(N-1), -4 * ones(N), ones(N-1))
#I = Diagonal(ones(N))

#A = 1/h^2 * Tridiagonal(I,T,I)

A = zeros((N+1)^2, (N+1)^2)
#A = sparse(A)

for i = 1:N+1
    A[i,i] = 1
    A[(N+1)*i,(N+1)*i] = 1
    A[N*(N+1)+i,N*(N+1)+i] = 1
    A[(N+1)*(i-1)+1,(N+1)*(i-1)+1] = 1
end

for i = 2:N, j = 2:N
    A[(N+1)*(i-1)+j,(N+1)*(i-1)+j] = -4/h^2
    A[(N+1)*(i-1)+j,(N+1)*(i-1)+j-1] = 1/h^2
    A[(N+1)*(i-1)+j,(N+1)*(i-1)+j+1] = 1/h^2
    A[(N+1)*(i-1)+j,(N+1)*i+j] = 1/h^2
    A[(N+1)*(i-1)+j,(N+1)*(i-2)+j] = 1/h^2

end

f = zeros((N+1)^2)
bounds(x,y) = cos.(π.*(x.+y.^2))
#f = sparse(f)

for i = 1:N+1
    f[i] = bounds(-1,(i-1)*h-1)
    f[(N+1)*i] = bounds((i-1)*h-1,1)
    f[N*(N+1)+i] = bounds(1,(i-1)*h-1)
    f[(N+1)*(i-1)+1] = bounds((i-1)*h-1,-1)
end

V = A^-1 * f

x = zeros(N+1)
for i = 1:N+1
    x[i] = lst[i]
end

y = x
K = reshape(V, (N+1, N+1))

surface(x, y, K)

# for here, my v(-0.5,0) should be K[50,25]
