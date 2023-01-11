using Plots
using SparseArrays
using LinearAlgebra

function sur_plot(n)
    h = 2/n
    p(x,y) = cos.(pi.*(x.+y.^2))
    A = zeros((n+1)*(n+1),(n+1)*(n+1))
    y = zeros((n+1)*(n+1),1)
    for i = 1:n+1
        A[i,i] = 1
        A[(n+1)*i,(n+1)*i] = 1
        A[n*(n+1)+i,n*(n+1)+i] = 1
        A[(n+1)*(i-1)+1,(n+1)*(i-1)+1] = 1
        y[i] = p(-1,(i-1)*h-1)
        y[(n+1)*i] = p((i-1)*h-1,1)
        y[n*(n+1)+i] = p(1,(i-1)*h-1)
        y[(n+1)*(i-1)+1] = p((i-1)*h-1,-1)
    end
    for i = 2:n,j = 2:n
        A[(n+1)*(i-1)+j,(n+1)*(i-1)+j] = -4/h^2
        A[(n+1)*(i-1)+j,(n+1)*(i-1)+j-1] = 1/h^2
        A[(n+1)*(i-1)+j,(n+1)*(i-1)+j+1] = 1/h^2
        A[(n+1)*(i-1)+j,(n+1)*i+j] = 1/h^2
        A[(n+1)*(i-1)+j,(n+1)*(i-2)+j] = 1/h^2
    end
    reshape(inv(A)*y,n+1,:)
end

sur_plot(100)