using Plots # add Plots.jl from the package manager if you have not already done so.

# HW 1 starting script (if you want): contains function computeLU() - to compute an LU-factorization of square
# matrix A, namely, A = LU, where L and U are lower and upper triangular matrices. 


"""
    computeLU(A)
Compute and return LU factorization `LU = A` of square matrix `A`.  
Might not work on all matrices, since no pivoting is done!

# Examples (don't need examples, but fine to include)
'''
julia> A = [6 -2 2;12 -8 6;3 -13 3]
3×3 Array{Int64,2}:
  6   -2  2
 12   -8  6
  3  -13  3
julia> (L, U) = computeLU(A)
([1.0 0.0 0.0; 2.0 1.0 0.0; 0.5 3.0 1.0], [6.0 -2.0 2.0; 0.0 -4.0 2.0; 0.0 0.0 -4.0])
julia> norm(A - L*U)
0.0
'''
"""


function computeLUP(A)
    println("\n", "****ComputeLUP(A):")
    N = size(A)[1]

    #Id = Matrix{Float64}(I, N, N) # N x N identity matrix
    Id = create_identity(N)

    L = copy(Id)   # initialize
    U = copy(Id)   # initialize
    P = copy(Id)   # initialize
    Ã  = copy(A) # initialize. Ã corresponds to A as it goes under elimination stages

    for k = 1:N-1 # march across columns
        println("\n", "Current k: ", k, " with pivot value: ", Ã[k,k] )
        
        println("\n", "\n", "---->ComputeLUP(A) calling Compute_Pk")
        (j, p) = compute_Pk(Ã, k)
        
        #update Ã
        println("\n", "Updating Ã.= p*Ã ")
        Ã .= p * Ã
        show(stdout, "text/plain", Ã)

        println("\n", "Updating P .= p * P ")
        #update permutation matrix
        P .= p * P
        show(stdout, "text/plain", P)

        println("\n", "\n", "---->ComputeLUP(A) calling Compute_Lk")
        (Lk, Lk_inv) = compute_Lk(Ã, k)

        println("\n", "Ã .= Lk * Ã ")
        Ã .= Lk * Ã
        show(stdout, "text/plain", Ã)

        println("\n", "L .= L * Lk_inv ")
        L .= L * Lk_inv
        show(stdout, "text/plain", L)

        if k>1 
            #do the swap
            (L[k, 1:k-1], L[j, 1:k-1]) = (L[j, 1:k-1], L[k, 1:k-1])
        end

    end

    U .= Ã
    

    return (L, U, P)

end

""" 
Finds row j, the row containing largest value in column k, including and below Ãkk.
Also compute permutation matrix pK,J
"""
function compute_Pk(A, k)
    println("\n", "--Inside compute_Pk()--")
    N = size(A)[1]

    Pk = create_identity(N) 
    
    max = A[k, k]
    max_row = k

    # search for the largest value in column
    for row = k:N
        #compare the current column value to the current max value
        if abs(A[row, k]) > abs(max)
                max = A[row, k]
                max_row = row
        end
    end

    #update identity matrix
    (Pk[k, :], Pk[max_row, :]) = (Pk[max_row, :], Pk[k, :]) 
    println("\n", "Row index of max value in column ", k, " is ", max_row, ". And max value is ", max)
    println("\n", "P(", k, ",", max_row, ") = ")
    show(stdout, "text/plain", Pk)
    
    return (max_row, Pk)

end

"""
    compute_Lk(A, k)
Compute Lk and its inverse from A, assuming first k-1 columns have undergone elimination.

"""
function compute_Lk(A, k)
    println("\n", "--Inside compute_Lk()--")
    N = size(A)[1]

    Lk = create_identity(N) # Matrix{Float64}(I, N, N)       # initialize as identity matrix
    Lk_inv = create_identity(N)# Matrix{Float64}(I, N, N)   # initialize as identity matrix

    println("\n", "--FOR LOOP--")
    # now modify column k, strictly below diagonal (i = k+1:N)
    for i = k+1:N
        println("\n", "-----------------")
        println("\n", "Row (i): ", i, " Col k: ", k )

        println("\n", "BEGGINNING MATRIX")
        show(stdout, "text/plain", A)

        println("\n", "Lk[i,k] = " , -A[i,k], " / ", A[k,k])
        println("\n", "Lk_inv[i,k] = " , A[i,k], " / ", A[k,k])

        Lk[i,k] = -A[i,k] / A[k,k]    # fill me in (compute elimination factors)
        Lk_inv[i,k] = A[i,k] / A[k,k]  # fill me in (compute elimination factors)

        println("\n", "-LK New value:")
        show(stdout, "text/plain", Lk)
        println("\n", "LK_INV New value:")
        show(stdout, "text/plain", Lk_inv)
        println("\n", "RESULTING MATRIX")
        show(stdout, "text/plain", A)
        println("\n", "-----------------")
    end
 
    return (Lk, Lk_inv)

end
"""
Needs to solve Ax=b by first computing an LUP factorization
Ax = b
PAx = Pb
LUx = Pb (Solve this by using forward and backward substitution) 
"""

function LUPsolve(A, b)
    (L, U, P) = computeLUP(A)
    
    N = size(L)[1]
    y = zeros(N, 1)
    #l[3,1] y[1] + l[3,2] y[2] + l[3.3] y[3] = b3
    #y[3] = b[3] - (l[3,1] y[1] + l[3,2] y[2]) / l[3,3]
    c = copy(b)
    tmpDot = ones(N, 1)

    y[1] = c[1] / L[1, 1]

    
    #solve Ly = ᵬ for y using forward substitution - dot(c[i], L[i, i]))/
    for i = 2:N
        tmpDot = 0
        for j = 1:N
            tmpDot .= tmpDot + L[i, 1:j] * y[1:j]
        end
        y[i] = (c[i] - tmpDot) / L[i, i]
    end
    println("\n", y)
    return

end


function create_identity(N)

    I = Matrix{Float64}(undef, N, N)
    I .= 0

    for i = 1:N
        I[i, i] = 1
    end

    return I
end


A = Matrix{Float64}(undef, 3, 3)
A .= [6 -2 2;12 -8 6;3 -13 3]
b = rand(3, 1)
println("BBBBB: ", b)

#(L, U, P) = computeLUP(A)
LUPsolve(A, b)

println("\n", "\n","*-*-*-*-*FINAL VALUES*-*-*-*-*")

println("L:")
show(stdout, "text/plain", L)
println("\n", "U:")
show(stdout, "text/plain", U)
println("\n", "P:")
show(stdout, "text/plain", P)


@assert L*U ≈ P*A "P*A is not equal to L*U"