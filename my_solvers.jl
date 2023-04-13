using Plots 

"""
    computeLUP(A)
Compute and return PLU decomposition `LU = PA` of square matrix `A`.  
It uses partial pivoting to be able to compute all square matrices.
Calls the functions:
--> compute_Pk(Ã, k): for computing the pivoting row, and the permutation matrix.
--> compute_Lk(Ã, k): for computing the Lk and Lk_inverse.
"""
function computeLUP(A)

    N = size(A)[1]

    #Id = Matrix{Float64}(I, N, N) # N x N identity matrix
    Id = create_identity(N)

    L = copy(Id)   # initialize
    U = copy(Id)   # initialize
    P = copy(Id)   # initialize
    Ã  = copy(A) # initialize. Ã corresponds to A as it goes under elimination stages

    for k = 1:N-1 # march across columns

        (j, p) = compute_Pk(Ã, k)
        
        #update Ã
        Ã .= p * Ã
  
        #update permutation matrix
        P .= p * P

        (Lk, Lk_inv) = compute_Lk(Ã, k)

        Ã .= Lk * Ã
        L .= L * Lk_inv

        if k>1 
            #do the swap
            (L[k, 1:k-1], L[j, 1:k-1]) = (L[j, 1:k-1], L[k, 1:k-1])
        end
    end

    U .= Ã

    return (L, U, P)

end

""" 
    compute_Pk(A, k)
Finds row j, the row containing largest value in column k, including and below Ã[k,k].
Also computes the permutation matrix p[k,j].
"""
function compute_Pk(A, k)

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
    
    return (max_row, Pk)

end

"""
    compute_Lk(A, k)
Computes the Lk and its inverse from A, 
assuming first k-1 columns have undergone elimination.
"""
function compute_Lk(A, k)

    N = size(A)[1]

    Lk = create_identity(N) 
    Lk_inv = create_identity(N)

    for i = k+1:N
        Lk[i,k] = -A[i,k] / A[k,k]    
        Lk_inv[i,k] = A[i,k] / A[k,k] 
    end
 
    return (Lk, Lk_inv)

end

"""
    LUPsolve(A, b)
Function that brings multiple functions together to assert the 
correctness of the PLU decomposition by comparing A*x and b.
Calls the functions:
--> computeLUP(A): to retrieve L, U, P matrices.
--> forward_sub(P, L, b, N): Function that solves Ly = b 
    for y by using forward substitution, and returns y.
--> backward_sub(U, y, N): Function that solves Ax = y 
    for x by using backward substitution, and returns x.
"""
function LUPsolve(A, b)
    (L, U, P) = computeLUP(A)
    @assert L*U ≈ P*A
    
    N = size(L)[1]

    y = forward_sub(P, L, b, N)
    x = backward_sub(U, y, N) 

    @assert A*x ≈ b
    return x

end

"""
    forward_sub(P, L, b, N)
Solves Ly = b for y by using forward substitution, and returns y.
"""
function forward_sub(P, L, b, N)
   
    y = zeros(N, 1)
    ᵬ = copy(b)

    ᵬ .= P*ᵬ

    y[1] = ᵬ[1] / L[1, 1]
    for k = 2:N
        tmp = 0
        for i = 1:k-1
            tmp += L[k, i] * y[i]
        end
        y[k] = (ᵬ[k] - (tmp))/L[k,k]
    end
    
    return y
end

"""
    backward_sub(U, y, N): 
Solves Ax = y for x by using backward substitution, and returns x.
"""
function backward_sub(U, y, N)
    u = copy(U)
    x = zeros(N, 1)

    x[N] = y[N] / u[N, N]

    for k = N-1:-1:1
        tmp = 0
        for i = N:-1:k+1
            tmp += u[k, i] * x[i]
        end
        x[k] = (y[k] - (tmp))/u[k,k]
    end
    
    return x
end

"""
    create_identity(N)
Creates an identity matrix with a size of N*N.
"""
function create_identity(N)

    I = Matrix{Float64}(undef, N, N)
    I .= 0

    for i = 1:N
        I[i, i] = 1
    end

    return I
end


"""
    test_one()
Tests the LPU decomposition with N = 10, N = 100, N = 1000
"""
function test_one()

    x = Matrix{Float64}(undef, 3, 1)
    x .= [10;100;1000]
    y = zeros(3, 1)
    N=10
    B = rand(N, N)
    I = create_identity(N)
    A = transpose(B) * B + I
    b = rand(N, 1)
    y[1] = @elapsed LUPsolve(A, b)
    N=100
    B = rand(N, N)
    I = create_identity(N)
    A = transpose(B) * B + I
    b = rand(N, 1)
    y[2] = @elapsed LUPsolve(A, b)
    N=1000
    B = rand(N, N)
    I = create_identity(N)
    A = transpose(B) * B + I
    b = rand(N, 1)
    y[3] = @elapsed LUPsolve(A, b)
    plot(x, y, seriestype = :scatter)
end

"""
    test_two()
Tests the LPU decomposition with incrementing N by 100 in each iteration starting from 10 up to 1000.
"""
function test_two()
    N = 1000
    increment = 100
    x = 1:increment:N
    y = zeros(div(N, increment), 1)
    index = 1
    for tmpN=10:increment:N
        B = rand(tmpN, tmpN)
        I = create_identity(tmpN)
        A = transpose(B) * B + I
        b = rand(tmpN, 1)
        y[index] = @elapsed LUPsolve(A, b)
        index = index + 1
    end
    plot(x, y)
end

################ HOMEWORK 2 #################

"""
    conj_grad(A, x, b, e, max_iterations)
Performs the conjugate gradient (CG) algorithm to obtain an approximate solution x to the linear system Ax = b
A is a positive definite matrix, x is a vector of initial guess, b is a RHS vector, e is the tolerance. 
    """
function conj_grad(A, x, b, e, max_iterations)
    N = size(A)[1]
    
    r = zeros(N,max_iterations+1)
    p = zeros(N,1)
    q = zeros(N,1)
    beta = zeros(N,1)
    delta = zeros(N,1)

    rho_old = 0
    rho_new = 0

    r[:, 1] = (A * x) - b
    
    for i = 2:max_iterations
        
        rho_new = transpose(r[:, i - 1]) * r[:, i - 1]
        if i==2
            p = r[:, 1]
        else
            beta = rho_new/rho_old
            p = r[:, i-1] + beta * p
        end

        q = A * p
        delta = rho_new/(transpose(p)*q)
        x = x - (delta * p)
        r[:, i] = r[:, i-1] - (delta * q)
        
        #calculate relative error
        err = sqrt((transpose(A*x - b) * (A*x - b))/ (transpose(x) * (x)))

        #check convergence
        if err[1,1] <= e
            return x
        end
        rho_old = rho_new
    end

    return 
end

"""
    calculate_mape(cg_x, lup_x, N)
Calculates the mean absolute percentage error between:
cg_x: x vector obtained by conjugate gradient
lup_x: x vector obtained by LUPsolve
"""
function calculate_mape(cg_x, lup_x, N)
    mape = 0
    for i=1:N
        mape += abs(((lup_x[i] - cg_x[i])/lup_x[i]))
    end
    mape = mape / N
    return mape
end

"""
    test_cg()
Runs the LPU decomposition, Conjugate Gradient with incrementing N by 100 in each iteration starting from 10 up to 1000
while keeping track of Mean Absolute Percentage Error in each iteration. Then Plots these three results on the same graph.
"""
function test_cg()
    N = 1000
    e = 10e-6
    increment = 100
    
    x_axis = 1:increment:N
    time_cg = zeros(div(N, increment), 1)
    time_lup = zeros(div(N, increment), 1)
    mape = zeros(div(N, increment), 1)
    
    index = 1
    for tmpN=10:increment:N
        
        max_iterations = tmpN
        x = zeros(tmpN, 1)
        B = rand(tmpN, tmpN)
        I = create_identity(tmpN)
        A = transpose(B) * B + I
        b = rand(tmpN, 1)
        
        time_lup[index] = @elapsed (lup_x = LUPsolve(A, b))
        time_cg[index] = @elapsed (cg_x = conj_grad(A, x, b, e, max_iterations))
        mape[index] = calculate_mape(cg_x, lup_x, tmpN)
        #println(" N: ", tmpN, " | MAPE: ", mape[index])
        index = index + 1
    end
    
    plot(x_axis,[time_lup, time_cg, mape], label = ["LUPsolve" "Conjugate Gradient" "Mean Absolute Percentage Error"], title = "Comparing CG and LUP")
    xlabel!("N")
    ylabel!("Time")
    
    #=
    plot(x_axis,[time_cg, mape], label = [ "Conjugate Gradient" "Mean Absolute Percentage Error"], title = "Conjugate Gradient time graph with MAPE")
    xlabel!("N")
    ylabel!("Time")
    =#
end

test_cg()

#=
N = 100
x = zeros(N, 1)
b = rand(N,1)
B = rand(N, N)
I = create_identity(N)
A = transpose(B) * B + I
e = 10e-6
max_iterations = N

lup_x = LUPsolve(A, b)
cg_x = conj_grad(A, x, b, e, max_iterations)
mape = calculate_mape(cg_x, lup_x, N)

println("Mean absolute percentage error of x calculated by Conjugate Gradient and x calculated by LUPsolve is: ", mape)
println("Where N: ", N, " | Max iterations: ", max_iterations, " | e: ", e)
=#





