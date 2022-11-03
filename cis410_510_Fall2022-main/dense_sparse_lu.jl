using LinearAlgebra
using SparseArrays
using Plots

function fill_in(A ,N)
    for i = 1:N
        A[i, i] = -2    
    end
    
    for i=1:N-1
        A[i+1, i] = 1
        A[i, i+1] = 1
    end

    return A
end

function compute_time(A, b)
    time_lu = @elapsed F = lu(A)
    time_x = @elapsed x = F\b 
    total_time = time_lu + time_x
    return total_time
end

function test()
    
    N = 1000
    increment = 10

    dense_time = zeros(Integer(N/increment))
    sparse_time = zeros(Integer(N/increment))

    index = 1
    for i = 10:increment:N
        b = ones(i)
        dense = zeros(i,i) #non-sparse
        sparse = spzeros(i,i) #sparse
        
        dense = fill_in(dense, i)
        sparse = fill_in(sparse, i)
        
        dense_time[index] = compute_time(dense, b)
        sparse_time[index] = compute_time(sparse, b)
        
        index = index + 1
    end

    x_axis = 10:increment:N
    plot(x_axis, [dense_time, sparse_time], label = ["dense_time" "sparse_time"], title = "Comparing Sparsed and dense matrix LU")
    xlabel!("N")
    ylabel!("Time")
    
end

test()
