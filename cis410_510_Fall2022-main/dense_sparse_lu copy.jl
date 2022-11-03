using LinearAlgebra
using SparseArrays

N = 1000
K = zeros(N,N)
J = spzeros(N,N)
b = ones(N)

for i = 1:N
    K[i, i] = -2
    
end

for i=1:N-1
    K[i+1, i] = 1
    K[i, i+1] = 1
end

time_non_spars_lu = @elapsed non_spars = lu(K)
time_non_spars_x = @elapsed x = non_spars\b 
time_non_spars_total = time_non_spars_lu + time_non_spars_x

println("Total Time for a non-sparse :", time_non_spars_total)


for i = 1:N
    J[i, i] = -2
    
end

for i=1:N-1
    J[i+1, i] = 1
    J[i, i+1] = 1
end


time_spars_lu = @elapsed spars = lu(J)
time_spars_x = @elapsed x = spars\b 
time_total_spars = time_spars_lu + time_spars_x

println("Total Time for a sparse :", time_total_spars)
println("That is almost ", Integer(round(time_non_spars_total/time_total_spars)), " times faster!")