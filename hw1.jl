function upper_triangulate(A, b)

    N = size(A)[1]
    M = [A b]

    for j = 1:N-1

        eliminate!(M, j)
        show(stdout, "text/plain", M)
        sleep(5)
    end

    return M
end

function eliminate!(M, j)

    pivot = M[j, j]
    pivot_row = M[j, :]

    for k = j+1:N
        fac = M[k, j]/pivot
        M[k, :] = M[k, :] - fac*pivot_row
    end

end



A = Matrix{Float64}(undef, 4, 4)
A[1, :] = [1 2 -1 1]
A[2, :] = [-1 1 2 -1]
A[3, :] = [2 -1 2 2]
A[4, :] = [1 1 -1 2]

b = [6; 3; 14; 8.0]

# Turn the augmented system [A | b] into upper triangular form. 
# March across columns 1:N-1

output = upper_triangulate(A, b)


#=
for j = 1:N-1
    # Now I'm on column j
    # Need to eliminate all entries below a_jj using row j
    # get pivot element
    pivot = M[j, j]
    pivot_row = M[j, :]
    for k = j+1:N
        fac = M[k, j]/pivot
        M[k, :] = M[k, :] - fac*pivot_row
    end
    show(stdout, "text/plain", M)
    sleep(5)
end
=#
