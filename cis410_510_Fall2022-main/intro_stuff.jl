# this is julia code

function add_stuff(x, y)
    w = 5
    return x + y
end


x = 5
y = 3

sum = add_stuff(x, y)
7
# add a comment

function traverse(A)
    N = size(A)[1]
    k = 3
    for row = k:N
        println("sa")
        println(A[row, k])
    end
end

A = Matrix{Float64}(undef, 3, 3)
A .= [6 -2 2;12 -8 6;3 -13 3]
show(stdout, "text/plain", A)
traverse(A)