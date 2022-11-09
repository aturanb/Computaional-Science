using Plots
include("euler_methods.jl")
#=
The 1-dimensional heat equation with an initial condition is given by
    u_t = κ*u_xx + F(t, x) on the domain 0 ≤ x ≤ 1, 0 ≤ t ≤ Tf
    u(x, 0) = f(x)
where k > 0 is the thermal diffusivity. 

Assume the temperature u(x, t) is zero at the boundaries:    
    u(0, t) = u(1, t) = 0

Using second-order finite difference approx in space, and forward Euler in time. 
This becomes y' = A*y + b(t), where A is the discrete Laplacian matrix, and 
y0 = f(interior_spatial_nodes)
=#

# λ = Δt/Δx^2, Stability → λ ≤ 1/(2κ)

κ = 2
Δx = 0.1
λ = 0.5*(1/(2κ))
# take ∆t = 10^−1, 10^−2 and 10^−3 and report on how forward Euler performs
Δt = round(λ*Δx^2, digits = 10)

Tf = 2
M = Integer(Tf/Δt) # how many time steps to take

N = Integer(1/Δx) # N+1 total spatial nodes


x = collect(range(0, 1, step = Δx)) # all of the nodes
x_int = x[2:end-1] #interior spatial nodes 

t = collect(range(0, Tf, step = Δt))

#Fill in A
A = zeros(N-1, N-1)
for i = 1:N-1
    A[i, i] = -2κ*(1/Δx^2)
end

for i = 1:N-2
    A[i+1, i] = 1κ*(1/Δx^2)
    A[i,i+1] = 1κ*(1/Δx^2)
end

function F(x, t)
    return 2*(π^2 - 1)*exp(-2t)*sin(π*x) # user specified source function
end

function f(x)
    return sin(π*x) # user specifed initial heat distribution
end


function b(t)
    return F.(x_int, t)
end


# Solve y' = Ay + b(t) with y0 = f(x_int)
y0 = f.(x_int)

y = Matrix{Float64}(undef,N+1, M+1)  # my entire solution at all nodes and all time steps.

# fill in initial data:
y[:, 1] = f.(x)


# Replace below with a call to my_forward_euler_linear()
# forward Euler: y_n+1 = y_n + Δt*(A*y_n + b(t_n))
# backward Euler: y_n+1 = y_n + Δt*(A*y_n+1 + b(t_n+1))



#Forward
for n = 1:M
    # Fill in boundary condition 
    y[1, n+1] = 0
    y[N+1, n+1] = 0
    # update remaining rows
    y[2:N, n+1] = y[2:N, n] + Δt*(A*y[2:N, n] + b(t[n]))
end

#=
#Backward
for n = 1:M
    # Fill in boundary condition 
    y[1, n+1] = 0
    y[N+1, n+1] = 0
    # update remaining rows
    slope = y[2:N, n] + Δt*(y[2:N, n])
    y[2:N, n+1] = y[2:N, n] + Δt*(slope + b(t[n+1]))
end
=#
#(t, y) = my_forward_euler_linear(0, f.(x), Tf, Δt, A, b)

plot(x, y[:, 1], legend = false)
for n = 2:M+1
    p = plot(x, y[:, n], Legend = false)
    display(p)

    #sleep(1)
end
