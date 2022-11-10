using Plots
include("euler_methods.jl")
#=
The 1-dimensional heat equation with an initial condition is given by
    u_t = κ*u_xx + F(t, x) on the domain 0 ≤ x ≤ 1, 0 ≤ t ≤ Tf
    u(x, 0) = f(x)
where k > 0 is the thermal diffusivity. 

Assume the temperature u(x, t) is zero at the boundaries:    
    u(0, t) = u(1, t) = 0
=#

κ = 2
Δx = 0.1
λ = 0.5*(1/(2κ)) # λ = Δt/Δx^2, Stability → λ ≤ 1/(2κ)
Δt = round(λ*Δx^2, digits = 10) # take ∆t = 10^−1, 10^−2 and 10^−3 and report on how forward Euler performs
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

function exact(t, x)
    return exp(-2t) * sin(π*x)
end

function my_RHS(t, y)
    return A*y + b(t)
end

# Solve y' = Ay + b(t) with y0 = f(x_int)
y0 = f.(x_int)

feY = zeros(N+1, M+1)  # my entire solution at all nodes and all time steps for forward euler.
beY = zeros(N+1, M+1)  # my entire solution at all nodes and all time steps for backward euler.

# fill in initial data:
feY[:, 1] = f.(x)
beY[:, 1] = f.(x)

# call my_backward_euler and my_forward_euler
(t, tmpy) = my_forward_euler(0, y0, Tf, Δt, my_RHS)
feY[2:N, :] = tmpy
(t, tmpy) = my_backward_euler(0, y0, Tf, Δt, my_RHS)
beY[2:N, :] = tmpy

xfine = collect(range(0, 1, step = Δx/100))
#=
for n = 1:M+1
    p = plot(x, feY[:, n], label = ["numerical_forward"])
    plot!(x, beY[:, n], label = ["numerical_backward"])
    plot!(xfine, exact.(t[n], xfine), label = ["exact"])
    ylims!((0, 1))
    display(p)
    
    sleep(0.1)
end
=#

