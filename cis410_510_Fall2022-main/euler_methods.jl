using Plots

# write forward Euler for the IVP system y' = f(t, y)
# where y is a vector in R^n (i.e. has n components)


function my_forward_euler_linear(t0, y0, Tf, Δt, A, b)
    
    # y0 has N components 
    N = size(y0)

    M = Integer(Tf/Δt)  # M+1 total temporal nodes

    t = Vector{Float64}(undef, M+1)
    y = Matrix{Float64}(undef, N, M+1)

    # fill in the initial condition:
    t[1] = t0
    y[:, 1] = y0

    for n = 1:M # take M time steps
        y[:, n+1] = y[:, n] + Δt*(A*y[:, n] + b)
        t[n+1] = t[n] + Δt
    end
    
    return (t, y)
end

function my_backward_euler_linear(t0, y0, Tf, Δt, f,  λ)
    N = Integer(Tf/Δt)

    t = Vector{Float64}(undef, N+1)
    y = Vector{Float64}(undef, N+1)

    #Fill in initial condition
    t[1] = t0
    y[1] = y0

    for n=1:N
        #FIXME
        slope = y[n] + Δt*f(t[n], y[n], λ)
        y[n+1] = y[n] + Δt*f(t[n+1], slope, λ)
        t[n+1] = t[n] + Δt
    end
    return (t, y)
end


function my_forward_euler(t0, y0, Tf, Δt, f,  λ)
    N = Integer(Tf/Δt)

    t = Vector{Float64}(undef, N+1)
    y = Vector{Float64}(undef, N+1)

    #Fill in initial condition
    t[1] = t0
    y[1] = y0

    for n=1:N
        y[n+1] = y[n] + Δt*f(t[n], y[n], λ)
        t[n+1] = t[n] + Δt
    end
    return (t, y)
end

function my_backward_euler(t0, y0, Tf, Δt, f,  λ)
    N = Integer(Tf/Δt)

    t = Vector{Float64}(undef, N+1)
    y = Vector{Float64}(undef, N+1)

    #Fill in initial condition
    t[1] = t0
    y[1] = y0

    for n=1:N
        slope = y[n] + Δt*f(t[n], y[n], λ)
        y[n+1] = y[n] + Δt*f(t[n+1], slope, λ)
        t[n+1] = t[n] + Δt
    end
    return (t, y)
end

function my_test_RHS(t, y, λ)
    return λ*y
end
#=
Tf = 4
Δt = 0.1
t0 = 0
y0 = 17
λ = -3

(T, Y) = my_forward_euler(t0, y0, Tf, Δt, my_test_RHS, λ)
(T2, Y2) = my_backward_euler(t0, y0, Tf, Δt, my_test_RHS, λ)

scatter(T, Y, label = "approx Forward Euler", shape = :circle, color = :green)

scatter!(T2, Y2, label = "approx Backward Euler", shape = :circle, color = :red)

tfine = 0:Δt/20:Tf
yexact = y0*exp.(λ*tfine)

plot!(tfine, yexact, label="exact", color = :red)
=#