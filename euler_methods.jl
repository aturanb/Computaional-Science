using Plots

# write forward Euler for the IVP system y' = f(t, y)
# where y is a vector in R^n (i.e. has n components)

function my_forward_euler(t0, y0, Tf, Δt, f)
    N = length(y0)

    M = Integer(Tf/Δt)  # M+1 total temporal nodes

    t = Vector{Float64}(undef, M+1)
    y = Matrix{Float64}(undef, N, M+1)

    #Fill in initial condition
    t[1] = t0
    y[:, 1] = y0

    for n=1:M
        y[:, n+1] = y[:, n] + Δt*f(t[n],y[:, n])
        t[n+1] = t[n] + Δt
    end
    return (t, y)
end

function my_backward_euler(t0, y0, Tf, Δt, f)
    N = length(y0)

    M = Integer(Tf/Δt)  # M+1 total temporal nodes

    t = Vector{Float64}(undef, M+1)
    y = Matrix{Float64}(undef, N, M+1)


    #Fill in initial condition
    t[1] = t0
    y[:, 1] = y0

    for n=1:M
        slope = y[:, n] + Δt*f(t[n], y[:, n])
        y[:, n+1] = y[:, n] + Δt*f(t[n], slope)
        t[n+1] = t[n] + Δt
    end
    return (t, y)
end
