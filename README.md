# heat_MOL.jl

This code provides the 1-dimensional heat equation with an initial condition using the forward and backward Euler methods. The heat equation describes the diffusion of heat in a medium over time, and it is given by:

u_t = κ*u_xx + F(t, x)

where κ is the thermal diffusivity, u is the temperature distribution, and F is a user-specified source function. The heat equation is subject to the initial condition:

u(x, 0) = f(x)

where f is the user-specified initial heat distribution. Additionally, the temperature is zero at the boundaries:

u(0, t) = u(1, t) = 0.

The code defines a grid of spatial nodes x and time steps t, and it fills in the matrix A that represents the second derivative operator using finite differences. It also defines the source function F, the initial heat distribution f, and the right-hand side of the ODE y' = Ay + b(t) with y0 = f(x_int), where x_int is the vector of interior spatial nodes.

The code then solves the ODE using the my_forward_euler and my_backward_euler functions, which implement the forward and backward Euler methods, respectively. It fills in the solution matrices feY and beY with the numerical solutions obtained by the forward and backward Euler methods, respectively. Finally, it plots the numerical solutions and the exact solution using the Plots package.

Overall, the code provides a simple example of how to solve the heat equation using finite differences and the forward and backward Euler methods, and it illustrates the importance of choosing an appropriate time step to ensure numerical stability.

```bash
julia

include("heat_MOL.jl")
```

## my_solvers.jl

This code is an implementation of the PLU decomposition algorithm with partial pivoting, which is a method for solving systems of linear equations.

The computeLUP(A) function computes the PLU decomposition of a given square matrix A using partial pivoting. It uses the compute_Pk(A, k) and compute_Lk(A, k) functions to perform the necessary operations. The PLU decomposition is returned as a tuple (L, U, P).

The compute_Pk(A, k) function finds the row j containing the largest value in column k of A, including and below A[k,k]. It also computes the permutation matrix p[k,j].

The compute_Lk(A, k) function computes the Lk and its inverse from A, assuming that the first k-1 columns have undergone elimination.

The LUPsolve(A, b) function brings together multiple functions to assert the correctness of the PLU decomposition by comparing A*x and b. It calls the computeLUP(A) function to retrieve the L, U, and P matrices, then uses the forward_sub(P, L, b, N) and backward_sub(U, y, N) functions to solve the system of linear equations.

The forward_sub(P, L, b, N) function solves Ly = b for y by using forward substitution.

The backward_sub(U, y, N) function solves Ax = y for x by using backward substitution.

The create_identity(N) function creates an identity matrix with a size of N*N.

```bash
julia

include("my_solvers.jl")
```
