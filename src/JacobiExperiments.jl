module JacobiExperiments

using LinearAlgebra
using .Threads

export textbook_solver
export nocopy_solver
export otf_norm_solver
export doconcurrent_solver
export vectorized_solver

#----------------------------------
#-----     JACOBI SOLVERS     -----
#----------------------------------

function textbook_solver(b::Matrix{T}, maxiter::Int) where {T}

    # Initialize variables
    nx = size(b, 1); ny = size(b, 2); dx = 1.0 / (nx - 1)
    tol = sqrt(eps(T))
    @assert (nx == ny) "Number of points in each direction need to be equal."
    u = zeros(T, nx, ny)
    v = zeros(T, nx, ny)
    l2_norm = T(1.0)
    iteration = 0
    while ((iteration < maxiter) && (l2_norm > tol))
        # Jacobi iteration.
        textbook_kernel!(nx, ny, v, u, b, dx)
        # Compute error norm.
        l2_norm = norm(u - v)
        # Update variable.
        u .= v
        # Update iteration counter.
        iteration += 1
    end
    println("Textbook solver :")
    println("    - Number of iterations : $iteration")
    println("    - l2-norm of the error : $l2_norm")
    return u
end

function nocopy_solver(b::Matrix{T}, maxiter::Int) where {T}

    # Initialize variables
    nx = size(b, 1); ny = size(b, 2); dx = 1.0 / (nx - 1)
    tol = sqrt(eps(T))
    @assert (nx == ny) "Number of points in each direction need to be equal."
    u = zeros(T, nx, ny)
    v = zeros(T, nx, ny)
    l2_norm = T(1.0)
    iteration = 0

    while ((iteration < maxiter) && (l2_norm > tol))
        # Jacobi kernel.
        textbook_kernel!(nx, ny, v, u, b, dx)
        # Update variables.
        textbook_kernel!(nx, ny, u, v, b, dx)
        # Compute error norm.
        l2_norm = norm(u - v)
        # Update iteration counter.
        iteration += 2
    end
    println("No-copy solver :")
    println("    - Number of iterations : $iteration")
    println("    - l2-norm of the error : $l2_norm")
    return u
end

function otf_norm_solver(b::Matrix{T}, maxiter::Int) where {T}

    # Initialize variables
    nx = size(b, 1); ny = size(b, 2); dx = 1.0 / (nx - 1)
    tol = sqrt(eps(T))
    @assert (nx == ny) "Number of points in each direction need to be equal."
    u = zeros(T, nx, ny)
    v = zeros(T, nx, ny)
    l2_norm = T(1.0)
    iteration = 0

    while ((iteration < maxiter) && (l2_norm > tol))
        # Jacobi iteration.
        textbook_kernel!(nx, ny, v, u, b, dx)
        # Update variable.
        textbook_kernel!(nx, ny, u, v, b, dx)
        # Compute error norm.
        if (mod(iteration, 1000) == 0)
            l2_norm = norm(u - v)
        end
        # Update iteration counter.
        iteration += 2
    end
    l2_norm = norm(u - v)
    println("On-the-fly solver :")
    println("    - Number of iterations : $iteration")
    println("    - l2-norm of the error : $l2_norm")
    return u
end

function doconcurrent_solver(b::Matrix{T}, maxiter; ::Int) where {T}

    # Initialize variables
    nx = size(b, 1); ny = size(b, 2); dx = 1.0 / (nx - 1)
    tol = sqrt(eps(T))
    @assert (nx == ny) "Number of points in each direction need to be equal."
    u = zeros(T, nx, ny)
    v = zeros(T, nx, ny)
    l2_norm = T(1.0)
    iteration = 0

    while ((iteration < maxiter) && (l2_norm > tol))
        # Jacobi kernel.
        doconcurrent_kernel!(nx, ny, v, u, b, dx)
        # Update variables.
        doconcurrent_kernel!(nx, ny, u, v, b, dx)
        # Compute error norm.
        if (mod(iteration, 1000) == 0)
            l2_norm = zero(T)
            for i in 2:(nx - 1), j in 2:(ny - 1)
                l2_norm += (u[i, j] - v[i, j])^2
            end
            l2_norm = sqrt(l2_norm)
        end
        # Update iteration counter.
        iteration += 2
    end
    l2_norm = norm(u - v)
    println("Do-concurrent solver :")
    println("    - Number of iterations : $iteration")
    println("    - l2-norm of the error : $l2_norm")
    return u
end

function vectorized_solver(b::Matrix{T}, maxiter::Int) where {T}

    # Initialize variables
    nx = size(b, 1); ny = size(b, 2); dx = 1.0 / (nx - 1)
    tol = sqrt(eps(T))
    @assert (nx == ny) "Number of points in each direction need to be equal."
    u = zeros(T, nx, ny)
    v = zeros(T, nx, ny)
    l2_norm = T(1.0)
    iteration = 0

    while ((iteration < maxiter) && (l2_norm > tol))
        # Jacobi kernel.
        vectorized_kernel!(nx, ny, v, u, b, dx)
        # Update variables.
        vectorized_kernel!(nx, ny, u, v, b, dx)
        # Compute error norm.
        if (mod(iteration, 1000) == 0)
            l2_norm = zero(T)
            for i in 2:(nx - 1), j in 2:(ny - 1)
                l2_norm += (u[i, j] - v[i, j])^2
            end
            l2_norm = sqrt(l2_norm)
        end
        # Update iteration counter.
        iteration += 2
    end
    l2_norm = norm(u - v)
    println("Vectorized solver :")
    println("    - Number of iterations : $iteration")
    println("    - l2-norm of the error : $l2_norm")
    return u
end

#----------------------------------
#-----     JACOBI KERNELS     -----
#----------------------------------

function textbook_kernel!(nx::Int, ny::Int, u::Matrix{T}, v::Matrix{T}, b::Matrix{T}, dx::T) where {T}
    for j in 2:(ny - 1)
        for i in 2:(nx - 1)
            u[i, j] = T(0.25) * (
                b[i, j] * dx^2 - v[i + 1, j] - v[i - 1, j]
                    - v[i, j + 1] - v[i, j - 1]
            )
        end
    end
end

function doconcurrent_kernel!(nx::Int, ny::Int, u::Matrix{T}, v::Matrix{T}, b::Matrix{T}, dx::T) where {T}
    @threads for j in 2:(ny - 1)
        for i in 2:(nx - 1)
            u[i, j] = T(0.25) * (
                b[i, j] * dx^2 - v[i + 1, j] - v[i - 1, j]
                    - v[i, j + 1] - v[i, j - 1]
            )
        end
    end
end

function vectorized_kernel!(nx::Int, ny::Int, u::Matrix{T}, v::Matrix{T}, b::Matrix{T}, dx::T) where {T}

    i = 2:(nx - 1)
    j = 2:(ny - 1)
    @views u[i, j] .= T(0.25) .* ( b[i, j] .* dx^2 - v[i .+ 1, j] .- v[i .- 1, j]
            .- v[i, j .+ 1] .- v[i, j .- 1]
    )
end

end
