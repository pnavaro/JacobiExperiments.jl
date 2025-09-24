using DelimitedFiles
using JacobiExperiments
using CUDA

function main()

    T = Float64
    nx = 512
    ny = nx
    maxiter = nx^2
    b = zeros(T, nx, ny)
    u = similar(b)

    #----- Initialize variables -----#
    x = LinRange(0, 1, nx); y = LinRange(0, 1, ny)

    # Create right-hand side vector.
    for i in eachindex(x), j in eachindex(y)
        b[i, j] = sin(2pi * x[i]) * sin(2pi * y[j])
    end

    writedlm("rhs_vector.txt", b, " ")

    #---------------------------------#
    #-----     JACOBI SOLVERS     ----#
    #---------------------------------#

    # Textbook Jacobi solver.
    time = @elapsed begin
        u .= textbook_solver(b, maxiter)
    end
    println(" - Time-to-solution     : $time")
    writedlm("textbook_solution.txt", u, " ")

    # No-copy Jacobi solver.
    time = @elapsed begin
        u .= nocopy_solver(b, maxiter)
    end
    println(" - Time-to-solution     : $time")
    writedlm("flipflop_solution.txt", u, " ")

    # On-the-fly Jacobi solver.
    time = @elapsed begin
        u .= otf_norm_solver(b, maxiter)
    end
    println(" - Time-to-solution     : $time")
    writedlm("flipflop_solution.txt", u, " ")

    # doconcurrent Jacobi solver.
    time = @elapsed begin
        u .= doconcurrent_solver(b, maxiter)
    end
    println(" - Time-to-solution     : $time")
    writedlm("doconcurrent_solution.txt", u, " ")

    # vectorized Jacobi solver.
    time = @elapsed begin
        u .= vectorized_solver(b, maxiter)
    end
    println(" - Time-to-solution     : $time")
    writedlm("vectorized_solution.txt", u, " ")
    
    # vectorized Jacobi solver on GPU
    if CUDA.functional()
        time = @elapsed begin
            u .= gpu_solver(b, maxiter)
        end
        println(" - Time-to-solution     : $time")
        writedlm("gpu_solution.txt", u, " ")
    end

end

main()
