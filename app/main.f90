program main
   use omp_lib
   use stdlib_constants, only: pi => pi_dp
   use stdlib_linalg_constants, only: ilp, dp, lk
   use stdlib_math, only: linspace
   use stdlib_io_npy, only: save_npy
   use Jacobi_Experiments
   implicit none
   integer(ilp), parameter :: nx = 512, ny = nx, maxiter = nx**2
   real(dp), dimension(nx, ny) :: b
   real(dp), dimension(nx, ny) :: u
   real(dp), dimension(nx) :: x
   real(dp), dimension(ny) :: y
   integer(ilp) :: i, j, k
   real(dp) :: start_time, end_time

   !----- Initialize variables -----!
   b = 0.0_dp; u = 0.0_dp
   x = linspace(0, 1, nx); y = linspace(0, 1, ny)

   ! Create right-hand side vector.
   do concurrent(i=1:nx, j=1:ny)
      b(i, j) = sin(2*pi*x(i))*sin(2*pi*y(j))
   end do

   call save_npy("rhs_vector.npy", b)

   !----------------------------------
   !-----     JACOBI SOLVERS     -----
   !----------------------------------

   !> Textbook Jacobi solver.
   start_time = omp_get_wtime()
   u = textbook_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   call save_npy("textbook_solution.npy", u)

   !> No-copy Jacobi solver.
   start_time = omp_get_wtime()
   u = nocopy_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   call save_npy("flipflop_solution.npy", u)

   !> On-the-fly Jacobi solver.
   start_time = omp_get_wtime()
   u = otf_norm_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   call save_npy("flipflop_solution.npy", u)

   !> doconcurrent Jacobi solver.
   start_time = omp_get_wtime()
   u = doconcurrent_solver(b, maxiter)
   end_time = omp_get_wtime()
   print *, "    - Time-to-solution     :", end_time - start_time
   call save_npy("doconcurrent_solution.npy", u)

end program main
