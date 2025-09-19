module Jacobi_Experiments
   use stdlib_linalg_constants, only: dp, ilp, lk
   implicit none(external)
   private

   real(dp), parameter :: eps = epsilon(1.0_dp)
   real(dp), parameter :: tol = sqrt(eps)

   interface
      module function textbook_solver(b, maxiter) result(u)
         implicit none(external)
         real(dp), intent(in) :: b(:, :)
         integer(ilp), intent(in) :: maxiter
         real(dp), allocatable :: u(:, :)
      end function textbook_solver

      module function nocopy_solver(b, maxiter) result(u)
         implicit none(external)
         real(dp), intent(in) :: b(:, :)
         integer(ilp), intent(in) :: maxiter
         real(dp), allocatable :: u(:, :)
      end function nocopy_solver

      module function otf_norm_solver(b, maxiter) result(u)
         implicit none(external)
         real(dp), intent(in) :: b(:, :)
         integer(ilp), intent(in) :: maxiter
         real(dp), allocatable :: u(:, :)
      end function otf_norm_solver

      module function doconcurrent_solver(b, maxiter) result(u)
         implicit none(external)
         real(dp), intent(in) :: b(:, :)
         integer(ilp), intent(in) :: maxiter
         real(dp), allocatable :: u(:, :)
      end function doconcurrent_solver
   end interface

   public :: textbook_solver
   public :: nocopy_solver
   public :: otf_norm_solver
   public :: doconcurrent_solver

contains

   !----------------------------------
   !-----     JACOBI SOLVERS     -----
   !----------------------------------

   module function textbook_solver(b, maxiter) result(u)
      implicit none(external)
      real(dp), intent(in) :: b(:, :)
      integer(ilp), intent(in) :: maxiter
      real(dp), allocatable :: u(:, :)
      ! Internal variables.
      integer(ilp) :: nx, ny, i, j, iteration
      real(dp), allocatable :: v(:, :)
      real(dp) :: dx, l2_norm

      ! Initialize variables
      nx = size(b, 1); ny = size(b, 2); dx = 1.0_dp/(nx - 1)
      if (nx /= ny) error stop "Number of points in each direction need to be equal."
      allocate (u(nx, ny), source=0.0_dp)
      allocate (v(nx, ny), source=0.0_dp)
      l2_norm = 1.0_dp
      iteration = 0
      do while ((iteration < maxiter) .and. (l2_norm > tol))
         ! Jacobi iteration.
         call textbook_kernel(nx, ny, v, u, b, dx)
         ! Compute error norm.
         l2_norm = norm2(u - v)
         ! Update variable.
         u = v
         ! Update iteration counter.
         iteration = iteration + 1
      end do
      print *, "Textbook solver :"
      print *, "    - Number of iterations :", iteration
      print *, "    - l2-norm of the error :", l2_norm
   end function

   module function nocopy_solver(b, maxiter) result(u)
      implicit none(external)
      real(dp), intent(in) :: b(:, :)
      integer(ilp), intent(in) :: maxiter
      real(dp), allocatable :: u(:, :)
      ! Internal variables.
      integer(ilp) :: nx, ny, i, j, iteration
      real(dp), allocatable :: v(:, :)
      real(dp) :: dx, dummy(2), l2_norm

      ! Initialize variables
      nx = size(b, 1); ny = size(b, 2); dx = 1.0_dp/(nx - 1)
      if (nx /= ny) error stop "Number of points in each direction need to be equal."
      allocate (u(nx, ny), source=0.0_dp)
      allocate (v(nx, ny), source=0.0_dp)
      l2_norm = 1.0_dp
      iteration = 0

      do while ((iteration < maxiter) .and. (l2_norm > tol))
         ! Jacobi kernel.
         call textbook_kernel(nx, ny, v, u, b, dx)
         ! Update variables.
         call textbook_kernel(nx, ny, u, v, b, dx)
         ! Compute error norm.
         l2_norm = norm2(u - v)
         ! Update iteration counter.
         iteration = iteration + 2
      end do
      print *, "No-copy solver :"
      print *, "    - Number of iterations :", iteration
      print *, "    - l2-norm of the error :", l2_norm
   end function

   module function otf_norm_solver(b, maxiter) result(u)
      implicit none(external)
      real(dp), intent(in) :: b(:, :)
      integer(ilp), intent(in) :: maxiter
      real(dp), allocatable :: u(:, :)
      ! Internal variables.
      integer(ilp) :: nx, ny, i, j, iteration
      real(dp), allocatable :: v(:, :)
      real(dp) :: dx, l2_norm

      ! Initialize variables
      nx = size(b, 1); ny = size(b, 2); dx = 1.0_dp/(nx - 1)
      if (nx /= ny) error stop "Number of points in each direction need to be equal."
      allocate (u(nx, ny), source=0.0_dp)
      allocate (v(nx, ny), source=0.0_dp)
      l2_norm = 1.0_dp
      iteration = 0

      do while ((iteration < maxiter) .and. (l2_norm > tol))
         ! Jacobi iteration.
         call textbook_kernel(nx, ny, v, u, b, dx)
         ! Update variable.
         call textbook_kernel(nx, ny, u, v, b, dx)
         ! Compute error norm.
         if (mod(iteration, 1000) == 0) then
            l2_norm = norm2(u - v)
         end if
         ! Update iteration counter.
         iteration = iteration + 2
      end do
      l2_norm = norm2(u - v)
      print *, "On-the-fly solver :"
      print *, "    - Number of iterations :", iteration
      print *, "    - l2-norm of the error :", l2_norm
   end function

   module function doconcurrent_solver(b, maxiter) result(u)
      implicit none(external)
      real(dp), intent(in) :: b(:, :)
      integer(ilp), intent(in) :: maxiter
      real(dp), allocatable :: u(:, :)
      ! Internal variables.
      integer(ilp) :: nx, ny, i, j, iteration
      real(dp), allocatable :: v(:, :)
      real(dp) :: dx, l2_norm

      ! Initialize variables
      nx = size(b, 1); ny = size(b, 2); dx = 1.0_dp/(nx - 1)
      if (nx /= ny) error stop "Number of points in each direction need to be equal."
      allocate (u(nx, ny), source=0.0_dp)
      allocate (v(nx, ny), source=0.0_dp)
      l2_norm = 1.0_dp
      iteration = 0

      do while ((iteration < maxiter) .and. (l2_norm > tol))
         ! Jacobi kernel.
         call doconcurrent_kernel(nx, ny, v, u, b, dx)
         ! Update variables.
         call doconcurrent_kernel(nx, ny, u, v, b, dx)
         ! Compute error norm.
         if (mod(iteration, 1000) == 0) then
            l2_norm = 0.0_dp
            do concurrent(i=2:nx - 1, j=2:ny - 1) reduce(+:l2_norm)
               l2_norm = l2_norm + (u(i, j) - v(i, j))**2
            end do
            l2_norm = sqrt(l2_norm)
         end if
         ! Update iteration counter.
         iteration = iteration + 2
      end do
      l2_norm = norm2(u - v)
      print *, "Do-concurrent solver :"
      print *, "    - Number of iterations :", iteration
      print *, "    - l2-norm of the error :", l2_norm
   end function

   !----------------------------------
   !-----     JACOBI KERNELS     -----
   !----------------------------------

   pure subroutine textbook_kernel(nx, ny, u, v, b, dx)
      implicit none(external)
      integer(ilp), intent(in) :: nx, ny
      real(dp), intent(out) :: u(nx, ny)
      real(dp), intent(in) :: v(nx, ny), b(nx, ny), dx
      integer(ilp) :: i, j
      do j = 2, ny - 1
         do i = 2, nx - 1
            u(i, j) = 0.25_dp*(b(i, j)*dx**2 - v(i + 1, j) - v(i - 1, j) &
                               - v(i, j + 1) - v(i, j - 1))
         end do
      end do
   end subroutine textbook_kernel

   pure subroutine doconcurrent_kernel(nx, ny, u, v, b, dx)
      implicit none(external)
      integer(ilp), intent(in) :: nx, ny
      real(dp), intent(out) :: u(nx, ny)
      real(dp), intent(in) :: v(nx, ny), b(nx, ny), dx
      integer(ilp) :: i, j
      do concurrent(i=2:nx - 1, j=2:ny - 1)
         ! Jacobi update.
         u(i, j) = 0.25_dp*(b(i, j)*dx**2 - v(i + 1, j) - v(i - 1, j) &
                            - v(i, j + 1) - v(i, j - 1))
      end do
   end subroutine doconcurrent_kernel

end module Jacobi_Experiments
