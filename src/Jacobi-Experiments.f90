module Jacobi_Experiments
   use stdlib_linalg_constants, only: dp, ilp, lk
   implicit none(external)
   private

   interface
      module function textbook_solver(b, maxiter) result(u)
         implicit none(external)
         real(dp), intent(in) :: b(:, :)
         integer(ilp), intent(in) :: maxiter
         real(dp), allocatable :: u(:, :)
      end function textbook_solver

      module function doconcurrent_solver(b, maxiter) result(u)
         implicit none(external)
         real(dp), intent(in) :: b(:, :)
         integer(ilp), intent(in) :: maxiter
         real(dp), allocatable :: u(:, :)
      end function doconcurrent_solver

      module function flipflop_solver(b, maxiter) result(u)
         implicit none(external)
         real(dp), intent(in) :: b(:, :)
         integer(ilp), intent(in) :: maxiter
         real(dp), allocatable :: u(:, :)
      end function flipflop_solver
   end interface

   public :: textbook_solver
   public :: doconcurrent_solver
   public :: flipflop_solver

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
      real(dp) :: dx

      ! Initialize variables
      nx = size(b, 1); ny = size(b, 2); dx = 1.0_dp/(nx - 1)
      if (nx /= ny) error stop "Number of points in each direction need to be equal."
      allocate (u(nx, ny), source=0.0_dp)
      allocate (v(nx, ny), source=0.0_dp)

      do iteration = 1, maxiter
         ! Jacobi iteration.
         call textbook_kernel(nx, ny, v, u, b, dx)
         ! Update variable.
         u = v
      end do
   end function

   module function flipflop_solver(b, maxiter) result(u)
      implicit none(external)
      real(dp), intent(in) :: b(:, :)
      integer(ilp), intent(in) :: maxiter
      real(dp), allocatable :: u(:, :)
      ! Internal variables.
      integer(ilp) :: nx, ny, i, j, iteration
      real(dp), allocatable :: v(:, :)
      real(dp) :: dx

      ! Initialize variables
      nx = size(b, 1); ny = size(b, 2); dx = 1.0_dp/(nx - 1)
      if (nx /= ny) error stop "Number of points in each direction need to be equal."
      allocate (u(nx, ny), source=0.0_dp)
      allocate (v(nx, ny), source=0.0_dp)

      do iteration = 1, maxiter, 2
         ! Jacobi kernel.
         call textbook_kernel(nx, ny, v, u, b, dx)
         ! Update variables.
         call textbook_kernel(nx, ny, u, v, b, dx)
      end do
   end function

   module function doconcurrent_solver(b, maxiter) result(u)
      implicit none(external)
      real(dp), intent(in) :: b(:, :)
      integer(ilp), intent(in) :: maxiter
      real(dp), allocatable :: u(:, :)
      ! Internal variables.
      integer(ilp) :: nx, ny, i, j, iteration
      real(dp), allocatable :: v(:, :)
      real(dp) :: dx

      ! Initialize variables
      nx = size(b, 1); ny = size(b, 2); dx = 1.0_dp/(nx - 1)
      if (nx /= ny) error stop "Number of points in each direction need to be equal."
      allocate (u(nx, ny), source=0.0_dp)
      allocate (v(nx, ny), source=0.0_dp)

      do iteration = 1, maxiter, 2
         ! Jacobi kernel.
         call doconcurrent_kernel(nx, ny, v, u, b, dx)
         ! Update variables.
         call doconcurrent_kernel(nx, ny, u, v, b, dx)
      end do
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
         u(i, j) = 0.25_dp*(b(i, j)*dx**2 - v(i + 1, j) - v(i - 1, j) &
                            - v(i, j + 1) - v(i, j - 1))/4
      end do
   end subroutine doconcurrent_kernel

end module Jacobi_Experiments
