!--------------------------------------------------------------------------
!       Distributed parallel solution of a linear system using  GMRES     -
!                                                                         -
!                    Submitted to: Professor Jack Dongarra                -                                 
!          By: Amirehsan Ghasemi   Email: aghasemi@vols.utk.edu           -
!                                                                         -
! Bredesen Center for Interdisciplinary Research and Graduate Education   -
!            University of Tennessee, Knoxville, TN 37996, USA            -
!     Copyright@ 2024 -  Amirehsan Ghasmei. All rights reserved.          -
!--------------------------------------------------------------------------
!                  Recommended Citation & Fianl report                    -
!                   DOI: 10.13140/RG.2.2.17501.42727                      -
!--------------------------------------------------------------------------
  ! suggest compile with something like:
  ! mpif90 -Og -g -fimplicit-none -fcheck=all -fbacktrace -pedantic -fbounds-check -Wall -Wextra -Wconversion -Wunderflow laplace_par.f90
module lap
  use mpi
  implicit none

contains

  ! test initial values 
  subroutine init_u(u, rnk)
    implicit none
    real*8, dimension(:), intent(inout) :: u
    integer, intent(in) :: rnk

    ! local vars
    integer :: i, nn

    nn = size(u)

    do i = 1, nn
       u(i) = dble(rnk * nn + i)
    end do

  end subroutine init_u

  ! inits the coord system for each process  
  subroutine init_x(xmin, xmax, n, sz, rnk, x, dx)
    implicit none
    real*8, intent(in)  :: xmin, xmax
    integer, intent(in) :: n, sz, rnk
    real*8, intent(out) :: x(:), dx

    ! local vars
    integer :: i, i1, i2, jj

    dx = (xmax - xmin) / dble(n+1)

    i1 = rnk     * (n/sz)
    i2 = (rnk+1) * (n/sz) - 1

    jj = 1
    do i = i1, i2
       x(jj) = dble(i+1) * dx
       jj = jj + 1
    end do
    
  end subroutine init_x

  
  ! this subroutine synchronizes the overlaps of array [u]
  ! on all processes.
  !
  subroutine sync(u, ul, ur, rnk, sz, world, u0, uend)
    implicit none
    real*8,  intent(in), target  :: u(:)
    real*8,  intent(out) :: ul, ur
    integer, intent(in)  :: rnk, sz, world
    real*8,  intent(in)  :: u0, uend

    ! local vars
    integer :: nn, ierr, stat0(MPI_STATUS_SIZE)
    real*8, pointer :: uu

    nn = size(u)

    ! FROM LEFT TO RIGHT ...
    ! beginning node
    if (rnk .eq. 0 ) then
       ul = u0
       if ( sz > 1 ) then
          uu => u(nn)
          call MPI_SEND(uu, 1, MPI_DOUBLE, (rnk+1), 0, world, ierr)
          call MPI_RECV(ur, 1, MPI_DOUBLE, (rnk+1), 0, world, stat0, ierr)
       end if
    end if

    ! middle nodes ...
    if ( (rnk > 0 ) .and. (rnk < (sz-1)) ) then
       ! recv/send left point
       !
       call MPI_RECV(ul, 1, MPI_DOUBLE, (rnk-1), 0, world, stat0, ierr)
       uu => u(1)
       call MPI_SEND(uu, 1, MPI_DOUBLE, (rnk-1), 0, world, ierr)

       ! send/recv last point
       uu => u(nn)
       call MPI_SEND(uu, 1, MPI_DOUBLE, (rnk+1), 0, world, ierr)
       call MPI_RECV(ur, 1, MPI_DOUBLE, (rnk+1), 0, world, stat0, ierr)       
    end if

    ! finally, last node
    if (rnk .eq. (sz-1) ) then
       ur = uend
       if ( sz > 1 ) then
          call MPI_RECV(ul, 1, MPI_DOUBLE, (rnk-1), 0, world, stat0, ierr)
          uu => u(1)
          call MPI_SEND(uu, 1, MPI_DOUBLE, (rnk-1), 0, world, ierr)
       end if
    end if


  end subroutine sync

  ! computes the residuals for heat equations in the form of
  ! ddt(u) = alpha * d2/dx2(u) + f = res
  !
  subroutine comp_residuals(u, f, alpha, dx, ul, ur, res)
    implicit none
    real*8, intent(in)   :: u(:), f(:), alpha, dx, ul, ur
    real*8, intent(out)  :: res(:)

    ! local vars
    integer :: i, nn

    nn = size(u)

    ! point 1
    res(1) = alpha * (u(2) - 2.0d0 * u(1) + ul) / (dx * dx) + f(1)

    ! points 2 .. (nn-1)
    do i = 2, (nn-1)
       res(i) = alpha * (u(i+1) - 2.0d0 * u(i) + u(i-1)) / (dx * dx) + f(i)
    end do

    ! point nn
    res(nn) = alpha * (ur - 2.0d0 * u(nn) + u(nn-1)) / (dx * dx) + f(nn)

  end subroutine comp_residuals

  subroutine mpi_dot_prod(u, v, world, dot)
    implicit none
    real*8, intent(in) :: u(:), v(:)
    integer, intent(in) :: world
    real*8, intent(out) :: dot

    ! local vars
    real*8  :: tmp
    integer :: ierr
    
    tmp = sum(u*v)
    
    call MPI_ALLREDUCE(tmp, dot, 1, MPI_DOUBLE, MPI_SUM, world, ierr)

  end subroutine mpi_dot_prod
  
end module lap

program laplace_par
  use mpi
  use lap
  implicit none

  ! local vars
  real*8, dimension(:), allocatable :: x, u, f, res
  !integer :: n = 120000000  ! first
  !integer  :: n = 1200000    ! second
  integer  :: n = 12000000   ! third
  !integer :: n = 12
  integer :: nn
  integer :: sz, rnk, ierr
  real*8 :: ul, ur, u0, uend, alpha, dx, dt
  integer :: itr, max_itr
  real*8 :: tt1, tt2

  ! MPI setup
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, sz, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rnk, ierr)


  nn = n / sz
  allocate(x(nn))
  !call init_u(u, rnk)
  call init_x(0.0d0, 1.0d0, n, sz, rnk, x, dx)

  !print *, 'gg = ', n, sz, rnk, x, dx

  allocate(u(nn), f(nn), res(nn))
  f = 0.0d0
  u = 0.0d0
  res = 0.0d0
  
  !
  u0   = -2.0d0
  uend = 10.0d0
  alpha = 1.0d-15
  max_itr = 10000
  !max_itr = 100
  dt = 1.0D-8
  !dt = 1.0D-3

  ! print *, 'Hello World from process: ', rnk, 'of ', sz, 'u = ', u, 'x = ', x 
  ! start timing
  tt1 = MPI_WTIME()
  
  res = 0.0d0
  do itr = 1, max_itr
     ! sync
     call sync(u, ul, ur, rnk, sz, MPI_COMM_WORLD, u0, uend)
     ! compute residuals
     call comp_residuals(u, f, alpha, dx, ul, ur, res)
     u = u + dt * res
  end do

  tt2 = MPI_WTIME()  

  print *, '2nd Hello World from process: ', rnk, 'of ', sz, 'u = ', u, 'ul = ', ul, 'ur = ', ur
  
  if ( rnk .eq. 0 ) then
     print *, 'Official MPI Wall time = ', (tt2 - tt1), ' second(s)!'
  end if

  ! MPI exit
  call MPI_FINALIZE(ierr)

  ! clean ups
  deallocate(x, u, f, res)

end program laplace_par

