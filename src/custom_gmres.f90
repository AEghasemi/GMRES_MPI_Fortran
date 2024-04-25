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
!              DOI: https://doi.org/10.13140/RG.2.2.17501.42727           -
!--------------------------------------------------------------------------
module custom_gmres
  use mpi
  implicit none

  ! type used only inside gmres
  ! for performance increase, we just allocate this once
  type gmres_dat

     real*8, dimension(:,:), allocatable :: v, h
     real*8, dimension(:), allocatable :: r0, g, yk, uk, s, c, zk, alpha, tmp

  end type gmres_dat

  ! type used for collecting solution related
  ! data and structures
  type sol_dat
     
     ! real*8, allocatable :: A(:,:)
     real*8 :: ubeg, uend
     real*8 :: ul, ur, ep, alpha, dx, dt
     real*8, allocatable :: u(:), u0(:), Fs(:), Fsp(:)
     
  end type sol_dat

  type mpi_dat
     ! MPI vars
     integer :: sz, rnk, ierr
  end type mpi_dat

contains

  subroutine init_gmres_dat(gd, num, nrows)
    implicit none
    type(gmres_dat), intent(inout) :: gd
    integer, intent(in) :: num, nrows

    ! local vars
    integer :: nvars

    nvars = nrows

    ! we don't use dynamic mem allocation in the gmres inner-outer loops 
    ! to increase speed. this might not be memory efficient though.
    allocate(gd%v(nvars,num+1), gd%r0(nrows), gd%g(num+1), gd%s(num), gd%c(num))
    allocate(gd%yk(nrows), gd%uk(nvars), gd%zk(nvars) )
    allocate(gd%h(num+1,num), gd%alpha(num))

  end subroutine init_gmres_dat

  subroutine free_gmres_dat(gd)
    implicit none
    type(gmres_dat), intent(inout) :: gd


    deallocate(gd%v, gd%r0, gd%g, gd%s, gd%c)
    deallocate(gd%yk, gd%uk, gd%zk)
    deallocate(gd%h, gd%alpha)

  end subroutine free_gmres_dat

  subroutine init_sol_dat(dat, n, ul, ur, ep, alpha, dx, dt)
    implicit none
    type(sol_dat), intent(inout) :: dat
    integer, intent(in)          :: n
    real*8, intent(in)           :: ul, ur, ep, alpha, dx, dt

    dat%ubeg   = ul
    dat%uend = ur
    dat%ep = ep
    dat%alpha = alpha
    dat%dx = dx
    dat%dt = dt

    allocate(dat%u(n), dat%u0(n), dat%Fs(n), dat%Fsp(n))

    dat%u  = 0.0d0
    dat%u0 = 0.0d0
    dat%Fs = 0.0d0
    dat%Fsp= 0.0d0

  end subroutine init_sol_dat

  subroutine free_sol_dat(dat)
    implicit none
    type(sol_dat), intent(inout) :: dat

    deallocate(dat%u, dat%u0, dat%Fs, dat%Fsp)


  end subroutine free_sol_dat

  subroutine comp_Fs(mp, dat, u, u0, alpha, dx, dt, Fs)
    implicit none
    type(mpi_dat), intent(in) :: mp
    type(sol_dat) :: dat    
    real*8, intent(in)  :: u(:), u0(:), alpha, dx, dt
    real*8, intent(out) :: Fs(:) 

    ! local vars
    integer :: nn, i
    !real*8 :: ul, ur

    ! sync
    call sync(u, dat%ul, dat%ur, mp%rnk, mp%sz, MPI_COMM_WORLD, dat%ubeg, dat%uend)
    
    nn = size(u)

    Fs(1) = -(u(1) - u0(1))/dt + alpha / (dx*dx) * (u(2) - 2.0d0*u(1) + dat%ul)

    do i = 2, (nn-1)
       Fs(i) = -(u(i) - u0(i))/dt + alpha / (dx*dx) * (u(i+1) - 2.0d0*u(i) + u(i-1))
    end do

    Fs(nn) = -(u(nn) - u0(nn))/dt + alpha / (dx*dx) * (dat%ur - 2.0d0*u(nn) + u(nn-1))

  end subroutine comp_Fs
  
  subroutine comp_Ax(mp, dat, x, Ax)
    implicit none
    type(mpi_dat)       :: mp
    type(sol_dat), intent(inout) :: dat
    real*8, intent(in)           :: x(:)
    real*8, intent(out)          :: Ax(:)

    ! can be over-loaded
    !Ax = matmul(dat%A, x)

    call comp_Fs(mp, dat, dat%u, dat%u0, dat%alpha, dat%dx, dat%dt, dat%Fs)
    call comp_Fs(mp, dat, (dat%u + dat%ep * x), dat%u0, dat%alpha, dat%dx, dat%dt, dat%Fsp)

    Ax = (1.0d0/dat%ep) * (dat%Fsp - dat%Fs)
    
  end subroutine comp_Ax
    
  subroutine gmres_orig(mp, dat, gd, b, epsil, num, nrst, sol)
    implicit none
    ! inputs
    type(mpi_dat) :: mp
    type(sol_dat), intent(inout) :: dat
    type(gmres_dat), intent(inout) :: gd
    !real*8, dimension(:,:), intent(in) :: A
    real*8, dimension(:), intent(in) :: b
    real*8, intent(in) :: epsil
    integer, intent(in) :: num, nrst

    !outputs
    real*8, dimension(:), intent(inout) :: sol

    ! locals
    integer :: k, j, i
    integer :: finish_flag, jrst
    real*8 :: norm_tmp, delta, gamma, resi
    integer :: nrows, nvars, len_res


    ! resetting the flag
    finish_flag = 0
    resi = 0.0d0
    nrows = size(b)
    nvars = size(b)
    len_res = 0


    ! do j number of restarts
    do jrst = 1, nrst    
       ! reset ------------
       gd%v = 0.0d0; gd%h = 0.0d0
       gd%c = 0.0d0; gd%s = 0.0d0
       gd%zk = 0.0d0; gd%g = 0.0d0
       gd%alpha = 0.0d0
       gd%yk = 0.0d0
       ! ------------------
       ! ***********************
       ! call mf%Ax(sol, r0)
       ! r0 = b - r0
       call comp_Ax(mp, dat, sol, gd%r0)
       !gd%r0 = b - matmul(A,sol)
       gd%r0 = b - gd%r0 
       ! ***********************
       gd%v(:,1) = gd%r0
       call mpi_dot_prod(u = gd%v(:,1), v = gd%v(:,1), world = MPI_COMM_WORLD, dot = norm_tmp)
       norm_tmp = sqrt(norm_tmp)
       !norm_tmp = norm(gd%v(:,1))
       gd%v(:,1) = gd%v(:,1)/norm_tmp
       gd%g(1) = norm_tmp
       ! Arnoldi iterations    
       do k = 1, num
          gd%g(k+1) = 0.0d0
          ! ***********************
          gd%yk = 0.0d0
          gd%yk = gd%v(:,k)
          ! call mf%Ax(yk, r0)
          ! uk = r0(mf%free)
          !gd%yk = gd%v(:,k)
          !gd%uk = matmul(A, gd%yk)
          call comp_Ax(mp, dat, gd%yk, gd%uk)
          ! ***********************
          do j = 1, k
             call mpi_dot_prod(u = gd%v(:,j), v = gd%uk, world = MPI_COMM_WORLD, dot = gd%h(j,k))             
             !gd%h(j,k) = dot_product(gd%v(:,j),gd%uk)
             !grahm-schmidt orthogonalization
             gd%uk = gd%uk - gd%h(j,k) * gd%v(:,j) 
          end do
          !
          call mpi_dot_prod(u = gd%uk, v = gd%uk, world = MPI_COMM_WORLD, dot = gd%h(k+1,k))
          gd%h(k+1,k) = sqrt(gd%h(k+1,k))
          !gd%h(k+1,k) = norm(gd%uk)        
          gd%v(:,k+1) = gd%uk/gd%h(k+1,k)
          ! applying givens
          do j = 1, (k-1)
             delta = gd%h(j,k)
             gd%h(j,k) = gd%c(j)*delta + gd%s(j)*gd%h(j+1,k)
             gd%h(j+1,k) = -gd%s(j)*delta + gd%c(j) * gd%h(j+1,k)
          end do
          gamma = sqrt(gd%h(k,k)**(2.0d0) + gd%h(k+1,k)**(2.0d0))
          gd%c(k) = gd%h(k,k) / gamma
          gd%s(k) = gd%h(k+1,k) / gamma
          gd%h(k,k) = gamma
          gd%h(k+1,k) = 0.0d0
          delta = gd%g(k)
          gd%g(k) = gd%c(k) * delta + gd%s(k) * gd%g(k+1)
          gd%g(k+1) = -gd%s(k) * delta + gd%c(k) * gd%g(k+1)
          resi = abs(gd%g(k+1))
          ! add to residual history
          len_res = len_res + 1
          ! allocate(tmp(len_res))
          ! if(len_res > 1) tmp(1:(len_res-1)) = res(1:(len_res-1)) 
          ! tmp(len_res) = resi
          ! call move_alloc(tmp, res)
          !
          if( resi <= epsil) then
             finish_flag = 1
             goto 100
          end if

       end do
       k = k - 1
       ! solving backward for alpha
100    gd%alpha(k) = gd%g(k) / gd%h(k,k)
       do i = k-1, 1, -1
          gd%alpha(i) = gd%g(i) / gd%h(i,i)
          do j = i+1, k      
             gd%alpha(i) = gd%alpha(i) - gd%h(i,j)/gd%h(i,i) * gd%alpha(j)     
          end do
       end do
       ! compute the final directional vector using
       ! the combination of search vectors 
       do j = 1, k
          gd%zk = gd%zk + gd%alpha(j)*gd%v(:,j)
       end do
       ! updating solution
       ! ***********************
       ! zk = zk/N
       ! r0 = 0.0d0
       gd%r0 = gd%zk
       sol = sol + gd%r0     
       ! sol = sol + zk/N    
       ! ***********************
       if(finish_flag == 1) then 
          goto 200        
       end if
    end do
    ! subroutine final clean-ups
200 return !done <here>

  end subroutine gmres_orig

  ! function norm(x)
  !   implicit none
  !   real*8, dimension(:), intent(in) :: x
  !   real*8 :: norm

  !   ! locals
  !   integer :: i

  !   norm = 0.0d0 !reset
  !   do i = 1, size(x)
  !      norm = norm + x(i)* x(i)
  !   end do

  !   norm = sqrt(norm)

  !   ! done
  ! end function norm

  subroutine heat_eq_matrix(xmin, xmax, umin, umax, f, A, b)
    implicit none
    real*8, intent(in) :: xmin, xmax, umin, umax
    real*8, dimension(:), intent(in) :: f
    real*8, dimension(:,:), intent(out) :: A
    real*8, dimension(:), intent(out) :: b
    

    ! local vars
    integer :: n, i, j
    real*8  :: dx

    n = size(A, 1)    
    dx = (xmax - xmin)/ dble(n-1)
    
    !hard reset
    A = 0.0d0
    b = 0.0d0
    
    ! at point 1
    A(1,1) = -2.0d0
    A(1,2) = 1.0d0
    b(1)   = f(1) * dx * dx - umin 
    
    ! between 
    do i = 2, (n-1)
       b(i) = f(i) * dx * dx
       do j = 2, (n-1)
          if ( i .eq. j ) then
             A(i,j)   = -2.0d0
             A(i,j-1) = 1.0d0
             A(i,j+1) = 1.0d0
          end if
       end do
    end do

    ! at point n
    A(n,n) = -2.0d0
    A(n,n-1) = 1.0d0
    b(n)   = f(n) * dx * dx - umax 

  end subroutine heat_eq_matrix

  ! this subroutine synchronizes the overlaps of array [u]
  ! on all processes.
  !
  subroutine sync(u, ul, ur, rnk, sz, world, ubeg, uend)
    implicit none
    real*8,  intent(in), target  :: u(:)
    real*8,  intent(out) :: ul, ur
    integer, intent(in)  :: rnk, sz, world
    real*8,  intent(in)  :: ubeg, uend

    ! local vars
    integer :: nn, ierr, stat0(MPI_STATUS_SIZE)
    real*8, pointer :: uu

    nn = size(u)

    ! FROM LEFT TO RIGHT ...
    ! beginning node
    if (rnk .eq. 0 ) then
       ul = ubeg
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
  
end module custom_gmres

program tester
  use custom_gmres
  implicit none

  ! gmres specific types 
  type(gmres_dat) :: gd
  type(sol_dat) :: dat

  integer :: n = 120
  integer :: nn, ii, s, num, nrst
  real*8, allocatable :: b(:), x0(:), sol(:)
  real*8 :: epsil
  type(mpi_dat) :: mp

  
  ! MPI setup
  call MPI_INIT(mp%ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mp%sz, mp%ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mp%rnk, mp%ierr)

  
  nn = n/mp%sz
  num = 10
  allocate(b(nn), x0(nn), sol(nn))
  call init_sol_dat(dat, nn, -1.0d0, 10.0d0, 1.0d-11, 1.0d0, (1.0d0/dble(nn+1)), 1.0d-3)
  call init_gmres_dat(gd, num, nn)



  do ii = 1, 10000
     !newton iterations
     do s = 1, 3

        !compute Fs
        call comp_Fs(mp, dat, dat%u, dat%u0, dat%alpha, dat%dx, dat%dt, dat%Fs)

        !Compute Du = -A\b

        ! gmres specific constants
        b = -dat%Fs
        epsil = 1.0D-14
        x0 = dat%u
        nrst = 5

        ! call GMRES
        sol = dat%u
        ! subroutine gmres_orig(dat, gd, b, epsil, num, nrst, sol)    
        call gmres_orig(mp, dat, gd, b, epsil, num, nrst, sol)


        ! update
        dat%u = dat%u + sol 

     end do
     dat%u0 = dat%u
  end do

  print *, 'rank = ', mp%rnk, 'u = ', dat%u

  ! MPI exit
  call MPI_FINALIZE(mp%ierr)
  
  ! clean ups
  call free_gmres_dat(gd)
  call free_sol_dat(dat)
  deallocate( b, x0, sol)
  

  
end program tester

! program tester
!   use custom_gmres

!   integer :: n, i, j
!   real*8, dimension(:), allocatable :: bb, ff, sol 

!   ! gmres specific types 
!   type(gmres_dat) :: gd
!   type(sol_dat) :: dat

!   n = 100

!   call init_sol_dat(dat, n)
  
!   allocate(bb(n), ff(n), sol(n))

!   ! no heat source (you can change it to include heat source)
!   ff = 0.0d0
!   call heat_eq_matrix(xmin = 0.0d0, xmax = 1.0d0, umin = 1.0d0, umax = 2.0d0, f = ff, A = dat%A, b = bb)
  

!   do i = 1, n 
!      print *, dat%A(i,:)
!   end do

!   do i = 1, n 
!      print *, bb(i)
!   end do

!   call init_gmres_dat(gd, n, n)
    
!   sol = 0.0d0
!   call gmres_orig(dat, gd, bb, dble(1.0e-12), n, 1, sol)

!   print *, sol

!   ! clean ups
!   call free_gmres_dat(gd)
!   call free_sol_dat(dat)
!   deallocate( bb, ff, sol)
  
  
! end program tester


  

! program tester
!   use custom_gmres
!   implicit none

!   real*8, dimension(4,4) :: A
!   real*8, dimension(4)   :: b, sol

  
!   A(1, :) = (/ 1.0d0 , 4.0d0, 2.0d0 , 0.0d0 /)
!   A(2, :) = (/ -1.0d0, 4.0d0, 0.0d0 , 0.0d0 /)
!   A(3, :) = (/ 0.0d0, -4.0d0, 2.0d0 , 3.0d0 /)
!   A(4, :) = (/ 0.0d0, 0.0d0, 2.0d0 , 5.0d0 /)

!   b = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0 /)

!   sol = 0.0d0
!   call gmres_orig(A, b, dble(1.0e-12), 4, 1, sol)

!   print *, sol

! end program tester

