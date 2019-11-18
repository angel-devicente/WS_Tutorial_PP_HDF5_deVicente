! Main solver routines for heat equation solver
module core
  use heat

contains

  ! Exchange the boundary data between MPI tasks
  subroutine exchange(field0, parallel)
    use mpi
    implicit none

    type(field), intent(inout) :: field0
    type(parallel_data), intent(inout) :: parallel

    integer :: ierr
    integer, dimension(MPI_STATUS_SIZE) :: status

    ! Send to left, receive from right
    call mpi_sendrecv(field0%data(:,1),field0%nx,MPI_DOUBLE_PRECISION,parallel%nleft,11, &
         field0%data(:,field0%ny+1),field0%nx,MPI_DOUBLE_PRECISION,parallel%nright,11,&
         parallel%comm,status,ierr)
    
    ! Send to right, receive from left
    call mpi_sendrecv(field0%data(:,field0%ny),field0%nx,MPI_DOUBLE_PRECISION,parallel%nright,12, &
         field0%data(:,0),field0%nx,MPI_DOUBLE_PRECISION,parallel%nleft,12,&
         parallel%comm,status,ierr)

    ! Send to up receive from down
    call mpi_sendrecv(field0%data(1,:),field0%ny,MPI_DOUBLE_PRECISION,parallel%nup,13, &
         field0%data(field0%nx+1,:),field0%ny,MPI_DOUBLE_PRECISION,parallel%ndown,13,&
         parallel%comm,status,ierr)

    ! Send to the down, receive from up
    call mpi_sendrecv(field0%data(field0%nx,:),field0%ny,MPI_DOUBLE_PRECISION,parallel%ndown,14, &
         field0%data(0,:),field0%ny,MPI_DOUBLE_PRECISION,parallel%nup,14,&
         parallel%comm,status,ierr)
  end subroutine exchange

  ! Compute one time step of temperature evolution
  ! Arguments:
  !   curr (type(field)): current temperature values
  !   prev (type(field)): values from previous time step
  !   a (real(dp)): update equation constant
  !   dt (real(dp)): time step value
  subroutine evolve_interior(curr, prev, a, dt)

    implicit none

    type(field), intent(inout) :: curr, prev
    real(dp) :: a, dt
    integer :: i, j, nx, ny

    nx = curr%nx
    ny = curr%ny

    do j = 2, ny - 1
       do i = 2, nx - 1
          curr%data(i, j) = prev%data(i, j) + a * dt * &
               & ((prev%data(i-1, j) - 2.0 * prev%data(i, j) + &
               &   prev%data(i+1, j)) / curr%dx**2 + &
               &  (prev%data(i, j-1) - 2.0 * prev%data(i, j) + &
               &   prev%data(i, j+1)) / curr%dy**2)
       end do
    end do
  end subroutine evolve_interior

  ! Compute one time step of temperature evolution
  ! Arguments:
  !   curr (type(field)): current temperature values
  !   prev (type(field)): values from previous time step
  !   a (real(dp)): update equation constant
  !   dt (real(dp)): time step value
  ! Update only the border-dependent part
  subroutine evolve_edges(curr, prev, a, dt)

    implicit none

    type(field), intent(inout) :: curr, prev
    real(dp) :: a, dt
    integer :: i, j, nx, ny

    nx = curr%nx
    ny = curr%ny

    j = 1
    do i = 1, nx
       curr%data(i, j) = prev%data(i, j) + a * dt * &
            & ((prev%data(i-1, j) - 2.0 * prev%data(i, j) + &
            &   prev%data(i+1, j)) / curr%dx**2 + &
            &  (prev%data(i, j-1) - 2.0 * prev%data(i, j) + &
            &   prev%data(i, j+1)) / curr%dy**2)
    end do
    j = ny
    do i = 1, nx
       curr%data(i, j) = prev%data(i, j) + a * dt * &
            & ((prev%data(i-1, j) - 2.0 * prev%data(i, j) + &
            &   prev%data(i+1, j)) / curr%dx**2 + &
            &  (prev%data(i, j-1) - 2.0 * prev%data(i, j) + &
            &   prev%data(i, j+1)) / curr%dy**2)
    end do
    i = 1
    do j = 1, ny
       curr%data(i, j) = prev%data(i, j) + a * dt * &
            & ((prev%data(i-1, j) - 2.0 * prev%data(i, j) + &
            &   prev%data(i+1, j)) / curr%dx**2 + &
            &  (prev%data(i, j-1) - 2.0 * prev%data(i, j) + &
            &   prev%data(i, j+1)) / curr%dy**2)
    end do
    i = nx
    do j = 1, ny
       curr%data(i, j) = prev%data(i, j) + a * dt * &
            & ((prev%data(i-1, j) - 2.0 * prev%data(i, j) + &
            &   prev%data(i+1, j)) / curr%dx**2 + &
            &  (prev%data(i, j-1) - 2.0 * prev%data(i, j) + &
            &   prev%data(i, j+1)) / curr%dy**2)
    end do

  end subroutine evolve_edges


end module core
