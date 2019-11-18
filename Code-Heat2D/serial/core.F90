! Main solver routines for heat equation solver
module core
  use heat

contains

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
