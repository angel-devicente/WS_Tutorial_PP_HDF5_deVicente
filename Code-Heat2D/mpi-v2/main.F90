! Heat equation solver in 2D.

program heat_solve
  use heat
  use core
  use io
  use setup
  use utilities
  use mpi

  implicit none

  real(dp), parameter :: a = 0.5                  ! Diffusion constant
  type(field) :: current, previous                ! Current and previus temperature fields

  real(dp) :: dt                                  ! Time step
  integer :: nsteps                               ! Number of time steps
  integer, parameter :: image_interval = 5        ! Image output interval

  type(parallel_data) :: parallelization
  integer :: ierr

  integer :: iter, iter0, saved

  call mpi_init(ierr)

  call initialize(current, previous, nsteps, parallelization, iter0)

  ! Draw the picture of the initial state
  saved = 0
  call write_field(current, saved, parallelization)
  iter0 = iter0 + 1

  ! Largest stable time step
  dt = current%dx**2 * current%dy**2 / &
       & (2.0 * a * (current%dx**2 + current%dy**2))

  ! Main iteration loop, save a picture every
  ! image_interval steps

  do iter = iter0, iter0 + nsteps
     call exchange(previous, parallelization)
     call evolve_interior(current, previous, a, dt)
     call evolve_edges(current,previous, a, dt)
     if (mod(iter, image_interval) == 0) then
        saved = saved + 1
        call write_field(current, saved, parallelization)
     end if
     call swap_fields(current, previous)
  end do

  if (parallelization % rank == 0) then
     write(*,'(A,G0)') 'Reference value at 5,5: ', previous % data(5,5)
  end if

  call finalize(current, previous)

  call mpi_finalize(ierr)

end program heat_solve
