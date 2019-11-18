! Field metadata for heat equation solver
module heat
  use iso_fortran_env, only : REAL64
  implicit none

  integer, parameter :: dp = REAL64
  real(dp), parameter :: DX = 0.01, DY = 0.01  ! Fixed grid spacing

  type :: field
     integer :: nx          ! local dimension of the field
     integer :: ny
     integer :: nx_full     ! global dimension of the field
     integer :: ny_full
     real(dp) :: dx
     real(dp) :: dy
     real(dp), dimension(:,:), allocatable :: data
  end type field

  type :: parallel_data
     integer :: size
     integer :: rank
     integer :: nup, ndown, nleft, nright  ! Ranks of neighbouring MPI tasks
  end type parallel_data

contains
  ! Initialize the field type metadata
  ! Arguments:
  !   field0 (type(field)): input field
  !   nx, ny: field dimensions 
  subroutine set_field_dimensions(field0, nx, ny)
    implicit none

    type(field), intent(out) :: field0
    integer, intent(in) :: nx, ny

    integer :: nx_local, ny_local
    integer :: dims(2) = (/2, 2/)   ! Hardcoded to 2x2 topology

    nx_local = nx / dims(1)
    ny_local = ny / dims(2)

    field0%dx = DX
    field0%dy = DY
    field0%nx = nx_local
    field0%ny = ny_local
    field0%nx_full = nx
    field0%ny_full = ny
  end subroutine set_field_dimensions

  subroutine parallel_setup(parallel, nx, ny)
    use mpi

    implicit none

    type(parallel_data), intent(out) :: parallel
    integer, intent(in), optional :: nx, ny

    integer :: nx_local, ny_local
    integer :: dims(2) = (/2, 2/)   ! Hardcoded to 2x2 topology
    integer :: ierr

    call mpi_comm_size(MPI_COMM_WORLD, parallel%size, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, parallel%rank, ierr)
    
    ! Set grid dimensions
    nx_local = nx / dims(1)
    ny_local = ny / dims(2)

    ! Ensure that the grid is divisible to the MPI tasks
    if (nx_local * dims(1) /= nx) then
       write(*,*) 'Cannot divide grid evenly to processors in x-direction', &
            & nx_local, dims(1), nx
       call mpi_abort(MPI_COMM_WORLD, -2, ierr)
    end if
    if (ny_local * dims(2) /= ny) then
       write(*,*) 'Cannot divide grid evenly to processors in y-direction', &
            & ny_local, dims(2), ny
       call mpi_abort(MPI_COMM_WORLD, -2, ierr)
    end if

    ! Setting neighbours
    ! NOTE: This is hardcoded to work with 4 processes
    if (mod(parallel%rank, 2) == 1) then
       ! right boundary
       parallel%nright = MPI_PROC_NULL
       parallel%nleft = parallel%rank - 1
    else
       ! left boundary
       parallel%nright = parallel%rank + 1
       parallel%nleft = MPI_PROC_NULL
    end if


    if ((parallel%rank / 2) == 1) then
       ! bottom boundary
       parallel%ndown = MPI_PROC_NULL
       parallel%nup = parallel%rank - 2
    else
       ! top boundary
       parallel%ndown = parallel%rank + 2
       parallel%nup = MPI_PROC_NULL
    end if

    if (parallel%rank == 0) then
       write(*,'(A, I3, A3, I3)') 'Using domain decomposition', dims(1), 'x', & 
         & dims(2)
       write(*,'(A, I5, A3, I5)') 'Local domain size', nx_local,  'x', ny_local
    end if

  end subroutine parallel_setup

end module heat
