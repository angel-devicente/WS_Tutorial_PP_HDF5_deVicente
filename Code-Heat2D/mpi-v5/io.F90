! I/O routines for heat equation solver
module io
  use heat
  use mpi
  use h5lt

  implicit none
  
contains

  ! Output routine, saves the temperature distribution as a .h5 file
  ! Arguments:
  !   curr (type(field)): variable with the temperature data
  !   iter (integer): index of the time step
  subroutine write_field(curr, iter, parallel)
    implicit none

    type(field), intent(in) :: curr
    integer, intent(in) :: iter
    type(parallel_data), intent(in) :: parallel

    character(len=85) :: filename
    real(dp), dimension(:,:), allocatable, target :: full_data

    integer :: coords(2)
    integer :: ix, jy
    integer :: p, ierr

    ! HDF5-related variables
    integer(hid_t) :: file_id        ! File identifier
    integer(hsize_t), dimension(2) :: total_dims
    integer(size_t)             :: s_a = 1
    integer,dimension(1)        :: att
    
    if (parallel%rank == 0) then
       allocate(full_data(curr%nx_full, curr%ny_full))
       ! Copy rand #0 data to the global array
       full_data(1:curr%nx, 1:curr%ny) = curr%data(1:curr%nx, 1:curr%ny)

       ! Receive data from other ranks
       do p = 1, parallel%size - 1
          call mpi_cart_coords(parallel%comm, p, 2, coords, ierr)
          ix = coords(1) * curr%nx + 1
          jy = coords(2) * curr%ny + 1
          call mpi_recv(full_data(ix, jy), 1, parallel%subarraytype, p, 22, &
               & parallel%comm, MPI_STATUS_IGNORE, ierr)
       end do

       write(filename,'(A5,I4.4,A3)')  'heat_', iter, '.h5'

       ! Open HDF5 Fortran interface
       call h5open_f(ierr)

       ! Create the file 
       call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)

       ! Write attributes nx and ny
       att = curr%nx_full ; call h5ltset_attribute_int_f(file_id,"/","nx", att,s_a, ierr)
       att = curr%ny_full ; call h5ltset_attribute_int_f(file_id,"/","ny", att,s_a, ierr)

       ! Write T dataset
       total_dims = (/ curr%nx_full, curr%ny_full /)
       call h5ltmake_dataset_f(file_id,"T",2,total_dims,H5T_NATIVE_DOUBLE,full_data,ierr) 
    
       call h5fclose_f(file_id, ierr)   ! Close
       call h5close_f(ierr)

       deallocate(full_data)
    else
       ! Send data
       call mpi_send(curr%data, 1, parallel%subarraytype, 0, 22, &
            & parallel%comm, ierr)
    end if

  end subroutine write_field


  ! Reads the temperature distribution from an input file
  ! Arguments:
  !   field0 (type(field)): field variable that will store the
  !                         read data
  !   filename (char): name of the input file
  ! Note that this version assumes the input data to be in C memory layout
  subroutine read_field(field0, filename, parallel)
    type(field), intent(out) :: field0
    character(len=85), intent(in) :: filename
    type(parallel_data), intent(out) :: parallel

    integer :: nx, ny, p, ierr
    real(dp), dimension(:,:), allocatable :: full_data
    integer :: coords(2)
    integer :: ix, jy

    ! HDF5-related variables
    integer(hid_t) :: file_id        ! File identifier
    integer(hsize_t), dimension(2) :: total_dims
    integer,dimension(1)        :: att

    ! Reading of the attributes is done by all processes, as in the
    ! previous versions with .dat files.
    
    ! Open HDF5 Fortran interface
    call h5open_f(ierr)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)
    
    ! Read attributes nx and ny
    call h5ltget_attribute_int_f(file_id,"/","nx", att, ierr) ; nx = att(1)
    call h5ltget_attribute_int_f(file_id,"/","ny", att, ierr) ; ny = att(1)
    total_dims = (/ nx, ny /)

    call h5fclose_f(file_id, ierr)   ! Close

    
    call parallel_setup(parallel, nx, ny)
    call set_field_dimensions(field0, nx, ny, parallel)

    ! The arrays for temperature field contain also a halo region
    allocate(field0%data(0:field0%nx+1, 0:field0%ny+1))

    if (parallel%rank == 0) then
       allocate(full_data(nx, ny))
       ! Read the data
       call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)
       call h5ltread_dataset_double_f(file_id, "T", full_data, total_dims, ierr)
       call h5fclose_f(file_id, ierr)   ! Close
       call h5close_f(ierr)
       
       ! Copy own local part
       field0%data(1:field0%nx, 1:field0%ny) = full_data(1:field0%nx, 1:field0%ny)
       ! Send to other process
       do p=1, parallel%size - 1
          call mpi_cart_coords(parallel%comm, p, 2, coords, ierr)
          ix = coords(1) * field0%nx + 1
          jy = coords(2) * field0%ny + 1
          call mpi_send(full_data(ix, jy), 1, parallel%subarraytype, p, 44, &
               &   parallel%comm, ierr)
       end do
    else
       ! Receive data
       call mpi_recv(field0%data, 1, parallel%subarraytype, 0, 44, &
            &        parallel%comm, MPI_STATUS_IGNORE, ierr)
    end if

    ! Set the boundary values
    field0%data(1:field0%nx, 0) = field0%data(1:field0%nx, 1)
    field0%data(1:field0%nx, field0%ny+1) = field0%data(1:field0%nx, field0%ny)
    field0%data(0, 0:field0%ny+1) = field0%data(1, 0:field0%ny+1)
    field0%data(field0%nx+1, 0:field0%ny+1) = field0%data(field0%nx, 0:field0%ny+1)

    close(10)
    if (parallel%rank == 0) then
       deallocate(full_data)
    end if

  end subroutine read_field

end module io
