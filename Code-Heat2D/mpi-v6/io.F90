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
    type(field), intent(in) :: curr
    integer, intent(in) :: iter
    type(parallel_data), intent(in) :: parallel

    character(len=85) :: filename
    integer :: coords(2), ierr

    ! PHDF5-related variables
    integer(hid_t) :: fapl_id        ! file access identifier
    integer(hid_t) :: file_id        ! File identifier
    integer(hid_t) :: dxpl_id        ! dataset transfer property list
    integer(hid_t) :: fspace_id      ! Dataspace identifier in file 
    integer(hid_t) :: mspace_id      ! Memspace identifier in memory
    integer(hid_t) :: dset_id        ! Dataset identifier
    integer(hsize_t), dimension(2) :: total_dims
    integer(hsize_t), dimension(2) :: dims,offset
    integer(size_t)             :: s_a = 1
    integer,dimension(1)        :: att
    
    write(filename,'(A5,I4.4,A3)')  'heat_', iter, '.h5'
      
    ! Open HDF5 Fortran interface
    call h5open_f(ierr)

    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)

    ! Create the file collectively.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, access_prp = fapl_id)

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Write attributes nx and ny
    att = curr%nx_full ; call h5ltset_attribute_int_f(file_id,"/","nx", att,s_a, ierr)
    att = curr%ny_full ; call h5ltset_attribute_int_f(file_id,"/","ny", att,s_a, ierr)

    ! Write T dataset
    total_dims = (/ curr%nx_full, curr%ny_full /)
    call h5screate_simple_f(2, total_dims, fspace_id, ierr) ! Dataset in file

    call mpi_cart_coords(parallel%comm, parallel%rank, 2, coords, ierr)
    dims = (/ curr%nx, curr%ny /)
    offset = (/ coords(1) * dims(1) , coords(2) * dims(2) /)
    call h5screate_simple_f(2, dims, mspace_id, ierr)      ! Dataset in memory
    call h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, dims, ierr)

    call h5dcreate_f(file_id, "T", H5T_NATIVE_DOUBLE, fspace_id, dset_id, ierr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE,curr%data(1:curr%nx,1:curr%ny),total_dims, ierr, &
         file_space_id=fspace_id, mem_space_id = mspace_id, xfer_prp = dxpl_id)
    call h5dclose_f(dset_id, ierr)   
    ! --------
    
    call h5sclose_f(fspace_id, ierr)
    call h5sclose_f(mspace_id, ierr)
    call h5pclose_f(fapl_id, ierr)   ! Collective access to the file
    call h5pclose_f(dxpl_id, ierr)   ! Collective access to the dataset
    call h5fclose_f(file_id, ierr)   ! The file
    call h5close_f(ierr)
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

    integer :: nx, ny
    integer :: coords(2),ierr

    ! PHDF5-related variables
    integer(hid_t) :: fapl_id        ! file access identifier
    integer(hid_t) :: file_id        ! File identifier
    integer(hid_t) :: dxpl_id        ! dataset transfer property list
    integer(hid_t) :: fspace_id      ! Dataspace identifier in file 
    integer(hid_t) :: mspace_id      ! Memspace identifier in memory
    integer(hid_t) :: dset_id        ! Dataset identifier
    integer(hsize_t), dimension(2) :: dims,offset
    integer,dimension(1)        :: att

    ! Open HDF5 Fortran interface
    call h5open_f(ierr)

    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, ierr)
    call h5pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)

    ! Open the file collectivelly.
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr, access_prp = fapl_id)     

    ! Create property list for collective dataset write
    call h5pcreate_f(H5P_DATASET_XFER_F, dxpl_id, ierr)
    call h5pset_dxpl_mpio_f(dxpl_id, H5FD_MPIO_COLLECTIVE_F, ierr)

    ! Read attributes nx and ny
    call h5ltget_attribute_int_f(file_id,"/","nx", att, ierr) ; nx = att(1)
    call h5ltget_attribute_int_f(file_id,"/","ny", att, ierr) ; ny = att(1)

    ! Set things up based on nx,ny
    call parallel_setup(parallel, nx, ny)
    call set_field_dimensions(field0, nx, ny, parallel)
    
    call mpi_cart_coords(parallel%comm, parallel%rank, 2, coords, ierr)
    dims = (/ field0%nx, field0%ny /)
    offset = (/ coords(1) * dims(1) , coords(2) * dims(2) /)
    call h5screate_simple_f(2, dims, mspace_id, ierr)      ! Dataset in memory

    ! Read T dataset
    ! The arrays for temperature field contain also a halo region
    allocate(field0%data(0:field0%nx+1, 0:field0%ny+1))

    call h5dopen_f(file_id, "T", dset_id, ierr)
    call h5dget_space_f(dset_id, fspace_id, ierr)          ! Dataset in file
    call h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, offset, dims, ierr)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, field0%data(1:field0%nx,1:field0%ny), dims, ierr, &
         file_space_id=fspace_id, mem_space_id=mspace_id, xfer_prp=dxpl_id)
    call h5dclose_f(dset_id, ierr)   
    ! ---------
    
    call h5sclose_f(fspace_id, ierr)
    call h5sclose_f(mspace_id, ierr)
    call h5pclose_f(fapl_id, ierr)   ! Collective access to the file
    call h5pclose_f(dxpl_id, ierr)   ! Collective access to the dataset
    call h5fclose_f(file_id, ierr)   ! The file
    call h5close_f(ierr)

    ! Set the boundary values
    field0%data(1:field0%nx, 0) = field0%data(1:field0%nx, 1)
    field0%data(1:field0%nx, field0%ny+1) = field0%data(1:field0%nx, field0%ny)
    field0%data(0, 0:field0%ny+1) = field0%data(1, 0:field0%ny+1)
    field0%data(field0%nx+1, 0:field0%ny+1) = field0%data(field0%nx, 0:field0%ny+1)
  end subroutine read_field

end module io
