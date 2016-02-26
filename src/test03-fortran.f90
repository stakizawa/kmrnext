module test03
  use iso_c_binding
  implicit none
contains

  integer(c_int) function loader(ds, file) bind(c) result(zz)
    use iso_c_binding
    use kmrnextf
    implicit none
    type(c_ptr), intent(in), value :: ds
    !type(c_ptr), intent(in), value :: file_ptr
    character(c_char), intent(in) :: file(:)

    !integer(c_size_t) :: len_file
    !character(c_char), pointer :: file(:)

    !len_file = C_strlen(file_ptr)
    !call C_F_POINTER(file_ptr, file, [len_file])
    write (*,*) 'file: ', file
    zz = 0
  end function loader

  ! subroutine print_data_store(ds, space, count, rank)
  !   implicit none
  !   type(c_ptr)  :: ds
  !   character(*) :: space
  !   integer      :: count
  !   integer      :: rank
  !   !-------------------------------------------------------------------------

  !   if (rank /= 0) then
  !      return
  !   end if
  ! end subroutine print_data_store

end module test03

program main
  use iso_c_binding
  use kmrnextf
  use test03
  implicit none
#ifdef BACKEND_KMR
  include "mpif.h"
#endif
  !---------------------------------------------------------------------------
  logical(1), parameter :: kPrint = .TRUE.
  integer,    parameter :: kDumpCount = 20

  integer(8), parameter :: kDimension3 = 3
  integer(8), parameter :: kDim3_0 = 10
  integer(8), parameter :: kDim3_1 = 10
  integer(8), parameter :: kDim3_2 = 10

  integer                    :: rank = 0
  integer                    :: ierr
  type(c_ptr)                :: next
  type(c_ptr)                :: ds1
  integer(8)                 :: sizes3(Max_Dimension_Size)
  character(6), dimension(1) :: files

  character(c_char), pointer :: str(:)

  data sizes3/kDim3_0,kDim3_1,kDim3_2,0,0,0,0,0/
  data files/'dummy1'/

  !---------------------------------------------------------------------------

  next = kmrnext_init()
#ifdef BACKEND_KMR
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
#endif

  !------ Create a DataStore
  ds1 = kmrnext_create_ds(next, kDimension3)
  ierr = kmrnext_ds_set_size(ds1, sizes3)
  if (kPrint .and. rank == 0) then
     write (*,*) '0. Create a DataStore'
     call kmrnext_ds_string(ds1, str)
     write (*,*) '  DataStore: ', str
     write (*,*)

     !call print_data_store(ds1, '  ', kDumpCount, rank)

     ! write (*,fmt='(a)', advance='no') 'this is test.' ! TODO
     ! write (*,*) 'another string'
  end if

  !------ Load data contents from a file
  ierr = kmrnext_ds_load_files(ds1, files, 1, loader)


  ierr = kmrnext_free_ds(ds1)
  ierr = kmrnext_finalize()
end program main
