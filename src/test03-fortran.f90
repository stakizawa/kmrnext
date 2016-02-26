module test03_constant
  implicit none

  logical(1), parameter :: kPrint = .TRUE.
  integer,    parameter :: kDumpCount = 20

  integer(8), parameter :: kDimension3 = 3
  integer(8), parameter :: kDim3_0 = 10
  integer(8), parameter :: kDim3_1 = 10
  integer(8), parameter :: kDim3_2 = 10
end module test03_constant

module test03
  use iso_c_binding
  use kmrnextf
  use test03_constant
  implicit none
contains

  integer(c_int) function loader(ds, file_ptr) bind(c) result(zz)
    implicit none
    type(c_ptr), intent(in), value :: ds
    type(c_ptr), intent(in), value :: file_ptr

    integer                    :: ierr, i, j, k
    integer(c_size_t)          :: len_file
    character(c_char), pointer :: file(:)
    type(c_ptr)                :: key, d
    integer(c_long), target    :: val

    len_file = C_strlen(file_ptr)
    call C_F_POINTER(file_ptr, file, [len_file])

    key = kmrnext_create_key(kDimension3)
    do i = 1, kDim3_0
       ierr = kmrnext_key_set(key, 1, int(i,kind=c_long))
       do j = 1, kDim3_1
          ierr = kmrnext_key_set(key, 2, int(j,kind=c_long))
          do k = 1, kDim3_2
             ierr = kmrnext_key_set(key, 3, int(k,kind=c_long))
             val = (i-1) * (j-1) * (k-1)
             ! 8 is sizeof(c_long)
             d = kmrnext_create_data(C_LOC(val), int(8,kind=c_long))
             ierr = kmrnext_ds_add(ds, key, d)
             ierr = kmrnext_free_data(d)
          end do
       end do
    end do
    ierr = kmrnext_free_key(key)
    zz = 0
  end function loader

  type(c_ptr) function dumper(dp) bind(c) result(zz)
    implicit none
    type(c_ptr), intent(in), value :: dp

    character(16, c_char), target :: string = 'ABCD'

    string = string//C_NULL_CHAR
    zz = C_LOC(string(1:1))
  end function dumper

  subroutine print_data_store(ds, space, count, rank)
    implicit none
    type(c_ptr),  intent(in) :: ds
    character(*), intent(in) :: space
    integer,      intent(in) :: count
    integer,      intent(in) :: rank

    integer(c_long)            :: ds_count
    character(c_char), pointer :: dumped_str(:)

    ds_count = kmrnext_ds_count(ds)
    write (*,*) 'ds_count: ', ds_count
    call kmrnext_ds_dump(ds, dumped_str, dumper)
    write (*,*) 'dumped: ', dumped_str

    if (rank /= 0) then
       return
    end if
  end subroutine print_data_store

end module test03

program main
  use iso_c_binding
  use kmrnextf
  use test03_constant
  use test03
  implicit none
#ifdef BACKEND_KMR
  include "mpif.h"
#endif
  !---------------------------------------------------------------------------
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

     ! write (*,fmt='(a)', advance='no') 'this is test.' ! TODO
     ! write (*,*) 'another string'
  end if

  !------ Load data contents from a file
  ierr = kmrnext_ds_load_files(ds1, files, 1, loader)
  if (kPrint) then
     if (rank == 0) then
        write (*,*) '1. Load data to a DataStore'
     end if
     call print_data_store(ds1, '  ', kDumpCount, rank)
     if (rank == 0) then
        write (*,*) ''
     end if
  end if


  ierr = kmrnext_free_ds(ds1)
  ierr = kmrnext_finalize()
end program main
