program main
  use iso_c_binding
  use kmrnextf
  ! use testfn
  implicit none
#ifdef BACKEND_KMR
  include "mpif.h"
#endif
  !---------------------------------------------------------------------------
  logical(1), parameter :: kPrint = .TRUE.

  integer(8), parameter :: kDimension3 = 3
  integer(8), parameter :: kDim3_0 = 10
  integer(8), parameter :: kDim3_1 = 10
  integer(8), parameter :: kDim3_2 = 10

  integer :: rank = 0
  integer :: ierr
  type(c_ptr) :: next
  type(c_ptr) :: ds1
  integer(8)  :: sizes3(Max_Dimension_Size)
  data sizes3/kDim3_0,kDim3_1,kDim3_2,0,0,0,0,0/

  !---------------------------------------------------------------------------

  next = kmrnext_init()
#ifdef BACKEND_KMR
  call mpi_comm_rank(mpi_comm_world, rank, ierr)
#endif

  !------ Create a DataStore
  ds1 = kmrnext_create_ds(next, kDimension3)
  ierr = kmrnext_ds_set_size(ds1, sizes3)
  if (kPrint .and. rank == 0) then
     print *, 'this is test.'
  end if


  ierr = kmrnext_free_ds(ds1);
  ierr = kmrnext_finalize()
end program main
