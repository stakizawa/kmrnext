program main
  use iso_c_binding
  use kmrnextf
  ! use testfn
  implicit none
#ifdef BACKEND_KMR
  include "mpif.h"
#endif
  !---------------------------------------------------------------------------

!!!!  integer :: ierr
  type(c_ptr) :: next

  !---------------------------------------------------------------------------

  next = kmrnext_init();
  print *, "this is test."
end program main
