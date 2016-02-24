module kmrnextf
  use iso_c_binding
  implicit none

  public

  interface kmrnext_init_ff
     type(c_ptr) function kmrnext_init_ff() bind(c, name='kmrnext_init_ff')
       use iso_c_binding
     end function kmrnext_init_ff
  end interface kmrnext_init_ff
  !--------------------------------------------------------------------------

contains

  type(c_ptr) function kmrnext_init() result(zz)
    zz = kmrnext_init_ff()
  end function kmrnext_init

end module kmrnextf
