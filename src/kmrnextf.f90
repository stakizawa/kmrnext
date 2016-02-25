module kmrnextf
  use iso_c_binding
  implicit none

  public

  integer(8), parameter :: Max_Dimension_Size = 8

  interface
     type(c_ptr) function C_kmrnext_init() bind(c, name='KMRNEXT_init_ff')
       use iso_c_binding
       implicit none
     end function C_kmrnext_init

     subroutine C_kmrnext_finalize() bind(c, name='KMRNEXT_finalize')
       use iso_c_binding
       implicit none
     end subroutine C_kmrnext_finalize

     type(c_ptr) function C_kmrnext_create_ds(next, size) &
          bind(c, name='KMRNEXT_create_ds')
       use iso_c_binding
       implicit none
       type(c_ptr),     value, intent(in) :: next
       integer(c_long), value, intent(in) :: size
     end function C_kmrnext_create_ds

     subroutine C_kmrnext_free_ds(ds) bind(c, name='KMRNEXT_free_ds')
       use iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: ds
     end subroutine C_kmrnext_free_ds

     subroutine C_kmrnext_ds_set_size(ds, size) &
          bind(c, name='KMRNEXT_ds_set_size')
       use iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: ds
       type(c_ptr), value, intent(in) :: size
     end subroutine C_kmrnext_ds_set_size

     ! kmrnext_ds_load_files
     ! kmrnext_ds_add
     ! kmrnext_ds_get
     ! kmrnext_ds_get_view
     ! kmrnext_ds_map
     ! kmrnext_ds_count
     ! kmrnext_ds_dump
     ! kmrnext_ds_string

     ! kmrnext_create_key
     ! kmrnext_free_key
     ! kmrnext_key_set_size
     ! kmrnext_key_set
     ! kmrnext_key_string

     ! kmrnext_create_data
     ! kmrnext_free_data
     ! kmrnext_data_value
     ! kmrnext_data_size

     ! kmrnext_create_dp
     ! kmrnext_free_dp
     ! kmrnext_dp_key
     ! kmrnext_dp_data

     ! kmrnext_create_view
     ! kmrnext_free_view
     ! kmrnext_view_set
     ! kmrnext_view_string

     ! kmrnext_free_datapacks

  end interface

  !--------------------------------------------------------------------------

contains

  type(c_ptr) function kmrnext_init() result(zz)
    zz = C_kmrnext_init()
  end function kmrnext_init

  integer function kmrnext_finalize() result(zz)
    call C_kmrnext_finalize()
    zz = 0
  end function kmrnext_finalize

  type(c_ptr) function kmrnext_create_ds(next, size) result(zz)
    type(c_ptr),     value, intent(in) :: next
    integer(c_long), value, intent(in) :: size
    zz = C_kmrnext_create_ds(next, size)
  end function kmrnext_create_ds

  integer(c_int) function kmrnext_free_ds(ds) result(zz)
    type(c_ptr), value, intent(in) :: ds
    call C_kmrnext_free_ds(ds);
    zz = 0
  end function kmrnext_free_ds

  integer(c_int) function kmrnext_ds_set_size(ds, size) result(zz)
    type(c_ptr), value, intent(in)         :: ds
    integer(c_long),         intent(in), target :: size(Max_Dimension_Size)
    call C_kmrnext_ds_set_size(ds, C_LOC(size))
    zz = 0
  end function kmrnext_ds_set_size

  ! kmrnext_ds_load_files
  ! kmrnext_ds_add
  ! kmrnext_ds_get
  ! kmrnext_ds_get_view
  ! kmrnext_ds_map
  ! kmrnext_ds_count
  ! kmrnext_ds_dump
  ! kmrnext_ds_string

  ! kmrnext_create_key
  ! kmrnext_free_key
  ! kmrnext_key_set_size
  ! kmrnext_key_set
  ! kmrnext_key_string

  ! kmrnext_create_data
  ! kmrnext_free_data
  ! kmrnext_data_value
  ! kmrnext_data_size

  ! kmrnext_create_dp
  ! kmrnext_free_dp
  ! kmrnext_dp_key
  ! kmrnext_dp_data

  ! kmrnext_create_view
  ! kmrnext_free_view
  ! kmrnext_view_set
  ! kmrnext_view_string

  ! kmrnext_free_datapacks

end module kmrnextf
