! The backend runtime (SERIAL, KMR)
#define BACKEND_KMR 1

module kmrnextf
  use iso_c_binding
#ifdef BACKEND_KMR
  use mpi
#endif
  implicit none

  integer(8), parameter :: Max_Dimension_Size = 8

  type, bind(c) :: datapacks
     integer(c_size_t) :: count
     type(c_ptr)       :: data
  end type datapacks

  type, bind(c) :: mapenvf
     integer(c_int) :: rank
#ifdef BACKEND_KMR
     integer(c_int) :: mpi_comm
#endif
     type(c_ptr)    :: p
  end type mapenvf

  abstract interface
     integer(c_int) function kmrnext_loadfn(ds, file_ptr) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
       type(c_ptr), intent(in), value :: file_ptr
       !character(c_char), intent(in) :: file(:)
     end function kmrnext_loadfn

     integer(c_int) function kmrnext_mapfn(ids, ods, key, dps, env) bind(c)
       use iso_c_binding
       import datapacks
       import mapenvf
       implicit none
       type(c_ptr),     intent(in), value :: ids
       type(c_ptr),     intent(in), value :: ods
       type(c_ptr),     intent(in), value :: key
       type(datapacks), intent(in), value :: dps
       type(mapenvf),   intent(in), value :: env
     end function kmrnext_mapfn

     type(c_ptr) function kmrnext_dumpfn(dp) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: dp
     end function kmrnext_dumpfn
  end interface

  interface
     type(c_ptr) function C_kmrnext_init() bind(c, name='KMRNEXT_init_ff')
       use iso_c_binding
       implicit none
     end function C_kmrnext_init

     subroutine C_kmrnext_finalize() bind(c, name='KMRNEXT_finalize')
       use iso_c_binding
       implicit none
     end subroutine C_kmrnext_finalize

     subroutine C_kmrnext_enable_profile(next) &
          bind(c, name='KMRNEXT_enable_profile')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: next
     end subroutine C_kmrnext_enable_profile

     subroutine C_kmrnext_disable_profile(next) &
          bind(c, name='KMRNEXT_disable_profile')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: next
     end subroutine C_kmrnext_disable_profile

     logical(c_bool) function C_kmrnext_profile(next) &
          bind(c, name='KMRNEXT_profile')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: next
     end function C_kmrnext_profile

#ifdef BACKEND_KMR
     integer(c_long) function C_kmrnext_nprocs(next) &
          bind(c, name='KMRNEXT_nprocs')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: next
     end function C_kmrnext_nprocs

     integer(c_long) function C_kmrnext_rank(next) &
          bind(c, name='KMRNEXT_rank')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: next
     end function C_kmrnext_rank
#endif

     type(c_ptr) function C_kmrnext_create_ds(next, size) &
          bind(c, name='KMRNEXT_create_ds')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: next
       integer(c_size_t), intent(in), value :: size
     end function C_kmrnext_create_ds

     subroutine C_kmrnext_free_ds(ds) bind(c, name='KMRNEXT_free_ds')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
     end subroutine C_kmrnext_free_ds

     subroutine C_kmrnext_ds_set_size(ds, size) &
          bind(c, name='KMRNEXT_ds_set_size')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
       type(c_ptr), intent(in), value :: size
     end subroutine C_kmrnext_ds_set_size

     subroutine C_kmrnext_ds_load_files(ds, files, nfiles, l) &
          bind(c, name='KMRNEXT_ds_load_files')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: ds
       type(c_ptr),       intent(in), value :: files
       integer(c_size_t), intent(in), value :: nfiles
       type(c_funptr),    intent(in), value :: l
     end subroutine C_kmrnext_ds_load_files

     subroutine C_kmrnext_ds_add(ds, key, dat) bind(c, name='KMRNEXT_ds_add')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
       type(c_ptr), intent(in), value :: key
       type(c_ptr), intent(in), value :: dat
     end subroutine C_kmrnext_ds_add

     type(c_ptr) function C_kmrnext_ds_get(ds, key) &
          bind(c, name='KMRNEXT_ds_get')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
       type(c_ptr), intent(in), value :: key
     end function C_kmrnext_ds_get

     type(datapacks) function C_kmrnext_ds_get_view(ds, key, view) &
          bind(c, name='KMRNEXT_ds_get_view')
       use iso_c_binding
       import datapacks
       implicit none
       type(c_ptr), intent(in), value :: ds
       type(c_ptr), intent(in), value :: key
       type(c_ptr), intent(in), value :: view
     end function C_kmrnext_ds_get_view

     type(c_ptr) function C_kmrnext_ds_remove(ds, key) &
          bind(c, name='KMRNEXT_ds_remove')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
       type(c_ptr), intent(in), value :: key
     end function C_kmrnext_ds_remove

     subroutine C_kmrnext_ds_map(ids, ods, view, m, p) &
          bind(c, name='KMRNEXT_ds_map_ff')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: ids
       type(c_ptr),       intent(in), value :: ods
       type(c_ptr),       intent(in), value :: view
       type(c_funptr),    intent(in), value :: m
       type(c_ptr),       intent(in), value :: p
     end subroutine C_kmrnext_ds_map

     integer(c_long) function C_kmrnext_ds_count(ds) &
          bind(c, name='KMRNEXT_ds_count')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
     end function C_kmrnext_ds_count

     type(c_ptr) function C_kmrnext_ds_dump(ds, d) &
          bind(c, name='KMRNEXT_ds_dump')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: ds
       type(c_funptr),    intent(in), value :: d
     end function C_kmrnext_ds_dump

     type(c_ptr) function C_kmrnext_ds_string(ds) &
          bind(c, name='KMRNEXT_ds_string')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
     end function C_kmrnext_ds_string

#ifdef BACKEND_KMR
     subroutine C_kmrnext_ds_set_split(ds, split) &
          bind(c, name='KMRNEXT_ds_set_split')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
       type(c_ptr), intent(in), value :: split
     end subroutine C_kmrnext_ds_set_split

     type(c_ptr) function C_kmrnext_ds_get_split(ds) &
          bind(c, name='KMRNEXT_ds_get_split')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
     end function C_kmrnext_ds_get_split

     subroutine C_kmrnext_ds_collate(ds) &
          bind(c, name='KMRNEXT_ds_collate')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
     end subroutine C_kmrnext_ds_collate

     logical(c_bool) function C_kmrnext_ds_collated(ds) &
          bind(c, name='KMRNEXT_ds_collated')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
     end function C_kmrnext_ds_collated
#endif

     type(c_ptr) function C_kmrnext_ds_duplicate(ds) &
          bind(c, name='KMRNEXT_ds_duplicate')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: ds
     end function C_kmrnext_ds_duplicate

     type(c_ptr) function C_kmrnext_create_key(size) &
          bind(c, name='KMRNEXT_create_key')
       use iso_c_binding
       implicit none
       integer(c_size_t), intent(in), value :: size
     end function C_kmrnext_create_key

     subroutine C_kmrnext_free_key(key) bind(c, name='KMRNEXT_free_key')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: key
     end subroutine C_kmrnext_free_key

     subroutine C_kmrnext_key_set_size(key, size) &
          bind(c, name = 'KMRNEXT_key_set_size')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: key
       type(c_ptr), intent(in), value :: size
     end subroutine C_kmrnext_key_set_size

     subroutine C_kmrnext_key_set(key, dim, val) &
          bind(c, name='KMRNEXT_key_set')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: key
       integer(c_size_t), intent(in), value :: dim
       integer(c_size_t), intent(in), value :: val
     end subroutine C_kmrnext_key_set

     type(c_ptr) function C_kmrnext_key_string(key) &
          bind(c, name='KMRNEXT_key_string')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: key
     end function C_kmrnext_key_string

     type(c_ptr) function C_kmrnext_create_data(val, size) &
          bind(c, name='KMRNEXT_create_data')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: val
       integer(c_size_t), intent(in), value :: size
     end function C_kmrnext_create_data

     subroutine C_kmrnext_free_data(dat) bind(c, name='KMRNEXT_free_data')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: dat
     end subroutine C_kmrnext_free_data

     type(c_ptr) function C_kmrnext_data_value(dat) &
          bind(c, name='KMRNEXT_data_value')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: dat
     end function C_kmrnext_data_value

     integer(c_size_t) function C_kmrnext_data_size(dat) &
          bind(c, name='KMRNEXT_data_size')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: dat
     end function C_kmrnext_data_size

     type(c_ptr) function C_kmrnext_create_dp(key, dat) &
          bind(c, name='KMRNEXT_create_dp')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: key
       type(c_ptr), intent(in), value :: dat
     end function C_kmrnext_create_dp

     subroutine C_kmrnext_free_dp(dp) bind(c, name='KMRNEXT_free_dp')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: dp
     end subroutine C_kmrnext_free_dp

     type(c_ptr) function C_kmrnext_dp_key(dp) bind(c, name='KMRNEXT_dp_key')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: dp
     end function C_kmrnext_dp_key

     type(c_ptr) function C_kmrnext_dp_data(dp) bind(c, name='KMRNEXT_dp_data')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: dp
     end function C_kmrnext_dp_data

     type(c_ptr) function C_kmrnext_create_view(size) &
          bind(c, name='KMRNEXT_create_view')
       use iso_c_binding
       implicit none
       integer(c_size_t), intent(in), value :: size
     end function C_kmrnext_create_view

     subroutine C_kmrnext_free_view(view) bind(c, name='KMRNEXT_free_view')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: view
     end subroutine C_kmrnext_free_view

     subroutine C_kmrnext_view_set(view, size) &
          bind(c, name='KMRNEXT_view_set')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: view
       type(c_ptr), intent(in), value :: size
     end subroutine C_kmrnext_view_set

     type(c_ptr) function C_kmrnext_view_string(view) &
          bind(c, name='KMRNEXT_view_string')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: view
     end function C_kmrnext_view_string

     subroutine C_kmrnext_free_datapacks(dps) &
          bind(c, name='KMRNEXT_free_datapacks')
       use iso_c_binding
       import datapacks
       implicit none
       type(datapacks), intent(in), value :: dps
     end subroutine C_kmrnext_free_datapacks

     integer(c_size_t) function C_strlen(string) bind(c, name='strlen')
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: string
     end function C_strlen

#if 0
     subroutine C_print_string(string) bind(c, name='print_string')
       use iso_c_binding
       implicit none
       character(c_char), intent(in) :: string(*)
     end subroutine C_print_string

     subroutine C_print_strings(strings, size) bind(c, name='print_strings')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: strings
       integer(c_size_t), intent(in), value :: size
     end subroutine C_print_strings
#endif
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

  subroutine kmrnext_enable_profile(next)
    type(c_ptr), intent(in), value :: next
    call C_kmrnext_enable_profile(next)
  end subroutine kmrnext_enable_profile

  subroutine kmrnext_disable_profile(next)
    type(c_ptr), intent(in), value :: next
    call C_kmrnext_disable_profile(next)
  end subroutine kmrnext_disable_profile

  logical(c_bool) function kmrnext_profile(next) result(zz)
    type(c_ptr), intent(in), value :: next
    zz = C_kmrnext_profile(next)
  end function kmrnext_profile

#ifdef BACKEND_KMR
  integer(c_long) function kmrnext_nprocs(next) result(zz)
    type(c_ptr), intent(in), value :: next
    zz = C_kmrnext_nprocs(next)
  end function kmrnext_nprocs

  integer(c_long) function kmrnext_rank(next) result(zz)
    type(c_ptr), intent(in), value :: next
    zz = C_kmrnext_rank(next)
  end function kmrnext_rank
#endif

  type(c_ptr) function kmrnext_create_ds(next, size) result(zz)
    type(c_ptr),       intent(in), value :: next
    integer(c_size_t), intent(in), value :: size
    zz = C_kmrnext_create_ds(next, size)
  end function kmrnext_create_ds

  integer function kmrnext_free_ds(ds) result(zz)
    type(c_ptr), intent(in), value :: ds
    call C_kmrnext_free_ds(ds)
    zz = 0
  end function kmrnext_free_ds

  integer function kmrnext_ds_set_size(ds, size) result(zz)
    type(c_ptr),     intent(in), value  :: ds
    integer(c_long), intent(in), target :: size(Max_Dimension_Size)
    call C_kmrnext_ds_set_size(ds, C_LOC(size))
    zz = 0
  end function kmrnext_ds_set_size

  integer function kmrnext_ds_load_files(ds, files, nfiles, l) result(zz)
    type(c_ptr),  intent(in),  value   :: ds
    integer,      intent(in)           :: nfiles
    character(*), intent(in)           :: files(nfiles)
    procedure(kmrnext_loadfn), bind(c) :: l

    character(c_char), allocatable, target :: cstring(:,:)
    type(c_ptr),       allocatable, target :: carray(:)
    integer                                :: i, j
    integer(c_size_t)                      :: array_len

    allocate(cstring(LEN(files(1))+1, nfiles))
    array_len = nfiles
    allocate(carray(nfiles))

    do i = 1, nfiles
       do j = 1, LEN_TRIM(files(i))
          cstring(j,i) = files(i)(j:j)
       end do
       cstring(LEN_TRIM(files(i))+1, i) = C_NULL_CHAR
       carray(i) = C_LOC(cstring(1,i))
    end do
    call C_kmrnext_ds_load_files(ds, C_LOC(carray), array_len, C_FUNLOC(l))

    deallocate(carray)
    deallocate(cstring)
    zz = 0
  end function kmrnext_ds_load_files

  integer function kmrnext_ds_add(ds, key, dat) result(zz)
    type(c_ptr), intent(in), value :: ds
    type(c_ptr), intent(in), value :: key
    type(c_ptr), intent(in), value :: dat
    call C_kmrnext_ds_add(ds, key, dat)
    zz = 0
  end function kmrnext_ds_add

  type(c_ptr) function kmrnext_ds_get(ds, key) result(zz)
    type(c_ptr), intent(in), value :: ds
    type(c_ptr), intent(in), value :: key
    zz = C_kmrnext_ds_get(ds, key)
  end function kmrnext_ds_get

  type(datapacks) function kmrnext_ds_get_view(ds, key, view) result(zz)
    type(c_ptr), intent(in), value :: ds
    type(c_ptr), intent(in), value :: key
    type(c_ptr), intent(in), value :: view
    zz = C_kmrnext_ds_get_view(ds, key, view)
  end function kmrnext_ds_get_view

  type(c_ptr) function kmrnext_ds_remove(ds, key) result(zz)
    type(c_ptr), intent(in), value :: ds
    type(c_ptr), intent(in), value :: key
    zz = C_kmrnext_ds_remove(ds, key)
  end function kmrnext_ds_remove

  integer function kmrnext_ds_map(ids, ods, view, m) result(zz)
    type(c_ptr), intent(in),  value    :: ids
    type(c_ptr), intent(in),  value    :: ods
    type(c_ptr), intent(in),  value    :: view
    procedure(kmrnext_mapfn), bind(c)  :: m
    call C_kmrnext_ds_map(ids, ods, view, C_FUNLOC(m), C_NULL_PTR)
    zz = 0
  end function kmrnext_ds_map

  integer(c_long) function kmrnext_ds_count(ds) result(zz)
    type(c_ptr), intent(in), value :: ds
    zz = C_kmrnext_ds_count(ds)
  end function kmrnext_ds_count

  subroutine kmrnext_ds_dump(ds, ostring, d)
    type(c_ptr),  intent(in),  value   :: ds
    character(c_char), intent(out), pointer :: ostring(:)
    procedure(kmrnext_dumpfn), bind(c) :: d

    type(c_ptr)       :: c_str
    integer(c_size_t) :: len_str

    c_str = C_kmrnext_ds_dump(ds, C_FUNLOC(d))
    len_str = C_strlen(c_str)
    call C_F_POINTER(c_str, ostring, [len_str])
  end subroutine kmrnext_ds_dump

  subroutine kmrnext_ds_string(ds, ostring)
    type(c_ptr),       intent(in),  value   :: ds
    character(c_char), intent(out), pointer :: ostring(:)
    type(c_ptr)       :: c_str
    integer(c_size_t) :: len_str
    c_str = C_kmrnext_ds_string(ds)
    len_str = C_strlen(c_str)
    call C_F_POINTER(c_str, ostring, [len_str])
  end subroutine kmrnext_ds_string

#ifdef BACKEND_KMR
  subroutine kmrnext_ds_set_split(ds, split)
    type(c_ptr), intent(in), value :: ds
    type(c_ptr), intent(in), value :: split
    call C_kmrnext_ds_set_split(ds, split)
  end subroutine kmrnext_ds_set_split

  type(c_ptr) function kmrnext_ds_get_split(ds) result(zz)
    type(c_ptr), intent(in), value :: ds
    zz = C_kmrnext_ds_get_split(ds)
  end function kmrnext_ds_get_split

  subroutine kmrnext_ds_collate(ds)
    type(c_ptr), intent(in), value :: ds
    call C_kmrnext_ds_collate(ds)
  end subroutine kmrnext_ds_collate

  logical(c_bool) function kmrnext_ds_collated(ds) result(zz)
    type(c_ptr), intent(in), value :: ds
    zz = C_kmrnext_ds_collated(ds)
  end function kmrnext_ds_collated
#endif

  type(c_ptr) function kmrnext_ds_duplicate(ds) result(zz)
    type(c_ptr), intent(in), value :: ds
    zz = C_kmrnext_ds_duplicate(ds)
  end function kmrnext_ds_duplicate

  type(c_ptr) function kmrnext_create_key(size) result(zz)
    integer(c_size_t), intent(in), value :: size
    zz = C_kmrnext_create_key(size)
  end function kmrnext_create_key

  integer function kmrnext_free_key(key) result(zz)
    type(c_ptr), intent(in), value :: key
    call C_kmrnext_free_key(key)
    zz = 0
  end function kmrnext_free_key

  integer function kmrnext_key_set_size(key, size) result(zz)
    type(c_ptr),     intent(in), value  :: key
    integer(c_long), intent(in), target :: size(Max_Dimension_Size)
    call C_kmrnext_key_set_size(key, C_LOC(size))
    zz = 0
  end function kmrnext_key_set_size

  integer function kmrnext_key_set(key, dim, val) result(zz)
    type(c_ptr),       intent(in), value :: key
    integer,           intent(in), value :: dim
    integer(c_size_t), intent(in), value :: val
    integer(c_size_t) :: dim8
    dim8 = dim - 1
    call C_kmrnext_key_set(key, dim8, val-1)
    zz = 0
  end function kmrnext_key_set

  subroutine kmrnext_key_string(key, ostring)
    type(c_ptr),       intent(in),  value   :: key
    character(c_char), intent(out), pointer :: ostring(:)
    type(c_ptr)       :: c_str
    integer(c_size_t) :: len_str
    c_str = C_kmrnext_key_string(key)
    len_str = C_strlen(c_str)
    call C_F_POINTER(c_str, ostring, [len_str])
  end subroutine kmrnext_key_string

  type(c_ptr) function kmrnext_create_data(val, size) result(zz)
    type(c_ptr),       intent(in), value :: val
    integer(c_size_t), intent(in), value :: size
    zz = C_kmrnext_create_data(val, size)
  end function kmrnext_create_data

  integer function kmrnext_free_data(dat) result(zz)
    type(c_ptr), intent(in), value :: dat
    call C_kmrnext_free_data(dat)
    zz = 0
  end function kmrnext_free_data

  type(c_ptr) function kmrnext_data_value(dat) result(zz)
    type(c_ptr), intent(in), value :: dat
    zz = C_kmrnext_data_value(dat)
  end function kmrnext_data_value

  integer(c_size_t) function kmrnext_data_size(dat) result(zz)
    type(c_ptr), intent(in), value :: dat
    zz = C_kmrnext_data_size(dat)
  end function kmrnext_data_size

  type(c_ptr) function kmrnext_create_dp(key, dat) result(zz)
    type(c_ptr), intent(in), value :: key
    type(c_ptr), intent(in), value :: dat
    zz = C_kmrnext_create_dp(key, dat)
  end function kmrnext_create_dp

  integer function kmrnext_free_dp(dp) result(zz)
    type(c_ptr), intent(in), value :: dp
    call C_kmrnext_free_dp(dp)
    zz = 0
  end function kmrnext_free_dp

  type(c_ptr) function kmrnext_dp_key(dp) result(zz)
    type(c_ptr), intent(in), value :: dp
    zz = C_kmrnext_dp_key(dp)
  end function kmrnext_dp_key

  type(c_ptr) function kmrnext_dp_data(dp) result(zz)
    type(c_ptr),       intent(in), value :: dp
    zz = C_kmrnext_dp_data(dp)
  end function kmrnext_dp_data

  type(c_ptr) function kmrnext_create_view(size) result(zz)
    integer(c_size_t), intent(in), value :: size
    zz = C_kmrnext_create_view(size)
  end function kmrnext_create_view

  integer function kmrnext_free_view(view) result(zz)
    type(c_ptr), intent(in), value :: view
    call C_kmrnext_free_view(view)
    zz = 0
  end function kmrnext_free_view

  integer function kmrnext_view_set(view, val) result(zz)
    type(c_ptr),     intent(in), value  :: view
    logical(c_bool), intent(in), target :: val(Max_Dimension_Size)
    call C_kmrnext_view_set(view, C_LOC(val))
    zz = 0
  end function kmrnext_view_set

  subroutine kmrnext_view_string(view, ostring)
    type(c_ptr),       intent(in),  value   :: view
    character(c_char), intent(out), pointer :: ostring(:)
    type(c_ptr)       :: c_str
    integer(c_size_t) :: len_str
    c_str = C_kmrnext_view_string(view)
    len_str = C_strlen(c_str)
    call C_F_POINTER(c_str, ostring, [len_str])
  end subroutine kmrnext_view_string

  integer function kmrnext_free_datapacks(dps) result(zz)
    type(datapacks), intent(in), value :: dps
    call C_kmrnext_free_datapacks(dps)
    zz = 0
  end function kmrnext_free_datapacks

#if 0
  subroutine print_string(string)
    ! Use like the following
    !   character(16) :: str = 'string0'
    !   call print_string(str)
    character(*), intent(in) :: string
    call C_print_string(string//C_NULL_CHAR)
  end subroutine print_string

  subroutine print_strings(strings, nstrings)
    ! Use like the following
    !   character(8), dimension(4) :: strs
    !   data strs/'AAAAAAAA', 'BBBBBBBB', 'CCCCCCCC', 'DDDDDDDD'/
    !   call print_strings(strs, 4)
    integer,      intent(in) :: nstrings
    character(*), intent(in) :: strings(nstrings)

    character(c_char), allocatable, target :: cstring(:,:)
    type(c_ptr),       allocatable, target :: ary(:)
    integer                                :: i, j
    integer(c_size_t)                      :: ary_len

    allocate(cstring(LEN(strings(1))+1, SIZE(strings)))
    ary_len = SIZE(strings)
    allocate(ary(ary_len))

    do i = 1, SIZE(strings)
       do j = 1, LEN_TRIM(strings(i))
          cstring(j,i) = strings(i)(j:j)
       end do
       cstring(LEN_TRIM(strings(i))+1, i) = C_NULL_CHAR
       ary(i) = C_LOC(cstring(1,i))
    end do
    call C_print_strings(ary, ary_len)

    deallocate(ary)
    deallocate(cstring)
  end subroutine print_strings
#endif

end module kmrnextf
