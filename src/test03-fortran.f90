module test03_constant
  implicit none

  logical(1), parameter :: kPrint = .TRUE.
  integer,    parameter :: kDumpCount = 20

  integer(8), parameter :: kDimension3 = 3
  integer(8), parameter :: kDim3_0 = 10
  integer(8), parameter :: kDim3_1 = 10
  integer(8), parameter :: kDim3_2 = 10
  integer(8), parameter :: kDimension2 = 2
  integer(8), parameter :: kDim2_0 = 10
  integer(8), parameter :: kDim2_1 = 10
end module test03_constant

module test03
  use iso_c_binding
  use kmrnextf
  use test03_constant
  implicit none

  !------ C helper functions
  interface
     type(c_ptr) function C_dumper_helper(key_str, val) &
          bind(c, name='dumper_helper')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: key_str
       integer(c_long),   intent(in), value :: val
     end function C_dumper_helper

     subroutine C_print_data_store_helper(str, ds_count, nspace, print_count) &
          bind(c, name='print_data_store_helper')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: str
       integer(c_long),   intent(in), value :: ds_count
       integer(c_int),    intent(in), value :: nspace
       integer(c_int),    intent(in), value :: print_count
     end subroutine C_print_data_store_helper

     subroutine C_print_get_result_helper(key_req_str, key_ans_str, &
          val, vsize) bind(c, name='print_get_result_helper')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: key_req_str
       type(c_ptr),       intent(in), value :: key_ans_str
       integer(c_long),   intent(in), value :: val
       integer(c_size_t), intent(in), value :: vsize
     end subroutine C_print_get_result_helper

     subroutine C_print_get_view_result_helper1(view_str, key_str, dp_count, &
          print_count) bind(c, name='print_get_view_result_helper1')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: view_str
       type(c_ptr),       intent(in), value :: key_str
       integer(c_size_t), intent(in), value :: dp_count
       integer(c_int),    intent(in), value :: print_count
     end subroutine C_print_get_view_result_helper1

     subroutine C_print_get_view_result_helper2(dp_key_str, dat_val) &
          bind(c, name='print_get_view_result_helper2')
       use iso_c_binding
       implicit none
       type(c_ptr),       intent(in), value :: dp_key_str
       integer(c_long),   intent(in), value :: dat_val
     end subroutine C_print_get_view_result_helper2
  end interface

contains

  integer(c_int) function loader(ds, file_ptr) bind(c) result(zz)
    implicit none
    type(c_ptr), intent(in), value :: ds
    type(c_ptr), intent(in), value :: file_ptr

    integer                    :: ierr, i, j, k
    integer(c_size_t)          :: len_file
    character(c_char), pointer :: file(:)
    type(c_ptr)                :: key, d
    integer(c_long),   target  :: val

    len_file = C_strlen(file_ptr)
    call C_F_POINTER(file_ptr, file, [len_file])

    key = kmrnext_create_key(kDimension3)
    do i = 1, kDim3_0
       ierr = kmrnext_key_set_dim(key, 1, int(i,kind=c_long))
       do j = 1, kDim3_1
          ierr = kmrnext_key_set_dim(key, 2, int(j,kind=c_long))
          do k = 1, kDim3_2
             ierr = kmrnext_key_set_dim(key, 3, int(k,kind=c_long))
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

    type(c_ptr)                :: key, dat, valptr
    type(c_ptr)                :: key_str
    integer(c_long),   pointer :: val

    key = kmrnext_dp_key(dp)
    dat = kmrnext_dp_data(dp)
    key_str = C_kmrnext_key_string(key)
    valptr = kmrnext_data_value(dat)
    call C_F_POINTER(valptr, val)
    zz = C_dumper_helper(key_str, val)
  end function dumper

  integer(c_int) function summarizer(ids, ods, key, dps, env) &
       bind(c) result(zz)
    implicit none
    type(c_ptr),     intent(in), value :: ids
    type(c_ptr),     intent(in), value :: ods
    type(c_ptr),     intent(in), value :: key
    type(datapacks), intent(in), value :: dps
    type(mapenvf),   intent(in), value :: env

    integer(c_long), target :: sum
    integer(c_size_t) :: i
    type(c_ptr), pointer       :: dplst(:)
    type(c_ptr) :: dp, dp_dat, valptr, dat
    integer(c_long), pointer :: val
    integer :: ierr
#ifdef BACKEND_KMR
    integer :: comm_size
#endif

    sum = 0
    call C_F_POINTER(dps%data, dplst, [dps%count])
    do i = 1, dps%count
       dp = dplst(i)
       dp_dat = kmrnext_dp_data(dp)
       valptr = kmrnext_data_value(dp_dat)
       call C_F_POINTER(valptr, val)
       sum = sum + val
    end do
    dat = kmrnext_create_data(C_LOC(sum), int(8,kind=c_long))
    ierr = kmrnext_ds_add(ods, key, dat)
    ierr = kmrnext_free_data(dat)
    zz = 0
#ifdef BACKEND_KMR
    if (env%rank == 0) then
       call mpi_comm_size(env%mpi_comm, comm_size, ierr)
       if (comm_size /= 1) then
          write (*,*) 'MPI_Comm size should be 1.'
       end if
    end if
#endif
  end function summarizer

  subroutine print_data_store(ds, nspace, count, rank)
    implicit none
    type(c_ptr),  intent(in) :: ds
    integer,      intent(in) :: nspace
    integer,      intent(in) :: count
    integer,      intent(in) :: rank

    integer(c_long)            :: ds_count
    type(c_ptr)                :: dumped_str
    !character(c_char), pointer :: dumped_str(:)

    ds_count = kmrnext_ds_count(ds)
    dumped_str = C_kmrnext_ds_dump(ds, C_FUNLOC(dumper))
    !call kmrnext_ds_dump(ds, dumped_str, dumper)
    if (rank /= 0) then
       return
    end if
    call C_print_data_store_helper(dumped_str, ds_count, nspace, count)
  end subroutine print_data_store

  subroutine print_get_result(key, dp, rank)
    implicit none
    type(c_ptr),  intent(in) :: key
    type(c_ptr),  intent(in) :: dp
    integer,      intent(in) :: rank

    type(c_ptr)               :: key_req_str, key_ans_str
    type(c_ptr)               :: key_ans, dat, valptr
    integer(c_long),  pointer :: val
    integer(c_size_t)         :: vsize

    if (rank /= 0) then
       return
    end if

    key_req_str = C_kmrnext_key_string(key)
    key_ans = kmrnext_dp_key(dp)
    key_ans_str = C_kmrnext_key_string(key_ans)
    dat = kmrnext_dp_data(dp)
    valptr = kmrnext_data_value(dat)
    call C_F_POINTER(valptr, val)
    vsize = kmrnext_data_size(dat)
    call C_print_get_result_helper(key_req_str, key_ans_str, val, vsize)
  end subroutine print_get_result

  subroutine print_get_view_result(dps, view, key, count, rank)
    implicit none
    type(datapacks), intent(in), value :: dps
    type(c_ptr),     intent(in)        :: view
    type(c_ptr),     intent(in)        :: key
    integer,         intent(in)        :: count
    integer,         intent(in)        :: rank

    type(c_ptr)                :: view_str, key_str
    type(c_ptr), pointer       :: dplst(:)
    integer(c_size_t)          :: i, cnt
    type(c_ptr)                :: dp, dp_key, dp_dat, dat_valptr, dp_key_str
    integer(c_long),   pointer :: dat_val

    if (rank /= 0) then
       return
    end if

    view_str = C_kmrnext_view_string(view)
    key_str = C_kmrnext_key_string(key)
    call C_print_get_view_result_helper1(view_str, key_str, dps%count, count)

    call C_F_POINTER(dps%data, dplst, [dps%count])
    cnt = 0
    do i = 1, dps%count
       if (count > 0 .and. cnt >= count) exit
       cnt = cnt + 1
       dp = dplst(i)
       dp_key = kmrnext_dp_key(dp)
       dp_key_str = C_kmrnext_key_string(dp_key)
       dp_dat = kmrnext_dp_data(dp)
       dat_valptr = kmrnext_data_value(dp_dat)
       call C_F_POINTER(dat_valptr, dat_val)
       call C_print_get_view_result_helper2(dp_key_str, dat_val)
    end do
    write (*,*) ''
  end subroutine print_get_view_result

end module test03


program main
  use iso_c_binding
  use kmrnextf
  use test03_constant
  use test03
#ifdef BACKEND_KMR
  !  include "mpif.h"
  use mpi
#endif
  implicit none
  !---------------------------------------------------------------------------
  integer                    :: rank = 0
  integer                    :: ierr
  type(c_ptr)                :: next
  type(c_ptr)                :: ds1, ds2
  integer(c_size_t)          :: sizes3(Max_Dimension_Size)
  integer(c_size_t)          :: sizes2(Max_Dimension_Size)
  character(6), dimension(1) :: files
  type(c_ptr)                :: key1, key2
  integer(c_size_t)          :: kval1(Max_Dimension_Size)
  integer(c_size_t)          :: kval2(Max_Dimension_Size)
  type(c_ptr)                :: dp1, dp2
  type(c_ptr)                :: v1, v2, v3, v4
  integer(c_long)            :: flags1(Max_Dimension_Size)
  integer(c_long)            :: flags2(Max_Dimension_Size)
  integer(c_long)            :: flags3(Max_Dimension_Size)
  integer(c_long)            :: flags4(Max_Dimension_Size)
  type(datapacks)            :: dps1, dps2, dps3, dps4, dps5, dps6, dps7, dps8

  character(c_char), pointer :: str(:)

  data sizes3/kDim3_0,kDim3_1,kDim3_2,0,0,0,0,0/
  data sizes2/kDim2_0,kDim2_1,0,0,0,0,0,0/
  data files/'dummy1'/
  data kval1/2,2,2,0,0,0,0,0/
  data kval2/2,2,3,0,0,0,0,0/
  data flags1/Split_All, Split_All, Split_All, &
       Split_None, Split_None, Split_None, Split_None, Split_None/
  data flags2/Split_All, Split_None, Split_All, &
       Split_None, Split_None, Split_None, Split_None, Split_None/
  data flags3/Split_All, Split_None, Split_None, &
       Split_None, Split_None, Split_None, Split_None, Split_None/
  data flags4/Split_None, Split_None, Split_None, &
       Split_None, Split_None, Split_None, Split_None, Split_None/

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
     ! deallocate(str) This will be an error on K
  end if

  !------ Load data contents from a file
  ierr = kmrnext_ds_load_files(ds1, files, 1, loader)
  if (kPrint) then
     if (rank == 0) then
        write (*,*) '1. Load data to a DataStore'
     end if
     call print_data_store(ds1, 2, kDumpCount, rank)
     if (rank == 0) then
        write (*,*) ''
     end if
  end if

  !------ Setup keys
  key1 = kmrnext_create_key(kDimension3)
  key2 = kmrnext_create_key(kDimension3)
  ierr = kmrnext_key_set_size(key1, kval1)
  ierr = kmrnext_key_set_size(key2, kval2)

  !------ Get a data from a DataStore
  dp1 = kmrnext_ds_get(ds1, key1)
  dp2 = kmrnext_ds_get(ds1, key2)
  if (kPrint) then
     if (rank == 0) then
        write (*,*) '2. Get a data from a DataStore by get()'
     end if
     call print_get_result(key1, dp1, rank)
     call print_get_result(key2, dp2, rank)
     if (rank == 0) then
        write (*,*) ''
     end if
  end if
  ierr = kmrnext_free_dp(dp1)
  ierr = kmrnext_free_dp(dp2);

  !------ Setup views
  v1 = kmrnext_create_view(kDimension3)
  v2 = kmrnext_create_view(kDimension3)
  v3 = kmrnext_create_view(kDimension3)
  v4 = kmrnext_create_view(kDimension3)
  ierr = kmrnext_view_set(v1, flags1)
  ierr = kmrnext_view_set(v2, flags2)
  ierr = kmrnext_view_set(v3, flags3)
  ierr = kmrnext_view_set(v4, flags4)

  !------ Get a data from a DataStore with a view
  dps1 = kmrnext_ds_get_view(ds1, key1, v1)
  dps2 = kmrnext_ds_get_view(ds1, key2, v1)
  dps3 = kmrnext_ds_get_view(ds1, key1, v2)
  dps4 = kmrnext_ds_get_view(ds1, key2, v2)
  dps5 = kmrnext_ds_get_view(ds1, key1, v3)
  dps6 = kmrnext_ds_get_view(ds1, key2, v3)
  dps7 = kmrnext_ds_get_view(ds1, key1, v4)
  dps8 = kmrnext_ds_get_view(ds1, key2, v4)
  if (kPrint) then
     if (rank == 0) then
        write (*,*) '3. Get data from a DataStore by get(view)'
     end if
     call print_get_view_result(dps1, v1, key1, kDumpCount, rank)
     call print_get_view_result(dps2, v1, key2, kDumpCount, rank)
     call print_get_view_result(dps3, v2, key1, kDumpCount, rank)
     call print_get_view_result(dps4, v2, key2, kDumpCount, rank)
     call print_get_view_result(dps5, v3, key1, kDumpCount, rank)
     call print_get_view_result(dps6, v3, key2, kDumpCount, rank)
     call print_get_view_result(dps7, v4, key1, kDumpCount, rank)
     call print_get_view_result(dps8, v4, key2, kDumpCount, rank)
     if (rank == 0) then
        write (*,*) ''
     end if
  end if
  ierr = kmrnext_free_datapacks(dps1)
  ierr = kmrnext_free_datapacks(dps2)
  ierr = kmrnext_free_datapacks(dps3)
  ierr = kmrnext_free_datapacks(dps4)
  ierr = kmrnext_free_datapacks(dps5)
  ierr = kmrnext_free_datapacks(dps6)
  ierr = kmrnext_free_datapacks(dps7)
  ierr = kmrnext_free_datapacks(dps8)

  !------ Apply map functions
  ds2 = kmrnext_create_ds(next, kDimension2)
  ierr = kmrnext_ds_set_size(ds2, sizes2)
  ierr = kmrnext_ds_map(ds1, ds2, v2, summarizer)
  if (kPrint) then
     if (rank == 0) then
        write (*,*) '4. Apply map to each data in a DataStore'
        write (*,*) '  Output DataStore'
     end if
     call print_data_store(ds2, 4, kDumpCount, rank)
     if (rank == 0) then
        write (*,*) ''
     end if
  end if
  ierr = kmrnext_free_ds(ds2)

  ierr = kmrnext_free_view(v1)
  ierr = kmrnext_free_view(v2)
  ierr = kmrnext_free_view(v3)
  ierr = kmrnext_free_view(v4)

  ierr = kmrnext_free_key(key1)
  ierr = kmrnext_free_key(key2)

  ierr = kmrnext_free_ds(ds1)
  ierr = kmrnext_finalize()
end program main
