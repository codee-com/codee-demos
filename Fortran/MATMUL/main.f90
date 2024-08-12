module checksum_module
contains
   subroutine checksum(C, chk)
      double precision, dimension(:, :), intent(in) :: C
      double precision, intent(out) :: chk
      integer :: i, j, N

      chk = 0.0
      do i = 1, size(C, 1)
         do j = 1, size(C, 2)
            chk = chk + C(i, j)
         end do
      end do
   end
end module checksum_module

program matmul_main

   use :: matmul_module
   use :: checksum_module
   use :: time_module

   implicit none
   integer, parameter :: DEBUG = 0
   integer ::iters, param_iters, N
   double precision, dimension(:, :), allocatable :: A, B, C
   integer :: tstart, tend
   double precision :: elapsed_time
   character(100) :: num
   integer :: i, j, k
   double precision :: chk

   ! Check number of inputs
   if (COMMAND_ARGUMENT_COUNT() == 1) then
      ! Read in the two values
      call GET_COMMAND_ARGUMENT(1, num)
      ! Convert them to integers
      read (num, *) N
   else
      print *, "Error: the command usage is: ./a.out SIZE"
      call EXIT(1)
   end if

   ! Initilizes data
   param_iters = 1
   allocate (A(N, N))
   allocate (B(N, N))
   allocate (C(N, N))

   ! Fill the matrix with incremental values
   do i = 1, N
      do j = 1, N
         A(i, j) = ((i - 1)*N) + (j - 1)
         B(i, j) = ((i - 1)*N) + (j - 1)
      end do
   end do

   ! calulate MATMUML
   call time_inizialization()
   call get_time(tstart)
   do iters = 0, param_iters
      call calculate_matmul(A, B, C)
   end do
   call get_time(tend)

   call get_real_time(tstart, tend, elapsed_time)

   ! print final output
   call checksum(C, chk)
   print *, "time (s)= ", elapsed_time
   print *, "size    = ", N
   if (param_iters > 1) then
      print *, "iters    = ", param_iters
   end if
   print *, "Chksum  = ", chk

   ! Debug Matmul
   if (DEBUG == 1) then
      C = MATMUL(A, B)
      call checksum(C, chk)
      print *, "ChksumOri = ", chk
   end if

   ! free memory
   deallocate (A)
   deallocate (B)
   deallocate (C)

end program matmul_main
