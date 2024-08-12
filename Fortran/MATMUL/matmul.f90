module matmul_module

   implicit none

contains
   subroutine calculate_matmul(A, B, C)
      double precision, dimension(:, :), intent(in) :: A, B
      double precision, dimension(:, :), intent(out) :: C
      INTEGER :: i, j, k
      ! Initialization
      do i = 1, size(C, 1)
         do j = 1, size(C, 2)
            C(i, j) = 0.0
         end do
      end do
      ! Accumulation
      do i = 1, size(C, 1)
         do j = 1, size(C, 2)
            do k = 1, size(A, 1)
               C(i, j) = C(i, j) + (A(i, k)*B(k, j))
            end do
         end do
      end do
   end subroutine calculate_matmul

end module matmul_module
