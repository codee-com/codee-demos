
subroutine contract_simple(dst, src, op)
    use, intrinsic :: iso_fortran_env, only : int64, real64

    real(real64), dimension(:,:,:,:), intent(inout) :: dst
    real(real64), dimension(:,:,:,:,:,:), intent(in) :: src
    real(real64), dimension(:,:,:,:), intent(in) :: op

    integer :: i, j, m, a, b, e, f
    real(real64) :: temp

    do j = 1, size(dst, 4)
        do i = 1, size(dst, 3)
            do b = 1, size(dst, 2)
                do a = 1, size(dst, 1)
                    temp = 0.0d0
                    do m = 1, size(op, 1)
                        do f = 1, size(op, 3)
                            do e = 1, size(op, 4)
                                temp = temp + op(m, a, e, f)*src(e, f, b, i, j, m)
                            end do
                        end do
                    end do
                    dst(a, b, i, j) = 0.5d0*temp
                end do
            end do
        end do
    end do
end subroutine contract_simple