
subroutine transpose_1_4_3_2_real64(src, dst)
    use, intrinsic :: iso_fortran_env, only : int64, real64

    real(real64), dimension(:,:,:,:), intent(in) :: src
    real(real64), dimension(:,:,:,:), intent(inout) :: dst

    integer(int64) :: d1, d2, d3, d4

    do d4 = 1, size(src, 4, kind=int64)
        do d3 = 1, size(src, 3, kind=int64)
            do d2 = 1, size(src, 2, kind=int64)
                do d1 = 1, size(src, 1, kind=int64)
                    dst(d1, d4, d3, d2) = src(d1, d2, d3, d4)
                end do
            end do
        end do
    end do
end subroutine transpose_1_4_3_2_real64