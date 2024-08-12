subroutine contract(t2, t3, v, f, np, nab, nc, nij, nk, cmap, kmap)
    use, intrinsic :: iso_fortran_env, only : int64, real64

    real(real64), dimension(:,:,:), intent(inout) :: t2
    real(real64), dimension(:,:,:,:), intent(in) :: t3
    real(real64), dimension(:,:,:), intent(in) :: v, f

    ! These were originally part of the mtc_patch class
    integer, intent(in) :: np, nab, nc, nij, nk        
    integer, dimension(:), allocatable, intent(in) :: cmap, kmap

    real(real64) :: temp1, temp2
    integer :: ij, bidx, b, a, midx, m, ef

#ifdef OPENMP
!$omp parallel do private(b,m,temp1,temp2)
#endif
    do ij = 1, nij
        do bidx = 1, nc
            b = cmap(bidx)
            do a = 1, np
                temp1 = 0.0d0
                do midx = 1, nk
                    m = kmap(midx)
                    temp2 = 0.0d0
                    do ef = 1, nab
                        temp2 = temp2+ t3(ef, bidx, ij, midx)*v(ef, a, m)
                    end do
                    temp1 = temp1 + temp2*f(midx, a, bidx)
                end do
                t2(a, b, ij) = t2(a, b, ij) + temp1
            end do
        end do
    end do
end subroutine contract