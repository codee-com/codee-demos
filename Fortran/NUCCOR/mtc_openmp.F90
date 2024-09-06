module mtc_openmp_module
    use, intrinsic :: iso_fortran_env, only : int64, real64
    use :: mtc_patch_module, only : mtc_patch

    implicit none
    private

    public :: mtc_openmp

    type :: mtc_openmp
    contains
        procedure :: contract => contract
        procedure :: contract_simple => contract_simple
        procedure :: transpose_1_4_3_2_real64 => transpose_1_4_3_2_real64
        procedure :: cleanup => cleanup
        procedure :: clear => clear
    end type mtc_openmp

    interface mtc_openmp
        module procedure constructor
    end interface mtc_openmp

contains
    function constructor() result(this)
        type(mtc_openmp) :: this

        call this%clear()
    end function constructor

    subroutine contract(this, t2, t3, v, f, p)
        class(mtc_openmp), intent(in) :: this
        real(real64), dimension(:,:,:), intent(inout) :: t2
        real(real64), dimension(:,:,:,:), intent(in) :: t3
        real(real64), dimension(:,:,:), intent(in) :: v, f
        type(mtc_patch), intent(in) :: p

        real(real64) :: temp1, temp2
        integer :: ij, bidx, b, a, midx, m, ef

        do ij = 1, p%nij
            do bidx = 1, p%nc
                do a = 1, p%np
                    temp1 = 0.0d0
                    do midx = 1, p%nk
                        m = p%kmap(midx)
                        temp2 = 0.0d0
                        do ef = 1, p%nab
                            temp2 = temp2+ t3(ef, bidx, ij, midx)*v(ef, a, m)
                        end do
                        temp1 = temp1 + temp2*f(midx, a, bidx)
                    end do
                    b = p%cmap(bidx)
                    t2(a, b, ij) = t2(a, b, ij) + temp1
                end do
            end do
        end do
    end subroutine contract

    subroutine transpose_1_4_3_2_real64(this, src, dst)
        class(mtc_openmp), intent(in) :: this
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

    subroutine contract_simple(this, dst, src, op)
        class(mtc_openmp), intent(in) :: this
        real(real64), dimension(:,:,:,:), intent(inout) :: dst
        real(real64), dimension(:,:,:,:,:,:), intent(in) :: src
        real(real64), dimension(:,:,:,:), intent(in) :: op

        integer :: nh, np, i, j, m, a, b, e, f
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

    subroutine cleanup(this)
        class(mtc_openmp), intent(inout) :: this

        call this%clear()
    end subroutine cleanup

    subroutine clear(this)
        class(mtc_openmp), intent(inout) :: this
    end subroutine clear
end module mtc_openmp_module
