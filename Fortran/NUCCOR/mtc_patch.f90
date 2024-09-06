module mtc_patch_module
    implicit none
    private

    public :: mtc_patch

    type :: mtc_patch
        integer :: nh, np, nab, nc, nij, nk
        integer, dimension(:), allocatable :: cmap, kmap
    contains
        procedure :: dump => dump
        procedure :: cleanup => cleanup
        procedure :: clear => clear
    end type mtc_patch

contains

    subroutine dump(this, dump_map)
        class(mtc_patch), intent(in) :: this
        logical, intent(in) :: dump_map

        integer :: idx

        write(*,*) "nh: ", this%nh
        write(*,*) "np: ", this%np
        write(*,*) "nab: ", this%nab
        write(*,*) "nc: ", this%nc
        write(*,*) "nij: ", this%nij
        write(*,*) "nk: ", this%nk
        write(*,*) "Allocated cmap: ", allocated(this%cmap)
        if ( allocated(this%cmap) .and. dump_map) then
            write(*,*) "Dump of cmap:"
            do idx = 1, size(this%cmap)
                write(*,*) idx, this%cmap(idx)
            end do
        end if
        write(*,*) "Allocated kmap: ", allocated(this%kmap)
        if ( allocated(this%kmap) .and. dump_map) then
            write(*,*) "Dump of kmap:"
            do idx = 1, size(this%kmap)
                write(*,*) idx, this%kmap(idx)
            end do
        end if
    end subroutine dump

    subroutine cleanup(this)
        class(mtc_patch), intent(inout) :: this

        if (allocated(this%cmap)) deallocate(this%cmap)
        if (allocated(this%kmap)) deallocate(this%kmap)

        call this%clear()
    end subroutine cleanup

    subroutine clear(this)
        class(mtc_patch), intent(inout) :: this

        this%np = 0
        this%nh = 0
        this%nab = 0
        this%nc = 0
        this%nij = 0
        this%nh = 0
    end subroutine clear
end module mtc_patch_module

