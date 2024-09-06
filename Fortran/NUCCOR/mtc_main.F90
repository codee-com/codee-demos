module config
    use, intrinsic :: iso_fortran_env, only : real64
    use :: mtc_patch_module, only : mtc_patch

    implicit none
    private

    public :: create_patch
contains
    type(mtc_patch) function create_patch(number_of_particles, basis_size, nonzero_fraction)
        integer, intent(in) :: number_of_particles, basis_size
        real(real64), intent(in) :: nonzero_fraction

        integer :: nh, np, nc, nab, nk, nij, idx
        integer, dimension(:), allocatable :: cmap, kmap
        logical, dimension(:), allocatable :: mask

        nh = number_of_particles
        np = basis_size - number_of_particles

        nc = max(1, ceiling(nonzero_fraction*np))
        nab = nc**2
        nk = max(1, ceiling(nonzero_fraction*nh))
        nij = nk**2

        mask = get_mask(nk, nh)
        kmap = pack([(idx, idx = 1, nh)], mask)

        mask = get_mask(nc, np)
        cmap = pack([(idx, idx = 1, np)], mask)
        create_patch = mtc_patch(nh, np, nab, nc, nij, nk, cmap, kmap)
    end function create_patch

    function get_mask(nsmall, nlarge) result(mask)
        integer, intent(in) :: nsmall, nlarge
        logical, dimension(:), allocatable :: mask

        integer :: counter, idx
        real :: rnd
        logical :: more

        allocate(mask(nlarge))
        mask = .false.

        more = .true.
        counter = 0
        do while ( more )
            call random_number(rnd)
            idx = ceiling(rnd*nlarge)
            if ( mask(idx) ) cycle

            mask(idx) = .true.
            counter = counter + 1
            more = counter < nsmall
        end do
    end function get_mask
end module config

program mtc_main
    use, intrinsic :: iso_fortran_env, only : real64
    use :: config, only : create_patch
    use :: mtc_patch_module, only : mtc_patch
    use :: mtc_module, only : mtc
    use :: mtc_openmp_module, only : mtc_openmp

    use, intrinsic :: omp_lib, only : omp_get_wtime

    implicit none

    real(real64), dimension(:,:,:,:,:,:), allocatable :: t3, t3_full
    real(real64), dimension(:,:,:), allocatable :: v, t2, f, t2_original
    real(real64), dimension(:,:,:,:), allocatable :: t3_tr_original, t3_tr, v_full, t2_full, t2_full_original
    type(mtc_patch) :: p
    type(mtc) :: anmtc
    type(mtc_openmp) :: anmtc_omp
    integer :: bra, ket, number_of_particles, basis_size, number_of_blocks
    real(real64) :: nonzero_fraction
    logical :: run, run_contract, run_transpose, run_contract_simple
    real(real64) :: memory
    real(real64) :: tstart, tend
    real(real64), parameter :: tolerance = 1.0d-8
    character(len=100) :: dummy


    if ( command_argument_count() < 4 ) then
        write(*,*) "Usage: ./mtc.x nparticles basissize nblocks nonzero_fraction [run kernel]"
        stop
    end if

    run = .true.

    call get_command_argument(1, dummy)
    read(dummy, *) number_of_particles
    call get_command_argument(2, dummy)
    read(dummy, *) basis_size
    call get_command_argument(3, dummy)
    read(dummy, *) number_of_blocks
    call get_command_argument(4, dummy)
    read(dummy, *) nonzero_fraction

    run = .true.
    if ( command_argument_count() > 4 ) then
        call get_command_argument(5, dummy)
        if ( dummy == "no" ) run = .false.
    end if

    run_contract = .false.
    run_contract_simple = .false.
    run_transpose = .false.
    if ( command_argument_count() > 5 ) then
        call get_command_argument(6, dummy)
        if ( dummy == "contract" ) run_contract = .true.
        if ( dummy == "contract_simple" ) run_contract_simple = .true.
        if ( dummy == "transpose" ) run_transpose = .true.
    else
        run_contract = .true.
    end if

    p = create_patch(number_of_particles, basis_size, nonzero_fraction)
    call p%dump(.false.)

    if ( run_contract ) then
        memory = real(number_of_blocks**2, real64)*p%nab*p%nc*p%nij*p%nk
        memory = memory + real(p%nab, real64)*p%np*p%nh
        memory = memory + p%np*p%np*p%nij
        memory = memory + p%nc*p%nk*p%np
        memory = memory + p%np*p%nab*p%nk + real(p%nab, real64)*p%nk*p%nij*p%nc
        memory = memory + p%np*p%nij*p%nc
        memory = memory*8
        write(*,'(a,f10.3,a)') "Memory usage for contract: ", memory/1024/1024/1024, " Gb"
    end if

    if ( run_transpose ) then
        memory = real(number_of_blocks**2, real64)*p%nab*p%nc*p%nij*p%nk
        memory = memory + real(p%nab, real64)*p%np*p%nh*p%nij
        memory = memory + 2*real(p%nab, real64)*p%np*p%nh*p%nij
        memory = memory*8
        write(*,'(a,f10.3,a)') "Memory usage for transpose: ", memory/1024/1024/1024, " Gb"
    end if

    if ( run_contract_simple ) then
        memory = real(p%np, real64)**3*p%nh**3
        memory = memory + real(p%np, real64)**3*p%nh
        memory = memory+ real(p%np, real64)**2*p%nh**2
        memory = memory*8
        write(*,'(a,f10.3,a)') "Memory usage for simple contract: ", memory/1024/1024/1024, " Gb"
    end if

    if (run) then
        if ( run_contract ) then
            allocate(t3(p%nab, p%nc, p%nij, p%nk, number_of_blocks, number_of_blocks))
            allocate(v(p%nab, p%np, p%nh))
            allocate(t2(p%np, p%np, p%nij))
            allocate(t2_original(p%np, p%np, p%nij))
            allocate(f(p%nk, p%np, p%nc))
            call random_number(t3)
            call random_number(v)
            call random_number(f)
            t2_original = 0.0d0
            t2 = 0.0d0

            tstart = omp_get_wtime()
            do ket = 1, number_of_blocks
                do bra = 1, number_of_blocks
                    call anmtc%contract(t2_original, t3(:,:,:,:,bra, ket), v, f, p)
                end do
            end do
            tend = omp_get_wtime()
            write(*,*) "Time spent in contraction reference: ", tend - tstart

            tstart = omp_get_wtime()
            do ket = 1, number_of_blocks
                do bra = 1, number_of_blocks
                    call anmtc_omp%contract(t2, t3(:,:,:,:,bra, ket), v, f, p)
                end do
            end do
            tend = omp_get_wtime()
            write(*,*) "Time spent in contraction Openmp version: ", tend - tstart

            if (.not. all(abs(t2 - t2_original) <= tolerance ) ) then
                 write(*,*) "Test Openmp contraction: ERROR"
            else
                write(*,*) "Test Openmp contraction: OK"
            end if

        end if

        if ( run_transpose ) then
            allocate(t3(p%nab, p%nc, p%nij, p%nk, number_of_blocks, number_of_blocks))
            allocate(t3_tr_original(p%nab, p%nk, p%nij, p%nc))
            allocate(t3_tr(p%nab, p%nk, p%nij, p%nc))

            call random_number(t3)
            tstart = omp_get_wtime()
            do ket = 1, number_of_blocks
                do bra = 1, number_of_blocks
                    call anmtc%transpose_1_4_3_2_real64(t3(:,:,:,:,bra, ket), t3_tr_original)
                end do
            end do
            tend = omp_get_wtime()
            write(*,*) "Time spent in transposed reference: ", tend - tstart

            tstart = omp_get_wtime()
            do ket = 1, number_of_blocks
                do bra = 1, number_of_blocks
                    call anmtc_omp%transpose_1_4_3_2_real64(t3(:,:,:,:,bra, ket), t3_tr)
                end do
            end do
            tend = omp_get_wtime()
            write(*,*) "Time spent in transposed OpenMP version: ", tend - tstart

            if (.not. all(abs(t3_tr - t3_tr_original) <= tolerance ) ) then
                 write(*,*) "Test transposed: ERROR"
            else
                write(*,*) "Test transposed: OK"
            end if

            deallocate(t3, t3_tr_original, t3_tr)
        end if

        if ( run_contract_simple ) then
            allocate(t3_full(p%np, p%np, p%np, p%nh, p%nh, p%nh))
            allocate(v_full(p%nh, p%np, p%np, p%np))
            allocate(t2_full_original(p%np, p%np, p%nh, p%nh))
            allocate(t2_full(p%np, p%np, p%nh, p%nh))
            call random_number(t3_full)
            call random_number(v_full)
            call random_number(t2_full_original)
            t2_full = t2_full_original

            tstart = omp_get_wtime()
            call anmtc%contract_simple(t2_full_original, t3_full, v_full)
            tend = omp_get_wtime()
            write(*,*) "Time spent in simple contraction reference: ", tend - tstart

            tstart = omp_get_wtime()
            call anmtc_omp%contract_simple(t2_full, t3_full, v_full)
            tend = omp_get_wtime()
            write(*,*) "Time spent in contraction Openmp version: ", tend - tstart

            if (.not. all(abs(t2_full - t2_full_original) <= tolerance ) ) then
                 write(*,*) "Test simple contraction: ERROR"
            else
                write(*,*) "Test simple contraction: OK"
            end if

            deallocate(t3_full, v_full, t2_full, t2_full_original)
        end if

    end if
end program mtc_main

