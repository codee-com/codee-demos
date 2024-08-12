PROGRAM pi_main
IMPLICIT NONE
    INTEGER :: niters
    DOUBLE PRECISION :: pi
    DOUBLE PRECISION :: tarray(2)
    DOUBLE PRECISION :: time_start,time_finish    
    INTEGER :: nb_ticks_sec, nb_ticks_max 
    INTEGER :: nb_ticks_initial, nb_ticks_final 
    INTEGER :: nb_ticks
    DOUBLE PRECISION :: elapsed_time
    
    niters = 900000000

    CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)

    CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)
    CALL calculate_pi(pi,niters)
    CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)


    nb_ticks = nb_ticks_final - nb_ticks_initial
    IF (nb_ticks_final < nb_ticks_initial) &
            nb_ticks = nb_ticks + nb_ticks_max
    elapsed_time   = REAL(nb_ticks) / ( nb_ticks_sec / 1000 )

            
    PRINT*,"Iteration count = ", niters
    PRINT*,"PI = ",pi
    PRINT*, "time (ms) = ", elapsed_time

END
