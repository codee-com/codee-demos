SUBROUTINE test_routine09(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    REAL(KIND=4) :: stopped

    IF (n > 0) THEN
        stopped = 1.0
    END IF

    IF (stopped > 0.5) THEN
        WRITE(*,*) 'Loop exited before step n'
    ELSE
        WRITE(*,*) 'No early exit'
    END IF
END SUBROUTINE
