SUBROUTINE test_routine01() 
    IMPLICIT NONE
    REAL(KIND=8) :: stopped

    IF (stopped > 0.5) THEN
        WRITE(*,*) 'Loop exited before step n'
    ELSE
        WRITE(*,*) 'No early exit'
    END IF
END SUBROUTINE
