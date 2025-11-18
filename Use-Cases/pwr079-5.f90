SUBROUTINE test_routine05()
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(100) :: A

    IF (A(3) > 0.5) THEN
        WRITE(*,*) 'Loop exited before step n'
    ELSE
        WRITE(*,*) 'No early exit'
    END IF
END SUBROUTINE
