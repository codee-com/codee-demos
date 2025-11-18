SUBROUTINE test_routine07()
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(100) :: A

    IF (COUNT(A > 0) == 0) THEN
        WRITE(*,*) 'Loop exited before step n'
    ELSE
        WRITE(*,*) 'No early exit'
    END IF
END SUBROUTINE
