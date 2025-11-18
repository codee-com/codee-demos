SUBROUTINE test_routine06()
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(100) :: A

    IF (COUNT(A(2:3) > 0) == 0) THEN
        WRITE(*,*) 'Loop exited before step n'
    ELSE
        WRITE(*,*) 'No early exit'
    END IF
END SUBROUTINE
