SUBROUTINE test_routine02(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    INTEGER(KIND=4) :: i
    REAL(KIND=4) :: sum

    DO i = 1, n
        sum = sum + i
    END DO
END SUBROUTINE
