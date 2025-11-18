SUBROUTINE test_routine11(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    INTEGER(KIND=4) :: i
    REAL(KIND=4) :: sum

    IF (n > 0) THEN
        sum = 1.0
    ELSE
        !  no initialization
    END IF

    DO i = 1, n
        sum = sum + i
    END DO
END SUBROUTINE
