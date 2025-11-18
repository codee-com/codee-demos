SUBROUTINE test_routine12(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    INTEGER(KIND=4) :: i
    REAL(KIND=4) :: sum_real

    IF (n > 0) THEN
        sum_real    = 1.0
    END IF

    DO i = 1, n
        sum_real    = sum_real + i
    END DO
END SUBROUTINE
