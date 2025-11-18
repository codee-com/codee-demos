SUBROUTINE test_routine12(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    INTEGER(KIND=4) :: i
    INTEGER(KIND=4) :: sum_integer

    IF (n > 0) THEN
        sum_integer = 1
    END IF

    DO i = 1, n
        sum_integer = sum_integer + i
    END DO
END SUBROUTINE
