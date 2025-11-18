SUBROUTINE test_routine13pointer_nullify(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    INTEGER(KIND=4) :: i
    INTEGER, POINTER :: sum_pointer

    IF (n > 0) THEN
        NULLIFY(sum_pointer)
    END IF

    DO i = 1, n
        sum_pointer = sum_pointer + 1
    END DO
END SUBROUTINE
