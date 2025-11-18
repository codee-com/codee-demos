SUBROUTINE test_routine12(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    INTEGER(KIND=4) :: i
    LOGICAL :: sum_logical

    IF (n > 0) THEN
        sum_logical = .true.
    END IF

    DO i = 1, n
        sum_logical = sum_logical .OR. .false.
    END DO
END SUBROUTINE
