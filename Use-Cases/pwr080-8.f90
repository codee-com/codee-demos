SUBROUTINE test_routine14(n, sum)
    IMPLICIT NONE
    INTEGER(KIND=4) :: n
    REAL(KIND=8), DIMENSION(100) :: A
    REAL(KIND=8) :: sum
    INTEGER(KIND=4) :: i

    IF (n > 0) THEN
        A(1) = 1.0
    END IF

    sum = 0.0
    DO i = 1, 100
        sum = sum + A(i)
    END DO
END SUBROUTINE
