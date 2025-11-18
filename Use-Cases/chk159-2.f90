SUBROUTINE test_routine17(n, sum)
    IMPLICIT NONE
    INTEGER(KIND=4) :: n
    REAL(KIND=8), DIMENSION(100) :: A
    REAL(KIND=8) :: sum
    INTEGER(KIND=4) :: i

    DO i = 4, 23
        A(i) = 1.0
    END DO

    sum = 0.0
    DO i = 1, n
        sum = sum + A(i)
    END DO
END SUBROUTINE
