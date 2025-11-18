SUBROUTINE test_routine04(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    REAL(KIND=8), DIMENSION(100) :: A
    REAL(KIND=8) :: sum
    INTEGER(KIND=4) :: i

    sum = 0
    DO i = 1, n
        sum = sum + A(i)
    END DO
END SUBROUTINE

