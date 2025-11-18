SUBROUTINE test_routine08(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    INTEGER(KIND=4) :: i
    REAL(KIND=4) :: sum
    REAL(KIND=8), DIMENSION(100) :: A
    REAL(KIND=8), DIMENSION(2) :: TMP

    TMP = A(2:3)

    DO i = 1, n
        sum = sum + i
    END DO
END SUBROUTINE
