SUBROUTINE test_routine03pointer(n)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    INTEGER(KIND=4) :: i 
    INTEGER, POINTER :: sum_pointer

    DO i = 1, n
        sum_pointer = sum_pointer + i
    END DO
END SUBROUTINE
