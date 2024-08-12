SUBROUTINE calculate_pi(pi, n)
   IMPLICIT NONE
   DOUBLE PRECISION :: x,sum,pi
   INTEGER :: i,n

   sum=0.0

   DO i= 0, n-1
      x = (i + 0.5) / n
      sum = sum + SQRT(1.0 - x * x) 
   ENDDO

   pi = 4.0 / n * sum
END SUBROUTINE
