!MODULE FOR MATHEMATICAL FUNCTIONS
	MODULE MATH
	USE COM
	IMPLICIT NONE
	
	CONTAINS
	
!--------------------------------------------------------------------------------------!	
	SUBROUTINE fact(f)
	!factorial subroutine
	IMPLICIT NONE
	INTEGER :: i
	REAL(8) :: f(0:50)
	
	f(0) = 1.0d0
	
	DO i = 1, 50
		f(i) = f(i-1)*REAL(i)
	END DO

	END SUBROUTINE
	
!--------------------------------------------------------------------------------------!	
	SUBROUTINE dfact(df)
	!double factorial subroutine
	IMPLICIT NONE
	REAL(8) :: df(0:10)
	INTEGER :: n, k

	DO n = 0, 10
	
       	df(n) = 1.0d0
       	
	IF (n > 1) THEN
	
		DO k = n, 1, -2
		df(n) = df(n) * REAL(k,8)
		END DO
		
	END IF
	
	END DO
	
	END SUBROUTINE
	
!--------------------------------------------------------------------------------------!
	SUBROUTINE expCOFF(AAA)
	IMPLICIT NONE
	
	REAL(8) :: AAA(1:3,1:3,1:3)
	
	AAA = 0.0D0
	AAA(1,1,1) = 1.0D0
	AAA(1,2,2) = 2.0D0
	AAA(1,3,2) = 2.0D0
	AAA(1,3,3) = 56.0D0
	AAA(2,2,1) = 1.0D0
	AAA(2,2,2) = -2.0D0
	AAA(2,3,2) = 10.0D0
	AAA(2,3,3) = -80.0D0
	AAA(3,3,1) = 1.0D0
	AAA(3,3,2) = -12.0D0
	AAA(3,3,3) = 24.0D0

	END SUBROUTINE
!--------------------------------------------------------------------------------------!	
    	FUNCTION legP(l, m, x) RESULT(P_lm)
	! Computes the associated Legendre polynomial P_l^m(x)
	!ref-wikipedia 4th recurrence relation (ChatGPT)
	IMPLICIT NONE
        INTEGER, INTENT(IN) :: l, m
        REAL(8), INTENT(IN) :: x
        REAL(8) :: P_lm

        INTEGER :: k
        REAL(8) :: P_l_prev, P_l_curr, P_l_next, factor

	! Check for valid input
        IF (m > l .OR. m < 0 .OR. ABS(x) > 1.0) THEN
        	PRINT *, "Invalid input: l >= m >= 0 and |x| <= 1"
        	STOP
        END IF

	! Handle the special case m = 0 (basic Legendre polynomial P_l)
        IF (m == 0) THEN
        	P_l_prev = 1.0D0
        	P_l_curr = x
        	IF (l == 0) THEN
                	P_lm = P_l_prev
                	RETURN
        	ELSE IF (l == 1) THEN
                	P_lm = P_l_curr
                	RETURN
                END IF

		DO k = 2, l
			P_l_next = ((2 * k - 1) * x * P_l_curr - (k - 1) * P_l_prev) / k
	                P_l_prev = P_l_curr
	                P_l_curr = P_l_next
		END DO
		
		P_lm = P_l_curr
	        RETURN
        END IF
		
	! Handle m > 0 using the recurrence relation
	factor = SQRT(1.0D0 - x * x)
	P_lm = 1.0D0

	DO k = 1, m
		P_lm = -P_lm * factor * (2.0D0 * k - 1.0D0)
	END DO

	IF (l == m) RETURN
	
	P_l_prev = P_lm
	P_l_curr = x * (2 * m + 1) * P_l_prev

	IF (l == m + 1) THEN
		P_lm = P_l_curr
	        RETURN
	END IF

        DO k = m + 2, l
        	P_l_next = ((2 * k - 1) * x * P_l_curr - (k + m - 1) * P_l_prev) / (k - m)
        	P_l_prev = P_l_curr
        	P_l_curr = P_l_next
        END DO

        P_lm = P_l_curr

	END FUNCTION legP	
	
!--------------------------------------------------------------------------------------!	
	REAL(8) FUNCTION TAUU(E)
	IMPLICIT NONE
	
	REAL(8), INTENT (IN) :: E
	REAL(8) :: m, c
	
	INTEGER :: i
	
	IF ((E.GT.1500.0D0).OR.(E.LT.30.0D0)) THEN
        	PRINT *, "Invalid input: 1500 >= E >= 30 "
        	STOP
        END IF
	
	DO i = 1, 8
	IF ((E.GE.EREF(i)).AND.(E.LT.EREF(i+1))) EXIT
	END DO
	m = (TAUREF(i+1)-TAUREF(i))/(EREF(i+1)-EREF(i))
	c = TAUREF(i)- m*EREF(i)
	
	TAUU = m*E + c

	END FUNCTION TAUU
		
		
!--------------------------------------------------------------------------------------!	
	REAL(8) FUNCTION AA(kk,k,q)
	!EXPANSION COEFFICIENTS a_kk(k,q) TO RELATE ASSOCIATED LEG-POL TO LEG-POL
	!REF [2]	
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: kk,k,q
	INTEGER :: i,j,l
	
	i = kk/2 + 1
	j = k/2 + 1
	l = q/2 + 1
	
	AA=AAA(i,j,l)	
	
	END FUNCTION AA
			
!--------------------------------------------------------------------------------------!	
	REAL(8) FUNCTION dDELTA(a,b)
	!DIRAC DELTA FUNCTION
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: a,b
	
	IF(a.EQ.b) THEN
	dDELTA = 1.0D0
	ELSE
	dDELTA = 0.0D0
	END IF

	
	END FUNCTION dDELTA	
	
		
	
	END MODULE MATH
