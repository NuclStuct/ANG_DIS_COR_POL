!MODULE FOR Angular Distribution and it's functions
	!REF [1]: NUCLEAR DATA Section A, Volume 3,Number 1, August 1967
	!TABLES OF COEFFICIENTS FOR ANGULAR DISTRIBUTION OF GAMMA RAYS FROM ALIGNED NUCLEI*
	!	BY : T. YAMAZAKI
	
	!REF [2]: NUCLEAR DATA TABLES, 11, 351-406 (1973)
	!DIRECTIONAL CORRELATIONS OF GAMMA RADIATIONS EMITTED FROM NUCLEAR STATES ORIENTED
	!BY NUCLEAR REACTION OR CRYOGENIC METHODS
	!	BY : K. S. KRANE, R. M. STEFFEN and R. M. WHEELER
	
	MODULE ANG_DIS
	USE COM
	USE QM_AM
	USE MATH
	IMPLICIT NONE

	CONTAINS
	
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION Fk(k,Jf,L1,L2,Ji)
	!FUNCTION FOR Fk
	!REF [1]:
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: k
	REAL(8), INTENT (IN) :: Jf,L1,L2,Ji
	REAL(8) :: t1,t2
	
	t1 = (-1)**(Jf-Ji-1)
	t2 = SQRT((2*L1+1)*(2*L2+1)*(2*Ji+1))
	
	Fk = t1*t2*cgc(L1,L2,REAL(k,8),1.0d0,-1.0d0)*racah_W(Ji,Ji,L1,L2,REAL(k,8),Jf)
	
	END FUNCTION Fk
		
!-----------------------------------------------------------------------------------------!	
	REAL(8) FUNCTION Pm(m,J)
	!FUNCTION FOR POPULATION PARAMETER Pm
	!REF [1]:
	IMPLICIT NONE

	REAL(8), INTENT (IN) :: m,J
	REAL(8) :: t1,t2,mm

	t1 = EXP(-1.0D0*((m/(sbyj*J))**2)/2.0D0)
	t2 = 0.0D0
	mm = -1.0D0*abs(J)

	DO
		t2 = t2+EXP(-1.0D0*((mm/(sbyj*J))**2)/2.0D0)
		mm = mm+1
		IF (mm.GT.ABS(J)) EXIT
	END DO
	
	Pm = t1/t2
	
	END FUNCTION Pm

!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION rho(k,J)
	!FUNCTION FOR STATISTICAL TENSOR
	!REF [1]
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: k
	REAL(8), INTENT (IN) :: J
	REAL(8) :: t1,t2,m
	
	m = -1.0d0*abs(J)
	t1 = SQRT(2.0d0*J+1.0D0)
	t2 = 0.0d0
	
	DO
		t2 = t2+(-1)**(J-m)*cgc(J,J,REAL(k,8),m,-m)*Pm(m,J)
		m = m+1
		IF (m.GT.ABS(J)) EXIT
	END DO
	
	rho = t1*t2
	
	END FUNCTION rho

!-----------------------------------------------------------------------------------------!		
	REAL(8) FUNCTION Bk(k,J)
	!FUNCTION FOR STATISTICAL TENSOR Bk FOR COMPLETE ALIGNMENT
	!REF [1]
	IMPLICIT NONE

	INTEGER, INTENT (IN) :: k
	REAL(8), INTENT (IN) :: J
	REAL(8) :: t1

	t1 = SQRT(2.0d0*J+1.0D0)

	IF (MOD(J,1.0d0).EQ.0.0d0) THEN
		Bk = t1*cgc(J,J,REAL(k,8),0.0d0,0.0d0)*(-1)**(J)
	ELSE
		Bk = t1*cgc(J,J,REAL(k,8),0.5d0,-0.5d0)*(-1)**(J-0.5d0)
	END IF

	END FUNCTION Bk

!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION Akmax(k,Ji,L1,L2,Jf,d)
	!FUNCTION FOR Ak WITH MAXIMUM ALIGNMENT
	!REF [1]
	IMPLICIT NONE

	INTEGER, INTENT (IN) :: k
	REAL(8), INTENT (IN) :: Ji,L1,L2,Jf,d
	REAL(8) :: L1L1,L1L2,L2L2
	
	L1L1 = Fk(k,Jf,L1,L1,Ji)
	L1L2 = Fk(k,Jf,L1,L2,Ji)
	L2L2 = Fk(k,Jf,L2,L2,Ji)

	Akmax=(1.0D0/(1.0D0+d*d))*Bk(k,Ji)*(L1L1+2.0D0*d*L1L2+d*d*L2L2)
	
	END FUNCTION Akmax
	
!-----------------------------------------------------------------------------------------!	
	REAL(8) FUNCTION alpha(k,J)
	!FUNCTION FOR ATTENUATION COEFFICIENT OF ALIGNMENT
	!FUNCTION FOR STATISTICAL TENSOR Bk FOR COMPLETE ALIGNMENT
	!REF [1]
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: k
	REAL(8), INTENT (IN) :: J

	alpha = rho(k,J)/Bk(k,J)
	
	END FUNCTION alpha

!-----------------------------------------------------------------------------------------!

	REAL(8) FUNCTION s_uk(k,Ji,L1,Jf)
	!FUNCTION FOR uk (small_uk) TO BE USED IN Uk
	!REF [1]	
	IMPLICIT NONE

	INTEGER, INTENT (IN) :: k
	REAL(8), INTENT (IN) :: Ji,L1,Jf
	REAL(8) :: t1
	
	t1 = SQRT((2.D0*Ji+1.0D0)*(2.D0*Jf+1.0D0))*(-1.0D0**(Ji+Jf-L1))
	
	s_uk = t1*racah_W(Ji,Ji,Jf,Jf,REAL(k,8),L1)
	
	END FUNCTION s_uk
	
!-----------------------------------------------------------------------------------------!	
	REAL(8) FUNCTION Uk(k,Ji,L1,L2,Jf,d)
	!FUNCTION FOR Uk WHICH TRANSLATES STATISTICAL TENSOR FROJ Ji TO JF
	!REF [1]
	IMPLICIT NONE

	INTEGER, INTENT (IN) :: k
	REAL(8), INTENT (IN) :: Ji,L1,L2,Jf,d
	REAL(8) :: u_L1,u_L2

	u_L1 = s_uk(k,Ji,L1,Jf)
	u_L2 = s_uk(k,Ji,L2,Jf)
	
	Uk = (1/(1+d*d))*(u_L1 + d*d*u_L2)
	
	END FUNCTION Uk
	
!-----------------------------------------------------------------------------------------!		
	REAL(8) FUNCTION W(Ji,L1,L2,Jf,d,teta)
	!FUNCTION FOR ANGULAR DISTRIBUTION W(TETA)
	!REF [1]
	IMPLICIT NONE

	REAL(8), INTENT (IN) :: Ji,L1,L2,Jf,d,teta
	REAL(8) :: t1
	INTEGER :: k
	
	t1 = 1.0D0
	
	DO k = 2, kmax, 2

		t1 = t1+alpha(k,Ji)*Akmax(k,Ji,L1,L2,Jf,d)*legP(k,0,DCOSD(teta))
		
	END DO
	
	W = t1
	
	END FUNCTION W
	
!-----------------------------------------------------------------------------------------!		
	REAL(8) FUNCTION W_exp(a2,a4,teta)
	!FUNCTION FOR EXPERIMENTAL ANGULAR DISTRIBUTION W(TETA)
	!REF [1]
	IMPLICIT NONE

	REAL(8), INTENT (IN) :: a2,a4,teta
	
	W_exp = 1.0D0 +a2*legP(2,0,DCOSD(teta))+a4*legP(4,0,DCOSD(teta))
		
	END FUNCTION W_exp
	
!-----------------------------------------------------------------------------------------!		
	REAL(8) FUNCTION dW_exp(a2,a4,da2,da4,teta,dteta)
	!FUNCTION FOR EXPERIMENTAL ERROR ANGULAR DISTRIBUTION W(TETA)
	!BY M ANSER ( CONSIDERING ERRORS IS A2 A4 AND TETA )
	IMPLICIT NONE

	REAL(8), INTENT (IN) :: a2,a4,da2,da4,teta,dteta
	REAL(8) :: t1, t2
		
	t1 = W_exp(da2,da4,teta)-1.0D0 
	t2=DSIND(teta)*(a2*(3.0D0*DCOSD(teta))+a4*0.5D0*(35.0D0*(DCOSD(teta)**3)-15.0D0*DCOSD(teta)))*dteta
		
	dW_exp = t1 - t2
		
	END FUNCTION dW_exp
	
!	t2 = W_exp(da2,da4,tet)-1.0d0 
!	t3=DSIND(teta(j))*(a2*(3.0d0*DCOSD(teta(j)))+a4*0.5d0*(35.0d0*(DCOSD(teta(j))**3)-15.0d0*DCOSD(teta(j))))*6.99676637d0

!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION Qk(TAU, K)
	IMPLICIT NONE
	! Input variables
	REAL(8) :: TAU
	INTEGER :: K

	! Local variables
	REAL(8) :: XJ(2), F(101), YL, YU, DL, XM, EX
	REAL(8) :: EV, OD, FINT
	INTEGER :: N, J, M, KA
	INTEGER :: MAX_DIVISIONS

	! Constants
	PARAMETER (MAX_DIVISIONS = 100)
	! Initialize
	XJ(1) = 0.0D0
	XJ(2) = 0.0D0

	! Loop over N
	DO N = 1, 2
		IF ((K.EQ.0).AND.(N.EQ.2)) CYCLE

		KA = K * (2 - N) + 1

	! Loop over angular ranges J
		DO J = 1, 3
			YL = B(J)
			YU = B(J + 1)
			DL = (YU - YL) / MAX_DIVISIONS

      ! Integration over YL to YU
			EV = 0.0D0
			OD = 0.0D0

			DO M = 1, 101
			XM = YL + DL * (M - 1)
			SELECT CASE (J)
			CASE (1)
			EX = TAU * (A - (D + XL) * TAN(XM)) / SIN(XM)
			CASE (2)
			EX = -TAU * XL / COS(XM)
			CASE (3)
			EX = TAU * (D * TAN(XM) - R) / SIN(XM)
			END SELECT

		        F(M) = SIN(XM) * (1.0D0 - EXP(EX))

!			IF (J.EQ.1) THEN
!			F(M) = F(M) * EXP(-TAU * (A / SIN(XM) - D / COS(XM)))
!			END IF

        ! Apply weighting based on KA
		        SELECT CASE (KA)
		        CASE (2)
		          F(M) = F(M) * COS(XM)
		        CASE (3)
		          F(M) = F(M) * (1.5D0 *(COS(XM)**2) - 0.5D0)
		        CASE (4)
		          F(M) = F(M) * (2.5D0 * (COS(XM)**3) - 1.5D0 * COS(XM))
		        CASE (5)
		          F(M) = F(M) * (4.375D0 * (COS(XM)**4) - 3.75D0 * (COS(XM)**2) + 0.375D0)
		        END SELECT

        ! Accumulate even and odd terms for Simpson's rule
		        IF (MOD(M, 2).EQ.0) THEN
		          EV = EV + F(M)
		        ELSE IF ((M.NE.1).AND.(M.NE.101)) THEN
		          OD = OD + F(M)
		        END IF
		      END DO

      ! Compute the integral using Simpson's rule
	      FINT = (DL / 3.0D0) * (F(1) + 4.0D0 * OD + 2.0D0 * EV + F(101))
	      XJ(N) = XJ(N) + FINT
		END DO
	END DO

  ! Compute final result
	Qk = XJ(1)
	IF (K.NE.0) Qk = Qk / XJ(2)

	RETURN
	END FUNCTION Qk
	
!-----------------------------------------------------------------------------------------!
!FUNCTIONS FOR ANGULAR CORRELATIONS ARE DEFINED BELOW
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION Fk_C(k,k2,k1,L1,L2,Jf,Ji)
	!FUNCTION FOR Fk_C {C FOR CORRELARION}
	!REF [2]
	!INPUT SCHEME = F^{k2 k1}_k (L1,L2,Jf,Ji)
	!		<-----a---||----b------>
	!follow a then b
	
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: k,k1,k2
	REAL(8), INTENT (IN) :: L1,L2,Ji,Jf
	REAL(8) :: t1,t2,t3,t4
	
	t1 = (-1.0D0)**(L2-1)
	t2 = SQRT((2.0D0*L1+1)*(2.0D0*L2+1)*(2.0D0*Ji+1)*(2.0D0*Jf+1))
	t3 = SQRT((2.0D0*REAL(k,8)+1.0D0)*(2.0D0*REAL(k2,8)+1.0D0)*(2.0D0*REAL(k1,8)+1.0D0))
	t4 = t1*t2*t3*three_J(L1,L2,REAL(k,8),1.0D0,-1.0D0,0.0D0)
	
	Fk_C = t4*nine_J(Jf,L1,Ji,Jf,L2,Ji,REAL(k2,8),REAL(k,8),REAL(k1,8))
	
	END FUNCTION Fk_C
	
	
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION Ak_C(k,k2,k1,L1,L2,Jf,Ji,dd)
	!FUNCTION FOR Ak_C {C FOR CORRELARION}
	!REF [2]
	!INPUT SCHEME = A^{k2 k1}_k (L1,L2,Jf,Ji,dd)
	!		<-----a---||----b------>
	!follow a then b
	
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: k,k1,k2
	REAL(8), INTENT (IN) :: L1,L2,Ji,Jf,dd
	REAL(8) :: L1L1,L1L2,L2L2
	
	L1L1 = Fk_C(k,k2,k1,L1,L1,Jf,Ji)
	L1L2 = Fk_C(k,k2,k1,L1,L2,Jf,Ji)
	L2L2 = Fk_C(k,k2,k1,L2,L2,Jf,Ji)
	
	Ak_C = (1/(1+dd*dd))*(L1L1+2*dd*L1L2+dd*dd*L2L2)
	
	END FUNCTION Ak_C	
	
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION g_N1(kk,k1,k,k2,phi)
	!FUNCTION FOR GEOMETRY-N1
	!REF [2]
	!INPUT SCHEME = g^{k1 k k2}_kk (phi)
	!		<-----a---||----b------>
	!follow a then b
	
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: kk,k1,k,k2
	REAL(8), INTENT (IN) :: phi
	REAL(8) :: t1,t2,t3,t4
	INTEGER :: q
	
	t4=0
!	OPEN(UNIT=99, FILE="debug.out", STATUS="UNKNOWN", ACTION="WRITE")

	
	DO q = 0, MIN(k,k2), 2
	t1 = (2.0D0-dDELTA(q,0))*cgc(REAL(k1,8),REAL(k,8),REAL(k2,8),0.0D0,REAL(q,8))
	t2 = SQRT(f(k-q)*f(k+q)*f(k2-q)*(2.0D0*REAL(k,8)+1.0D0)/((2.0D0*REAL(k2,8)+1.0D0)*f(k2+q)))
	t3 = ((-1.0D0)**((k+q)/2))*AA(kk,k2,q)*DCOSD(REAL(q,8)*phi)/((df(k-q))*(df(k+q)))	
	
	t4 = t4 + t1*t2*t3
!	WRITE(99,*) kk, k1, k, k2, phi, q, t1, t2, t3 ,t1*t2*t3
	END DO
!	CLOSE(99)
	g_N1 = t4
	
	
	END FUNCTION g_N1	
	
!-----------------------------------------------------------------------------------------!
!	REAL(8) FUNCTION g_N2(kk,k1,k,k2,phi)
!	!FUNCTION FOR GEOMETRY-N2
!	!REF [2]
!	!INPUT SCHEME = g^{k1 k k2}_kk (phi)
!	!		<-----a---||----b------>
!	!follow a then b
!	
!	IMPLICIT NONE
!	
!	INTEGER, INTENT (IN) :: k,k1,k2,kk
!	REAL(8), INTENT (IN) :: phi
!	REAL(8) :: t1,t2,t3,t4
!	INTEGER :: q
!	
!	t4=0
!	
!	DO q = 0, MIN(k,k2), 2
!	t1 = (2.0D0-dDELTA(q,0))*cgc(REAL(k1,8),REAL(k,8),REAL(k2,8),0.0D0,REAL(q,8))
!	t2 = SQRT(f(k2-q)*f(k2+q)*f(k2-q)*(2.0D0*k+1.0D0)/((2.0D0*k2+1.0D0)*f(k2+q)))
!	t3 = AA(kk,k,q)*DCOSD(q*phi)*(-1.0D0**((k2+q)/2))/((df(k2-q))*(df(k2+q)))	
!	
!	t4 = t4 + t1*t2*t3
!	END DO
!	
!	g_N2 = t4
!	
!	
!	END FUNCTION g_N2	
	
	
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION g_N2(kk,k1,k,k2,phi)
	!FUNCTION FOR GEOMETRY-N2 BASED ON SYMMETRY RELATIONS
	!REF [2]
	!INPUT SCHEME = g^{k1 k k2}_kk (phi)
	!		<-----a---||----b------>
	!follow a then b
	
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: k,k1,k2,kk
	REAL(8), INTENT (IN) :: phi
	REAL(8) :: t1
	

	t1 = SQRT((2.0D0*k+1.0D0)/(2.0D0*k2+1.0D0))
	
	g_N2 = t1*g_N1(kk,k1,k2,k,phi) 
	
	
	END FUNCTION g_N2	
!	
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION H_N1(k1,k,k2,teta,phi)
	!FUNCTION FOR GEOMETRY-N1
	!REF [2]
	!INPUT SCHEME = H_{k1 k k2} (tete, phi)
	!		--------------------->
	
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: k,k1,k2
	REAL(8), INTENT (IN) :: teta, phi
	REAL(8) :: t1
	INTEGER :: kk
	t1 = 0
	
	DO kk= 0 , k2 , 2
	t1 = t1 + g_N1(kk,k1,k,k2,phi)*legP(kk,0,DCOSD(teta))
	END DO
	
	H_N1 = t1 
	
	END FUNCTION H_N1
	
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION H_N2(k1,k,k2,teta,phi)
	!FUNCTION FOR GEOMETRY-N2
	!REF [2]
	!INPUT SCHEME = H_{k1 k k2} (tete, phi)
	!		--------------------->
	
	IMPLICIT NONE
	
	INTEGER, INTENT (IN) :: k,k1,k2
	REAL(8), INTENT (IN) :: teta, phi
	REAL(8) :: t1
	INTEGER :: kk
	t1 = 0
	
	DO kk= 0 , k , 2
	t1 = t1 + g_N2(kk,k1,k,k2,phi)*legP(kk,0,DCOSD(teta))
	END DO
	
	H_N2 = t1 
	
	END FUNCTION H_N2


!-----------------------------------------------------------------------------------------!	
	REAL(8) FUNCTION W_N1(X1,d1,X2,d2,teta,phi)
	!FUNCTION FOR ANGULAR CORRELATION FOR GEOMETRY-N1
	!REF [2]
	!INPUT SCHEME = X1 , X2 = {L1,L2,Jf,Ji}
	!		--------------------->
	
	IMPLICIT NONE
	
	REAL(8), INTENT (IN) :: X1(1:4), X2(1:4)
	REAL(8), INTENT (IN) :: d1, d2, teta, phi
	REAL(8) :: t1,t2
	INTEGER :: k,k1,k2
	
	t1 = 0
	
	DO k = 0 , kmax , 2
		DO k1 = 0 , kmax , 2
			DO k2 = 0 , kmax , 2
		t2 = 0
		t2 = Akmax(k2,X2(4),X2(1),X2(2),X2(3),d2)*H_N1(k1,k,k2,teta,phi)/Bk(k2,X2(4))
		t1 = t1 + rho(k1,X1(4))*Ak_C(k,k2,k1,X1(1),X1(2),X1(3),X1(4),d1)*t2
			END DO
		END DO
	END DO
	
	W_N1 = t1 
	
	END FUNCTION W_N1

!-----------------------------------------------------------------------------------------!	
	REAL(8) FUNCTION W_N2(X1,d1,X2,d2,teta,phi)
	!FUNCTION FOR ANGULAR CORRELATION FOR GEOMETRY-N1
	!REF [2]
	!INPUT SCHEME = X1 , X2 = {L1,L2,Jf,Ji}
	!		--------------------->
	
	IMPLICIT NONE
	
	REAL(8), INTENT (IN) :: X1(1:4), X2(1:4)
	REAL(8), INTENT (IN) :: d1, d2, teta, phi
	REAL(8) :: t1,t2
	INTEGER :: k,k1,k2
	
	t1 = 0
	
	DO k = 0 , kmax , 2
		DO k1 = 0 , kmax , 2
			DO k2 = 0 , kmax , 2
		t2 = 0
		t2 = Akmax(k2,X2(4),X2(1),X2(2),X2(3),d2)*H_N2(k1,k,k2,teta,phi)/Bk(k2,X2(4))
		t1 = t1 + rho(k1,X1(4))*Ak_C(k,k2,k1,X1(1),X1(2),X1(3),X1(4),d1)*t2
			END DO
		END DO
	END DO
	
	W_N2 = t1 
	
	END FUNCTION W_N2

	
	
	END MODULE ANG_DIS
