!MODULE FOR ANGULAR MOMENTUM ALGEBRA IN QUANTUM MECHANICS
	MODULE QM_AM
	USE COM
	IMPLICIT NONE
	
	CONTAINS
	
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION cgc(j1,j2,j,m1,m2)
	!FUNCTION FOR CLEBSCH-GORDAN COEFFICENTS	
	!<j1j2m1m2|j1j2jm> = [j1 j2 j]	(m1+m2=m)
	!		     [m1 m2 m] =INPUT SCHEME
	!REF: ANGULAR MOMENTUM TECHNIQUE IN QUANTUM MECHANICS 
	!	BY V DEVANATHAN, PAGE:13, EQUATION (2.21)
	IMPLICIT NONE
	
	REAL(8), INTENT (IN) :: j1,j2,j,m1,m2
	REAL(8) :: t1,t2,t3,m
	INTEGER :: a,b,c,d,e,imin,i

	m = m1+m2
	a = (j1+j2-j)
	b = (j1-m1)
	c = (j2+m2)
	d = (j-j2+m1)
	e = (j-j1-m2)
	
	t1 = SQRT((2*j+1)*f(INT(j+j1-j2))*f(INT(j-j1+j2))*f(INT(j1+j2-j))/f(INT(j+j1+j2+1)))
	t2 = SQRT(f(INT(j1+m1))*f(INT(j1-m1))*f(INT(j2+m2))*f(INT(j2-m2))*f(INT(j-m))*f(INT(j+m)))
	imin = min(d,e)

	IF (imin.LT.0) THEN 
		imin=IABS(imin)
	ELSE
		imin=0
	END IF
	
	t3=0.0
	
	DO i = imin ,MIN(a,b,c)
		t3 = t3+((-1)**i)/REAL(f(i)*f(a-i)*f(b-i)*f(c-i)*f(d+i)*f(e+i),8)
	END DO
	
	cgc=t1*t2*t3
	
	END FUNCTION cgc
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION three_J(j1,j2,j,m1,m2,m)
	!Date: 7 Aug 2025
	!FUNCTION FOR Wigner 3-j symbols	
	!(j1 j2 j)	(m1+m2=-m)
	!(m1 m2 m) =INPUT SCHEME
	!REF: ANGULAR MOMENTUM TECHNIQUE IN QUANTUM MECHANICS 
	!	BY V DEVANATHAN, PAGE:16, EQUATION (2.33)
	IMPLICIT NONE
	REAL(8), INTENT (IN) :: j1,j2,j,m1,m2,m
	REAL(8) :: t1,t2,t3
	
	t3 = -m1 - m2 
	
	IF (t3.EQ.m) THEN
	t1 = (-1.0D0)**(j1-j2-m)
	t2 = SQRT(2.0D0*j + 1.0D0)
	three_J = t1*cgc(j1,j2,j,m1,m2)/t2
	ELSE
	three_J = 0.0D0
	END IF
	
	END FUNCTION three_J
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION tri(a,b,c)
	!FUNCTION FOR TRIANGLE CONDITION AND OPERATION
	!USED IN RACAH-W, 6-j AND 9-j 
	!REF: ANGULAR MOMENTUM TECHNIQUE IN QUANTUM MECHANICS 
	!	BY V DEVANATHAN, PAGE:82
	IMPLICIT NONE
	
	REAL(8), INTENT (IN) :: a,b,c
	
!	IF ((c.LT.abs(a-b)).OR.(c.GT.(a+b))) THEN
!		PRINT *, "Invalid input: |a-b|<=c<=(a+b)"
!        	STOP
!	END IF
!ABOVE BLOCK IS GIVING ERROR NEEDS PROPER THINKING
	tri=SQRT(f(int(a+b-c))*f(int(a-b+c))*f(int(b+c-a))/f(int(a+b+c+1)))

	END FUNCTION tri

!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION racah_W(a,b,c,d,e,ff)
	!FUNCTION FOR RACAH-W 
	!W(abcd,ef) = INPUT SCHEME
	!REF: ANGULAR MOMENTUM TECHNIQUE IN QUANTUM MECHANICS 
	!	BY V DEVANATHAN, PAGE:81, EQUATION (7.26)
	IMPLICIT NONE
	
	REAL(8), INTENT (IN) :: a,b,c,d,e,ff
	INTEGER :: aa,bb,cc,dd,ee,fff,gg,imin,i
	REAL(8) :: t1,t2,tri1,tri2,tri3,tri4
	
	aa = INT(-a-b-e)
	bb = INT(-c-d-e)
	cc = INT(-a-c-ff)
	dd = INT(-b-d-ff)
	ee = INT(a+b+c+d)
	fff = INT(a+d+e+ff)
	gg = INT(b+c+e+ff)
	
	tri1 = tri(a,b,e)
	tri2 = tri(c,d,e)
	tri3 = tri(a,c,ff)
	tri4 = tri(b,d,ff)
	
	t1= tri1*tri2*tri3*tri4
	imin=MIN(aa,bb,cc,dd)
	
	IF (imin.LT.0) THEN
		imin = IABS(imin)
	ELSE
		imin = 0
	END IF

 	DO i=imin,MIN(ee,fff,gg)
		t2 = t2+((-1)**(ee+i))*f(i+1)/(f(aa+i)*f(bb+i)*f(cc+i)*f(dd+i)*f(ee-i)*f(fff-i)*f(gg-i))
	END DO

	racah_W = t1*t2

	END FUNCTION racah_W
	
!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION nine_J(a,b,c,d,e,ff,g,h,j)
	!FUNCTION FOR 9-j SYMBOLS
	!REF: QUANTUM THEORY OF ANGULAR MOMENTUM 
	!	BY D A Varshalovich et al., PAGE:340, EQUATION (21)
	!|a b c|
	!{d e ff} = INPUT SCHEME 
	!|g h j|
	IMPLICIT NONE
	
	REAL(8), INTENT (IN) :: a,b,c,d,e,ff,g,h,j
	REAL(8) :: r_W1,r_W2,r_W3,t1
	INTEGER :: i,imin,imax,aep,aem,chp,chm,fgp,fgm
	
	aep = NINT(a + e)
	aem = NINT(a - e)
	chp = NINT(c + h)
	chm = NINT(c - h)
	fgp = NINT(ff + g)
	fgm = NINT(ff - g)
	
	
	t1 = 0.0D0
	imin = MAX(IABS(aem),IABS(chm),IABS(fgm))
	imax = MIN(IABS(aep),IABS(chp),IABS(fgp))
	
	DO i = imin, imax
		r_W1 = racah_W(a,e,c,h,REAL(i,8),b)
		r_W2 = racah_W(a,e,g,ff,REAL(i,8),d)
		r_W3 = racah_W(c,h,ff,g,REAL(i,8),j)
	
		t1 = t1 + (2.0D0*REAL(i,8)+1.0D0)*r_W1*r_W2*r_W3
	END DO
	
	nine_J = t1
	
	END FUNCTION nine_J

!-----------------------------------------------------------------------------------------!
	REAL(8) FUNCTION djmm(j,m1,m2,b)
	!FUNCTION FOR WIGNER D-FUNCTION djmm'(b)
	!REF: QUANTUM THEORY OF ANGULAR MOMENTUM 
	!	BY D A Varshalovich et al., PAGE:77, EQUATION (5)
	IMPLICIT NONE
	
	REAL(8), INTENT (IN) :: j,m1,m2,b
	INTEGER :: i,imax
	REAL(8) :: t1,t2,t3,t4,t5
	
	!CHECKS IF THE INPUTS ARE VALID
	IF ((m1.GT.j).OR.(m1.GT.j).OR.(j.LT.0)) THEN
		PRINT *, "Invalid input: j>=m1,m2 OR j<0"
        	STOP
        ELSE IF ((m1.LT.(-1*j)).OR.(m1.LT.(-1*j))) THEN
		PRINT *, "Invalid input: -j<=m1,m2"
        	STOP        
	END IF
		
	imax= INT(MIN((j-m1),(j+m2)))
	
	t5 = 0.0D0
	t1 = SQRT(f(INT(j+m1))*f(INT(j-m1))*f(INT(j+m2))*f(INT(j-m2))) 

	DO i = MAX(0,INT(m2-m1)), imax
		t2 = ((-1.0d0)**(i))*(DCOSD(b/2.0D0)**(2.0D0*J-2.0D0*REAL(i,8)-m1+m2))
		t3 = DSIND(b/2.0D0)**(2.0D0*REAL(i,8)+m1-m2)
		t4 = f(i)*f(INT(j-m1-i))*f(INT(j+m2-i))*f(INT(m1-m2+i))
		t5 = t5 + t2*t3/t4 
	END DO
	
	djmm = ((-1.0d0)**(m1-m2))*t1*t5
	
	END FUNCTION djmm
	
	
	END MODULE QM_AM
