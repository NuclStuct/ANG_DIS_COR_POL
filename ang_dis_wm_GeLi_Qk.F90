	program ang_dist
	USE COM
	USE MATH
	USE QM_AM
	USE ANG_DIS
	implicit none

	character*15::f2,f3
	real*8::Jf,L1,L2,Ji,dd,a2,a4,chisq,wcalc,wexp,da2,da4,t1,t2,t3
	integer::i,j,k,l,n
	real(8)::dtet,teta,en
	real(8)::r1,r2,tet,q0,q2,q4,m,x,EX,r3,phi,EY
	REAL(KIND=8), PARAMETER :: PI = 3.141592653589793D0

	call fact(f)

	EREF = (/ 30.0D0, 50.0D0, 100.0D0, 150.0D0, 200.0D0, 300.0D0, 500.0D0, 1000.0D0, 1500.0D0 /)	
	TAUREF = (/ 70.6D0, 17.0D0, 2.59D0, 0.982D0, 0.521D0, 0.196D0, 0.0517D0, 0.00927D0, 0.00382D0 /)

	
	A=0.37D0
	D=5.0D0
	R=1.77D0
	XL=3.19D0
	!initializing angle limits
	B(1) = ATAN2(A, D + XL)
	B(2) = ATAN2(A, D)
	B(3) = ATAN2(R, D + XL)
	B(4) = ATAN2(R, D)
	
	WRITE(*,*) Qk(TAUU(60.0d0),2),Qk(TAUU(60.0d0),4),TAUU(60.0d0) 
	WRITE(*,*) Qk(70.6D0,2),Qk(70.6D0,4)
!i could add it to com	
!	teta = (/ 90.0d0, 115.0d0, 125.0d0, 135.0d0, 145.0d0 /)
	dtet = 	6.99676637d0
!	dtet = 0.0d0	
	do
	
	write(*,*) "Input file name (enter 'exit' if you want to exit)"
	read(*,*) f2,f3
	if (f2.EQ."exit") exit !could make a module out of it
	OPEN(UNIT=14,FILE=f2,STATUS='new')
	OPEN(UNIT=15,FILE=f3,STATUS='new')
	en=30.d0
	do n=1,14699
		
	teta = 0.0d0
!	t2=0.0d0
!	t3=0.0d0
	do j=0,1800
	
	t1=0.0d0
	q0=0.0d0
	q2=0.0d0
	q4=0.0d0
	
	
	do i=1, 10000
	CALL RANDOM_NUMBER(r1)
	CALL RANDOM_NUMBER(r2)
	CALL RANDOM_NUMBER(r3)
	
	tet = (B(4)-B(1))*r1 + B(1)
	phi = (B(4)-B(1))*r3 + B(1)
	if (r2.GE.0.5D0) then
	m=-1.0d0
	else
	m=1.0d0
	end if
	
	if ((tet.GE.B(4)).AND.(tet.LT.B(3))) then
	l=3
	else if ((tet.GE.B(3)).AND.(tet.LT.B(2))) then
	l=2
	else if ((tet.GE.B(2)).AND.(tet.LE.B(1))) then
	l=1
	end if
	
	SELECT CASE (l)
	CASE (1)
	EX =  (A - (D + XL) * TAN(tet)) / SIN(tet)
	CASE (2)
	EX = - XL / COS(tet)
	CASE (3)
	EX =  (D * TAN(tet) - R) / SIN(tet)
	END SELECT
	
	if ((phi.GE.B(4)).AND.(phi.LT.B(3))) then
	k=3
	else if ((phi.GE.B(3)).AND.(phi.LT.B(2))) then
	k=2
	else if ((phi.GE.B(2)).AND.(phi.LE.B(1))) then
	k=1
	end if
	
	SELECT CASE (k)
	CASE (1)
	EY = (A - (D + XL) * TAN(phi)) / SIN(phi)
	CASE (2)
	EY = -XL / COS(phi)
	CASE (3)
	EY = (D * TAN(phi) - R) / SIN(phi)
	END SELECT

	x=SQRT((EX**2) + (EY**2))
	t1=1.0d0-EXP(-TAUU(en)*x)		
			
!	IF (l.EQ.1) THEN
!		t1 = t1 * EXP(-TAUU(30.0D0) * (A / SIN(tet) - D / COS(tet)))
!	END IF
	
	
	q0=q0+t1
	q2=q2+t1*legP(2,0,COS((PI/180.0D0)*teta+m*tet))
	q4=q4+t1*legP(4,0,COS((PI/180.0D0)*teta+m*tet))
	end do
!	q2=q2/q0
!	q4=q4/q0
	q2=q2/(q0*legP(2,0,COS((PI/180.0D0)*teta)))
	q4=q4/(q0*legP(4,0,COS((PI/180.0D0)*teta)))
	
!	t2=t2+q2
!	t3=t3+q4
	Write(14,*) teta, en, q2
	Write(15,*) teta, en, q4
	teta = teta+0.10d0
	end do
	write(*,*) en
	en = en + 0.1d0
	end do	
!	t2=t2/18001.0D0
!	t3=t3/18001.0D0
!	Write(*,*)'avg q2,q4'
!	Write(*,*) t2,t3	
	CLOSE(14)
	CLOSE(15)	
	end do
	
	end program



