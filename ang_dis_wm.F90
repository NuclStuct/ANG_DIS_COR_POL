	program ang_dist
	USE COM
	USE MATH
	USE QM_AM
	USE ANG_DIS
	implicit none

	character*15::f2
	real*8::Jf,L1,L2,Ji,dd,a2,a4,chisq,wcalc,wexp,da2,da4,t1,t2
	integer::i,j,k
	real(8)::dtet,teta(5)

	call fact(f)
	call dfact(df)
	call expCOFF(AAA)

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
	
	WRITE(*,*) nine_J(1.0D0,1.0D0,2.0D0,1.0D0,1.0D0,2.0D0,2.0D0,2.0D0,0.0D0)
	WRITE(*,*) Fk_C(2,2,0,1.0D0,1.0D0,1.0D0,2.0D0)
	WRITE(*,*) Fk_C(2,2,0,1.0D0,2.0D0,1.0D0,2.0D0)
	WRITE(*,*) Fk_C(2,2,0,1.0D0,3.0D0,1.0D0,2.0D0)
	WRITE(*,*) Fk_C(2,2,0,2.0D0,2.0D0,1.0D0,2.0D0)
	
	
	
!	WRITE(*,*) Qk(TAUU(60.0d0),2),Qk(TAUU(60.0d0),4),TAUU(60.0d0) 
!	WRITE(*,*) Qk(70.6D0,2),Qk(70.6D0,4)
!i could add it to com	
	teta = (/ 90.0d0, 105.0d0, 120.0d0, 135.0d0, 150.0d0 /)
!	dtet = 	6.99676637d0
	dtet = 0.0d0	
	do
	
	write(*,*) "Input file name (enter 'exit' if you want to exit)"
	read(*,*) f2
	if (f2.EQ."exit") exit
	write(*,*) "Input Ji and Jf"
	read(*,*) Ji, Jf
	L1=1.0d0
!	L1=Jf-Ji
	L2=L1+1.0d0
	write(*,*) "Input experimental A2 and A4"
	read(*,*) a2, a4
	write(*,*) "Input experimental error in A2 and A4"
	read(*,*) da2, da4 
!	write(*,*) "Input teta"
!	read(*,*) teta 
	
	OPEN(UNIT=14,FILE=f2,STATUS='new')
	
	
	dd=-10.0d0

	do i=1,20000
	
	chisq=0.0d0
	
	do j=1,5
	
	

	
	wcalc = W(Ji,L1,L2,Jf,dd,teta(j))
	wexp = W_exp(a2,a4,teta(j))
	t2 = dW_exp(a2,a4,da2,da4,teta(j),dtet)
	

	t1=(wcalc-wexp)**2

	if (t2.NE.0.0d0) then
	
	chisq= chisq+t1/(t2**2.0d0)
	
	end if

	end do
	
	write(14,*) dd, chisq
	dd=dd+0.001d0
	
	end do
	
	CLOSE(14)
	
	end do
	
	end program



