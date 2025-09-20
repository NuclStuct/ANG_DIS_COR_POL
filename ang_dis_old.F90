	program ang_dist
	USE COM
	USE MATH
	USE QM_AM
	USE ANG_DIS
	implicit none
!	real*8,external::cgc,racahW,Fk,rho,Pm,Bk,alpha,Akmax,legP
!	real*8,external::Fk,rho,Pm,Bk,alpha,Akmax
	character*15::f2
	real*8::Jf,L1,L2,Ji,d,a2,a4,chisq,wcalc,wexp,da2,da4,t1,t2,t3,r
	integer::i,j,k
	real(8)::tet,dtet,ran,t4,t5,t6,teta(5)
	call fact(f)
	
	teta = (/ 90.0d0, 115.0d0, 125.0d0, 135.0d0, 145.0d0 /)
	dtet = 	6.99676637d0
	
	do
	
	t4 = 0
	t5 = 0
	t6 = 0
	
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
	
	
	d=-10.0d0

	do i=1,1000

!	teta=0.0d0
	chisq=0.0d0
	wcalc = 0.0d0
	wexp = 0.0d0
	t2 = 0.0d0
	
	do j=1,5
	
	do k=1,1000
	
	call random_number(r)
	
	ran = -1.0d0*dtet + 2.0d0*r*dtet
	tet = teta(j)+ran

	t4 = W(Ji,L1,L2,Jf,d,tet)
	t5 = W_exp(a2,a4,tet)
	t6 = W_exp(da2,da4,tet)-1.0d0
	
	wcalc = wcalc + t4
	wexp = wexp + t5
	t2 = t2 + t6
	
	end do	

!	wcalc=1.0d0+alpha(2,Ji)*Akmax(2,Ji,L1,L2,Jf,d)*legP(2,0,DCOSD(teta))+alpha(4,Ji)*Akmax(4,Ji,L1,L2,Jf,d)*legP(4,0,DCOSD(teta))
!	wexp=1+a2*legP(2,0,DCOSD(teta))+a4*legP(4,0,DCOSD(teta))
	
!	t1=(wcalc-wexp)**2
!	t2=legP(2,0,DCOSD(teta))*da2+legP(4,0,DCOSD(teta))*da4
!	t2 = W_exp(da2,da4,tet)-1.0d0 
!	t3=DSIND(teta(j))*(a2*(3.0d0*DCOSD(teta(j)))+a4*0.5d0*(35.0d0*(DCOSD(teta(j))**3)-15.0d0*DCOSD(teta(j))))*6.99676637d0
!	write(*,*) t3

	!error in teta is 6.99676637d0 degrees

	t1=(wcalc-wexp)**2
	t3=0.d0

	if (t2.NE.0.0d0) then
	chisq= chisq+t1/((t2+t3)**2.0d0)
	end if
!	teta=teta+1.0d0
	end do
	
	write(14,*) d, chisq
	d=d+0.01d0
	
	end do
	
	CLOSE(14)
	
	end do
	
	end program

!function for cgc
!	real*8 function cgc(j1,j2,j,m1,m2)
!	USE COM
!	implicit none
!	real*8, intent(in)::j1,j2,j,m1,m2
!	real*8::t1,t2,t3,m!,f(0:50)
!	integer::a,b,c,d,e,imin,i
!	m=m1+m2
!	a=(j1+j2-j)
!	b=(j1-m1)
!	c=(j2+m2)
!	d=(j-j2+m1)
!	e=(j-j1-m2)
!	call factorial(f)
!	t1=sqrt((2*j+1)*f(int(j+j1-j2))*f(int(j-j1+j2))*f(int(j1+j2-j))/f(int(j+j1+j2+1)))
!	t2=sqrt(f(int(j1+m1))*f(int(j1-m1))*f(int(j2+m2))*f(int(j2-m2))*f(int(j-m))*f(int(j+m)))
!	imin=min(d,e)
!	if(imin.lt.0)then
!	  imin=iabs(imin)
!	  else
!	    imin=0
!	end if
!	t3=0.0
!	do i=imin,min(a,b,c)
!	  t3=t3+((-1)**i)/real(f(i)*f(a-i)*f(b-i)*f(c-i)*f(d+i)*f(e+i))
!	end do
!	cgc=t1*t2*t3
!	end function


!delta  used in racah W
!	real*8 function delta(a,b,c)
!	implicit none
!	real*8, intent(in)::a,b,c
!	real*8::k,f(0:50)
!	call factorial(f)
!	k=(f(int(a+b-c))*f(int(a-b+c))*f(int(b+c-a)))/f(int(a+b+c+1))
!	delta=sqrt(k)
!	end function


!function for racah W 
!pg 81 devanathan
!	real*8 function racahW(a,b,c,d,e,f)
!	implicit none
!	real*8, intent(in)::a,b,c,d,e,f
!	real*8,external::delta
!	integer::aa,bb,cc,dd,ee,ff,gg,imin,i
!	real*8::t1,t2,k(0:50)
!	aa=(-a-b-e)
!	bb=(-c-d-e)
!	cc=(-a-c-f)
!	dd=(-b-d-f)
!	ee=(a+b+c+d)
!	ff=(a+d+e+f)
!	gg=(b+c+e+f)
!	call factorial(k)
!	t1=delta(a,b,e)*delta(c,d,e)*delta(a,c,f)*delta(b,d,f)
!	imin=min(aa,bb,cc,dd)
!	if(imin.lt.0)then
!	  imin=iabs(imin)
!	  else
!	    imin=0
!	end if
!	do i=imin,min(ee,ff,gg)
!	  t2=t2+((-1)**(ee+i))*k(i+1)/real(k(aa+i)*k(bb+i)*k(cc+i)*k(dd+i)*k(ee-i)*k(ff-i)*k(gg-i),8)
!	end do
!	racahW=t1*t2
!	end function
!
!
!fuction for Fk
!ref-1
!	real*8 function Fk(k,Jf,L1,L2,Ji)
!	USE QM_AM
!	implicit none
!	integer, intent(in)::k
!	real*8, intent(in)::Jf,L1,L2,Ji
!	real*8,external::cgc,racahW
!	real*8::t1,t2
!	t1=(-1)**(Jf-Ji-1)
!	t2=sqrt((2*L1+1)*(2*L2+1)*(2*Ji+1))
!	Fk=t1*t2*cgc(L1,L2,real(k,8),1.0d0,-1.0d0)*racah_W(Ji,Ji,L1,L2,real(k,8),Jf)
!	end function
	

!function for Population Parameter
!sigma might cause error 
!took sigm/j =.3
!not sure if external variable can just be a variable
!not a function 
!	real*8 function Pm(m,J)
!	USE COM
!	implicit none
!	real*8, intent(in)::m,J
!	!real*8,external::sigma
!	real*8::t1,t2,mm
!	t1=EXP(-1*((m/(sbyj*J))**2)/2)
!	t2=0.0
!	mm=-1.0d0*abs(J)
!	do
!	t2=t2+EXP(-1*((mm/(sbyj*J))**2)/2)
!	mm=mm+1
!	if(mm.GT.abs(J)) exit
!	end do
!	Pm=t1/t2
!	end function


!function for statistical tensor
!	real*8 function rho(k,J)
!	USE QM_AM
!	implicit none
!	integer,intent(in)::k
!	real*8,intent(in)::J
!	real*8::t1,t2,m
!	real*8,external::Pm
!	m=-1.0d0*abs(J)
!	t1=sqrt(2.0d0*J+1)
!	t2=0.0d0
!	do
!	t2=t2+(-1)**(J-m)*cgc(J,J,real(k,8),m,-m)*Pm(m,J)
!	m=m+1
!	if(m.GT.abs(J)) exit
!	end do
!	rho=t1*t2
!	end function

!function for Bk statistical tensor for complete alignment
!	real*8 function Bk(k,J)
!	USE QM_AM
!	implicit none
!	integer,intent(in)::k
!	real*8,intent(in)::J
!	real*8,external::cgc
!	real*8::t1
!	t1=sqrt(2.0d0*J+1)
!	if(mod(J,1.0d0).EQ.0.0d0) then
!	Bk=t1*cgc(J,J,real(k,8),0.0d0,0.0d0)*(-1)**(J)
!	else
!	Bk=t1*cgc(J,J,real(k,8),0.5d0,-0.5d0)*(-1)**(J-0.5d0)
!	end if
!	end function


!fuction for Ak with max alignment
!	real*8 function Akmax(k,Jf,L1,L2,Ji,d)
!	implicit none
!	integer,intent(in)::k
!	real*8,intent(in)::Jf,L1,L2,Ji,d
!	real*8,external::Fk, Bk
!	Akmax=(1/(1+d*d))*Bk(k,Ji)*(Fk(k,Jf,L1,L1,Ji)+2*d*Fk(k,Jf,L1,L2,Ji)+d*d*Fk(k,Jf,L2,L2,Ji))
!	end function
	
!	real*8 function alpha(k,J)
!	implicit none
!	integer,intent(in)::k
!	real*8,intent(in)::J
!	real*8,external::rho, Bk
!	alpha=rho(k,J)/Bk(k,J)
!	end function

!function for legendre polynomial with k=2 and 4
!	real*8 function legP(k,teta)
!	implicit none
!	integer,intent(in)::k
!	real*8,intent(in)::teta
!	real*8::cos
!	cos=dcosd(teta)
!	select case (k)
!		case (2)
!		legP=0.5d0*(3.0d0*(cos**2)-1.0d0)
!		case (4)
!		legP=(1/8.0d0)*(35.0d0*(cos**4)-30.0d0*(cos**2)+3.0d0)
!	end select
!	end function
!
!factorial subroutine
!	subroutine factorial(f)
!	implicit none
!	integer::i
!	real*8::f(0:50)
!	f(0)=1.0d0
!	do i=1,50
!	  f(i)=f(i-1)*real(i)
!	end do
!	end subroutine

