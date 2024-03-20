program ab
implicit none
integer*8 i,imax
real*8 tmax,t,dt,x1,x2,dx1,dx2,k,xs1,xs2,dxs1,dxs2
real*8 v1,v2,dv1,dv2,tau,z1,z2,u1,u2,delx1,delx2
parameter(tmax=1500.,dt=1.e-3, tau=100.)!gam = 0.1
parameter(imax=int(tmax/dt),k=10.)
real*8, dimension(imax) :: noise,adap1,adap2

!make noise
call random_seed()  

do i = 1, imax / 2
call random_number(u1)
call random_number(u2)

z1 = sqrt(-2. * log(u1)) * cos(2.*3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2.*3.14 * u2)

noise(2 * i - 1) = z1
noise(2 * i) = z2
end do

!make file
open(1,file='aa')
open(2,file='bb')

!!!!!!!!초기조건!!!!!!!!!!!!
t=0.
x1=1.;v1=1.;xs1=1. 
x2=1.;v2=1.;xs2=1. 

do i=1,imax
write(1,*) t,x2
write(2,*) t,xs2

!!!!!!!!!!!!!!!!!!!!!!adaptation force 생성!!!!!!!!!!!!!!!
dxs1=(x1-xs1)/tau
dxs2=(x2-xs2)/tau
xs1=xs1+dxs1*dt
xs2=xs2+dxs2*dt
delx1=x1-xs1
delx2=x2-xs2
adap1(i)=sign(1.0,real(delx1))
adap2(i)=sign(1.0,real(delx2))

! 자발 진동 상황
dv1=-0.1*v1-x1+adap1(i)+noise(i)-k*(x2-x1)
dv2=-0.1*v2-x2+adap2(i)+noise(i)-k*(x1-x2)
dx1=v1
dx2=v2
v1=v1+dv1*dt
v2=v2+dv2*dt
x1=x1+dx1*dt
x2=x2+dx2*dt
t=i*dt
end do
end program

