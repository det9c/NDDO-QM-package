subroutine finduhf
use indices
use constants
use tables
     use control
use scratch_array
implicit double precision (a-h,o-z)
interface

subroutine square(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::matrix
double precision,dimension(:),intent(in)::vector
integer,intent(in)::idim1
end subroutine square

 subroutine twoe_pair(i,j,t1)
 integer,intent(in)::i,j
 double precision,intent(inout)::t1
end subroutine twoe_pair

subroutine hcoregrad(loop)
integer,dimension(2),intent(in)::loop
end subroutine hcoregrad

subroutine f1grad(loop,h,fock,density,twoe)
double precision,dimension(:),intent(inout)::fock
integer,dimension(2),intent(in)::loop
double precision,dimension(:),intent(in)::density,twoe,h
end subroutine f1grad

subroutine f2graduhf(loop,fock,density,twoe)
integer,dimension(2),intent(in)::loop
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density,twoe
end subroutine f2graduhf

subroutine trace(S,N,out)
double precision,dimension(:,:),intent(in)::S
integer::N
double precision,intent(out)::out
end subroutine trace

subroutine hcoregradsp(i)
integer,intent(in)::i
end subroutine hcoregradsp

      subroutine pack1(indi,indj,ij)
integer,intent(in)::indi,indj
integer,intent(inout)::ij
end subroutine pack1


end interface


double precision,allocatable,dimension(:,:)::grad,grad2
!double precision,intent(inout),dimension(:)::g
double precision::delta,xold,yold,zold
double precision,dimension(:,:),allocatable::pone,ptwo,hlocal  & ! i.e. p-local,h-local etc
,temp,temp2,temp3,pthr
double precision,dimension(:),allocatable::focka,fockb
double precision,dimension(2)::force
integer,dimension(2)::loop


!  matrix s3 holds the frozen density from the SCF calculation
allocate(grad(3,numat))
if(.not.save_tree)then
print*,'Computing gradient using UHF determinant...'
end if
call cpusec(time1)
delta=.001
!delta=1d-10
eold=etotal

allocate(focka(ndim1))
allocate(fockb(ndim1))
allocate(scrvec(ndim1))
allocate(scrmat(num_basis,num_basis))
grad=zero
fock=zero
do i=1,numat-1

i1=ifirst(i)
i2=ilast(i)
aa=zero
bb=zero
cc=zero
dd=zero
 

 do j=i+1,numat
j1=ifirst(j)
j2=ilast(j)
ni=nbas(species(i))
nj=nbas(species(j))
loop(1)=i
loop(2)=j

xold=x(i)
x(i)=xold+delta

!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,focka,s3,twoe)
call  f2graduhf(loop,focka,adens,twoe)
call f1grad(loop,h,fockb,s3,twoe)
call  f2graduhf(loop,fockb,bdens,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+two*(s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij))
   end do
   end do



force(1)=half*energy+t1



x(i)=xold-delta
!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,focka,s3,twoe)
call  f2graduhf(loop,focka,adens,twoe)
call f1grad(loop,h,fockb,s3,twoe)
call  f2graduhf(loop,fockb,bdens,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+two*(s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij))
   end do
   end do


force(2)=half*energy+t1


der=(force(1)-force(2))/(two*delta)
grad(1,i)=grad(1,i)+der
grad(1,j)=grad(1,j)-der

x(i)=xold


yold=y(i)
y(i)=yold+delta

!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,focka,s3,twoe)
call  f2graduhf(loop,focka,adens,twoe)
call f1grad(loop,h,fockb,s3,twoe)
call  f2graduhf(loop,fockb,bdens,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+two*(s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij))
   end do
   end do




force(1)=half*energy+t1



y(i)=yold-delta
!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,focka,s3,twoe)
call  f2graduhf(loop,focka,adens,twoe)
call f1grad(loop,h,fockb,s3,twoe)
call  f2graduhf(loop,fockb,bdens,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+two*(s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij))
   end do
   end do

force(2)=half*energy+t1


der=(force(1)-force(2))/(two*delta)
grad(2,i)=grad(2,i)+der
grad(2,j)=grad(2,j)-der

y(i)=yold


zold=z(i)
z(i)=zold+delta

!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,focka,s3,twoe)
call  f2graduhf(loop,focka,adens,twoe)
call f1grad(loop,h,fockb,s3,twoe)
call  f2graduhf(loop,fockb,bdens,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+two*(s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij))
   end do
   end do




force(1)=half*energy+t1



z(i)=zold-delta
!vdbh=zero
call twoe_pair(i,j,t1)
call hcoregrad(loop)
!vdbfock=zero
call f1grad(loop,h,focka,s3,twoe)
call  f2graduhf(loop,focka,adens,twoe)
call f1grad(loop,h,fockb,s3,twoe)
call  f2graduhf(loop,fockb,bdens,twoe)


energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=j1,j2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij)
   end do
   end do

   do itemp=i1,i2
   do jtemp=j1,j2
      call pack1(jtemp,itemp,ij)
       energy=energy+two*(s3(ij)*two*h(ij)+adens(ij)*focka(ij)+bdens(ij)*fockb(ij))
   end do
   end do


force(2)=half*energy+t1


der=(force(1)-force(2))/(two*delta)
grad(3,i)=grad(3,i)+der
grad(3,j)=grad(3,j)-der

z(i)=zold







end do




end do



if(sparkles)then
do i=1,numat
i1=ifirst(i)
i2=ilast(i)

xold=x(i)
x(i)=xold+delta
enuc=zero
h=zero
call hcoregradsp(i)

energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do




force(1)=half*energy+enuc








x(i)=xold-delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do





force(2)=half*energy+enuc

der=(force(1)-force(2))/(two*delta)
grad(1,i)=grad(1,i)+der


x(i)=xold

yold=y(i)
y(i)=yold+delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do





force(1)=half*energy+enuc




y(i)=yold-delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do




force(2)=half*energy+enuc



der=(force(1)-force(2))/(two*delta)
grad(2,i)=grad(2,i)+der


y(i)=yold


zold=z(i)
z(i)=zold+delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do




force(1)=half*energy+enuc
z(i)=zold-delta
enuc=zero
h=zero
call hcoregradsp(i)
energy=zero
   do itemp=i1,i2
   do jtemp=i1,i2
      call pack1(jtemp,itemp,ij)
      energy=energy+s3(ij)*two*h(ij)
   end do
   end do



force(2)=half*energy+enuc



der=(force(1)-force(2))/(two*delta)
grad(3,i)=grad(3,i)+der


z(i)=zold

end do


end if










grad=grad*627.51d0/27.21d0


if(allocated(g))then
do i=1,numat
      ISTART=3*(I-1)+1
      g(istart)=grad(1,i)
      g(istart+1)=grad(2,i)
      g(istart+2)=grad(3,i)
end do
end if















rnorm=zero
if(.not.save_tree)then
print*,'                     Gradient'
print*,' Atom      dE/dx            dE/dy         dE/dz'
print*,'--------------------------------------------------------------'
do i=1,numat
write(*,10)i,grad(1,i),grad(2,i),grad(3,i)
rnorm=rnorm+grad(1,i)**2+grad(2,i)**2+grad(3,i)**2
end do 

rnorm=dsqrt(rnorm)
print*,'||F|| =',rnorm
print*,'--------------------------------------------------------------'
10 format(i4,f15.4,f15.4,f15.4)
call cpusec(time2)
write(*,202)time2-time1
202 format(' Elapsed time for gradient calculation (sec) =  ',F10.3)
end if
rnorm=0.0d0
do i=1,numat
rnorm=rnorm+grad(1,i)**2+grad(2,i)**2+grad(3,i)**2
end do

rnorm=dsqrt(rnorm)

open(unit=1,file='cartnorm')
write(1,*)rnorm
close(1)


if(.not.save_tree)print*,''
if((OPTIMIZE .or. GRADIENT .or. FREQUENCY) .and. (.not. XYZ))then
numint=3*ninput-6
numcart=3*ninput
rnorm2=zero
if(allocated(gint))deallocate(gint)
allocate(gint(numint))
gint=zero

icount=0
allocate(grad2(3,ninput))
do i=1,ninput
if(zstore(i).eq.0)then
grad2(1,i)=zero
grad2(2,i)=zero
grad2(3,i)=zero
else
icount=icount+1
grad2(1,i)=grad(1,icount)
grad2(2,i)=grad(2,icount)
grad2(3,i)=grad(3,icount)
end if
end do
call dgemv( 'N', numint,numcart,one,bmatrix,numint,grad2,1,zero,gint,1)
deallocate(grad2)
IF(.not. FREQUENCY .and. .not. TS .and. .not. MODES)then
write(*,*)'Internal coordinate gradients Kcal/(mol*Ang)'
write(*,*)'Type     Value       Gradient'
write(*,*)'-----------------------------'
open(unit=98,file='intgrad')
do i=1,numint
write(98,*)gint(i)
!write(*,*)inttype(i),q(intadd(i,1),intadd(i,2)),gint(i)
fac=1.0d0
if(inttype(i).ne.'Bond')fac=180.0d0/pi
write(*,40)inttype(i),q(intadd(i,1),intadd(i,2))*fac,gint(i)
rnorm2=rnorm2+gint(i)*gint(i)
end do
close(98)
rnorm2=dsqrt(rnorm2)
print*,'Internal Gradient Norm =',rnorm2
40 format(A,f12.5,f15.10)

open(unit=1,file='intnorm')
write(1,*)rnorm2
close(1)

end if
end if













deallocate(grad)
deallocate(scrvec)
deallocate(scrmat)



end subroutine finduhf



