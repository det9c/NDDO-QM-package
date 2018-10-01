subroutine make_hess
use constants
use tables
use scratch_array
use control
implicit double precision (a-h,o-z)
interface
subroutine scfopt(x0,y0,z0,gout,natoms,ff)
double precision,intent(inout) ::ff
double precision, dimension(natoms),intent(inout) ::x0,y0,z0
double precision, dimension(3*natoms),intent(inout) ::gout
end subroutine scfopt


subroutine projector(hessian,i3n)
integer,intent(in)::i3n
double precision,dimension(i3n,i3n),intent(inout)::hessian
end subroutine projector

subroutine inertia(reori)
logical,intent(in)::reori
end subroutine inertia

end interface

double precision,allocatable,dimension(:,:)::hessian,gplus,gminus,hessold
double precision,allocatable,dimension(:)::vec2,amass,oldfreq,scr2
!double precision,allocatable,dimension(:)::x0,y0,z0
factor=5140.36636949d0 ! conversion factor from mass weighted hessian to cm(-1)
! to derive: (H*mol/bohr**2 grams) are units of eigenvalues of mass weighted
! hessian.  Conversion is then:
! 627.51 kcal/mol * 1000 cal/kcal * 4.184 J/cal * 1bohr**2/.529177(-10)**2 m 
! * 1000 g/kg
!take square root of all the above and divide by 2*pi*c where c is speed of light in cm/s
delta=1D-3
i3n=3*numat

print*,'-----------------------------------------------'
print*,'|                                             |'
print*,'|           VIBRATIONAL FREQUENCY             |'
print*,'-----------------------------------------------'
print*,''
allocate(hessian(i3n,i3n))
allocate(hessold(i3n,i3n))
allocate(gplus(i3n,1))
allocate(gminus(i3n,1))
allocate(oldfreq(i3n))
optimize=.false. ! set this to fool program into not printing xyz all the 
!time.  see subroutine scfopt.f90 
!allocate(x0(numat))
!allocate(y0(numat))
!allocate(z0(numat))
!x0=x
!y0=y
!x0=z
call cpusec(time1)
write(*,*)'NORMAL MODES WILL BE EVALUATED AT CURRENT GEOMETRY'
write(*,*)'IRRESPECTIVE OF MAGNITUDE OF FIRST DERIVATIVES...'
hessian=zero
icol=0
allocate(amass(i3n))
! evaluate initial gradient

SAVE_TREE=.true.
print*,''
call inertia(.true.)
if(.not.xyz)call buildb
write(*,*)'COMPUTING SECOND DERIVATIVE MATRIX...'
do i=1,numat
!print*,'Working on components for atom',i
icol=icol+1
! compute d2E/[dx(i) dq(j)] where q=x,y,z
xold=x(i)
x(i)=x(i)+delta
call scfopt(x,y,z,gplus,numat,ff)
x(i)=x(i)-2.0d0*delta
call scfopt(x,y,z,gminus,numat,ff)
gplus=(gplus-gminus)/two/delta
x(i)=xold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)
amass(icol)=atmass(zeff(i))

icol=icol+1
! compute d2E/[dy(i) dq(j)] where q=x,y,z
yold=y(i)
y(i)=y(i)+delta
call scfopt(x,y,z,gplus,numat,ff)
y(i)=y(i)-2.0d0*delta
call scfopt(x,y,z,gminus,numat,ff)
gplus=(gplus-gminus)/two/delta
y(i)=yold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)
amass(icol)=atmass(zeff(i))


icol=icol+1
! compute d2E/[dz(i) dq(j)] where q=x,y,z
zold=z(i)
z(i)=z(i)+delta
call scfopt(x,y,z,gplus,numat,ff)
z(i)=z(i)-2.0d0*delta
call scfopt(x,y,z,gminus,numat,ff)
gplus=(gplus-gminus)/two/delta
z(i)=zold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)
amass(icol)=atmass(zeff(i))
end do


hessian=hessian*.529177d0*.529177d0/627.51d0 ! convert to atomic units

! symmetrize hessian
hessian=hessian+transpose(hessian)
hessian=hessian/two


! mass weight hessian
amass=1.0d0/dsqrt(amass)
do i=1,i3n
do j=1,i3n
hessian(j,i)=hessian(j,i)*amass(j)*amass(i)
end do
end do


!diagonalize hessian using tred3/tql3
if(allocated(scrmat))deallocate(scrmat)
if(allocated(scrvec))deallocate(scrvec)
allocate(scrmat(i3n,i3n))
allocate(scrvec(i3n))
allocate(vec2(i3n))

hessold=hessian
call tred3(i3n,i3n,hessian,scrvec,vec2,scrmat)
call tql3(i3n,i3n,scrvec,vec2,scrmat,iout) ! eigenvals are in scrvec
!call matprt(scrmat,i3n,i3n,i3n,i3n)



do i=1,i3n
scrvec(i)=dsqrt(abs(scrvec(i)))*factor*dsign(one,scrvec(i))
end do
oldfreq=scrvec

!scrvec=abs(scrvec)
!scrvec=dsqrt(scrvec)*factor



hessian=hessold
do i=1,i3n
do j=1,i3n
if(abs(hessian(j,i)).lt.1D-5)hessian(j,i)=0.0d0
end do
end do

call projector(hessian,i3n)

open(unit=1,file='hessian.out')
write(1,*)hessian
close(1)


if(allocated(scrmat))deallocate(scrmat)
if(allocated(scrvec))deallocate(scrvec)
allocate(scrmat(i3n,i3n))
allocate(scrvec(i3n))
call tred3(i3n,i3n,hessian,scrvec,vec2,scrmat)
call tql3(i3n,i3n,scrvec,vec2,scrmat,iout) ! eigenvals are in scrvec
!call eig(hessian,scrmat,i3n,i3n,0)
!do i=1,i3n
!scrvec(i)=hessian(i,i)
!print*,scrvec(i)
!end do
!print*,'eigenvectors of mass weighted hessian'
!call matprt(scrmat,i3n,i3n,i3n,i3n)


if(IRC)then
allocate(transvec2(i3n))
call dcopy(i3n,scrmat(1,1),1,transvec2,1)
transvec2=transvec2*autoang
rnorm=zero
rnorm=ddot(i3n,transvec2,1,transvec2,1)
rnorm=one/dsqrt(rnorm)
call dscal(i3n,rnorm,transvec2,1)
rnorm=ddot(i3n,transvec2,1,transvec2,1)
print*,'norm n her',dsqrt(rnorm)
end if






! convert vectors to cartesian displacement vectors in BOHR
do i=1,i3n
do j=1,i3n
scrmat(j,i)=scrmat(j,i)*amass(j)
end do
end do






!normalize vecs
do i=1,i3n
rnorm=zero
rnorm=ddot(i3n,scrmat(1,i),1,scrmat(1,i),1)
rnorm=one/dsqrt(rnorm)
call dscal(i3n,rnorm,scrmat(1,i),1)
end do





!print*,'Cartesian displacement vectors (atomic units)!!!'
!call matprt(scrmat,i3n,i3n,i3n,i3n)

!do i=1,i3n
!print*,scrmat(i,1)
!end do


! transform eigenvecs to internal coordinates
!this won't work with dummy atoms
numint=3*numat-6
if(allocated(hessold))deallocate(hessold)
allocate(hessold(numint,i3n))
if(.not.xyz)call dgemm( 'N', 'N', numint, i3n, i3n, one,bmat4int &
 ,numint,scrmat, i3n,zero,hessold,numint )

!normalize vecs
do i=1,i3n
rnorm=zero
rnorm=ddot(numint,hessold(1,i),1,hessold(1,i),1)
rnorm=one/dsqrt(rnorm)
call dscal(numint,rnorm,hessold(1,i),1)
end do


! print*,'eigevecs of hessian in internals'
!call matprt(hessold,numint,i3n,numint,i3n)


if(IRC)then
allocate(transvec(numint))
call dcopy(numint,hessold(1,1),1,transvec,1)

end if






do i=1,i3n
scrvec(i)=dsqrt(abs(scrvec(i)))*factor*dsign(one,scrvec(i))
end do

print*,'*********************************************************'
print*,'Table of vibrational frequencies (cm-1) before'
print*,'and after projection of translational/rotational modes'
print*,'*********************************************************'

write(*,*)'MODE      ','Before projection','  After projection'
write(*,*)'----      ','----------------','  ----------------'
do i=1,i3n
write(*,90)i,oldfreq(i),scrvec(i)
end do


open(unit=1,file='mode')
do i=1,i3n
write(1,*)scrvec(i)
end do
close(1)







90 format(i4,f20.4,f20.4)

call cpusec(time2)
write(*,*)'NORMAL MODE ANALYSIS REQUIRED',time2-time1,' seconds'

deallocate(scrmat)





end subroutine make_hess


