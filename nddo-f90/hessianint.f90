subroutine make_hess2
use constants
use tables
use scratch_array
use control
implicit double precision (a-h,o-z)

interface
subroutine scfq(q0,natoms,ff,stat1,stat2,gout,numint)
integer,intent(in)::natoms,numint
double precision,intent(inout) ::ff
double precision, dimension(numint),intent(inout) ::gout
double precision,dimension(natoms,3),intent(in)::q0
logical::stat1,stat2
end subroutine scfq


subroutine outerp(n,alpha,v1,v2,mat)
double precision,dimension(n,1),intent(in)::v1,v2
double precision,dimension(n,n),intent(inout)::mat
double precision,intent(in)::alpha
integer,intent(in)::n
end subroutine outerp




end interface
!parameter(stpsz=0.02d0)
double precision,allocatable,dimension(:,:)::hessian,gplus,gminus,hessold &
,scr3,scr4,gold,deltah,greal
double precision,allocatable,dimension(:)::vec2,amass,fvec,scale,step,vec3,oldvec 
!double precision,allocatable,dimension(:)::x0,y0,z0
logical::stat1,stat2
factor=5140.36636949d0 ! conversion factor from mass weighted hessian to cm(-1)
! to derive: (H*mol/bohr**2 grams) are units of eigenvalues of mass weighted
! hessian.  Conversion is then:
! 627.51 kcal/mol * 1000 cal/kcal * 4.184 J/cal * 1bohr**2/.529177(-10)**2 m 
! * 1000 g/kg
!take square root of all the above and divide by 2*pi*c where c is speed of light in cm/s

print*,'-----------------------------------------------'
print*,'|          Transition State Search             |'
print*,'|       Eigenvector Following Algorithm        |'
print*,'|                                              |'
print*,'|                                              |'
print*,'-----------------------------------------------'


if(.not. CALCALL)print*,'Exact hessian will be evaluated every ',UPDATE, ' cycles.'

numint=3*ninput-6
delta=1.0D-3
i3n=numint
stpsz=ste1
delta=del1
print*,'stpsz and delta',stpsz,delta

allocate(hessian(i3n,i3n))
allocate(hessold(i3n,i3n))
allocate(deltah(i3n,i3n))
allocate(gplus(i3n,1))
allocate(gold(i3n,1))
allocate(gminus(i3n,1))
allocate(greal(i3n,1))
allocate(fvec(i3n))
allocate(scale(i3n))
allocate(step(i3n))
allocate(oldvec(i3n))

optimize=.true. ! set this to fool program into not printing xyz all the 
call cpusec(time1)
hessian=zero
allocate(amass(i3n))

SAVE_TREE=.true.


ncycle=0
! compute gradient and see if below tolerance
10 ncycle=ncycle+1


if(ncycle.gt.100)then
print*,'Max number of optimization steps exceeded'
! for srp
open(unit=12,file='opt')
write(12,*)'2'
write(12,*)rnorm2
close(12)
open(unit=10,file='geometry')
rewind 10
do i=2,ninput
write(10,*)q(i,1)
end do
do i=3,ninput
write(10,*)q(i,2)*180.0d0/pi
end do
do i=4,ninput
write(10,*)q(i,3)*180.0d0/pi
end do
close(10)
! end for srp
!if(rnorm2.lt.15.0d0)then
!rsrp=rnorm2
!return
!else
stop
!end if
end if

if(.not.GRADIENT)then
PRINT*,'-----------------------------------------------'
PRINT*,'GEOMETRY OPTIMIZATION CYCLE NUMBER: ',ncycle
write(*,*)'GEOMETRY AT CURRENT ITERATE'
do i=1,ninput
write(*,301)zstore(i),q(i,1),opt(i,1),q(i,2)*180.0d0/pi,opt(i,2),q(i,3)*180.0d0/pi, &
opt(i,3),ref(i,1),ref(i,2),ref(i,3)
end do
301 format(i3,f20.10,i3,f20.10,i3,f20.10,i3,i3,i3,i3)
end if

! check for angles over 180 degrees
do i=1,ninput
degree=q(i,2)*180.0d0/pi
if(degree .gt. 180.0d0)then
print*,'An angle is greater than 180 degrees. Reconstruct Z-matrix'
stop
end if
end do



stat1=.false.
stat2=.true.

call scfq(q,ninput,ff,stat1,stat2,gplus,i3n)
if(ncycle .eq. 1)gold=gplus
greal=gplus
write(*,*)'Internal coordinate gradients (Kcal/Mol Ang):'
write(*,*)'Type     Value       Gradient'
write(*,*)'-----------------------------'
rnorm2=zero
do i=1,numint
!write(*,*)inttype(i),q(intadd(i,1),intadd(i,2)),gint(i)
!fac2=.529177/627.51d0
fac2=1.0d0
fac=1.0d0
if(inttype(i).ne.'Bond')fac=180.0d0/pi
!if(inttype(i).ne.'Bond')fac2=1.0d0/627.51d0
write(*,70)inttype(i),q(intadd(i,1),intadd(i,2))*fac,gplus(i,1)*fac2
rnorm2=rnorm2+gplus(i,1)*gplus(i,1)*fac2*fac2
end do
rnorm2=dsqrt(rnorm2)
print*,'Internal Coordinate Gradient Norm =',rnorm2
70 format(A,f12.5,f15.5)

! compute gradient norm of optimized coordinates
rnorm3=zero
do i=1,nopt
   do j=1,numint
      if(intadd(j,1).eq.internal(1,i) .and. intadd(j,2).eq.internal(2,i))exit
   end do
 rnorm3=rnorm3+gplus(j,1)**2
end do
rnorm3=dsqrt(rnorm3)
print*,'Optimized Coordinate Gradient Norm =',rnorm3



converge=OPTTOL ! convert to hartree/bohr
if(TS .and. rnorm3.lt.converge)then
print*,'-------------------------------------------------'
print*,'|Based on the optimized coordinate gradient norm |'
print*,'|being below the tolerance, a transition state   |'
print*,'|may have been found.                            |'
print*,'-------------------------------------------------'
rsrp=rnorm3
return
end if




icol=0
iexact=mod(ncycle,UPDATE)
if(ncycle .gt. 1   .and.   .not. CALCALL .and. iexact.ne.0)then ! do a powell update of hessian
term1 = ddot(i3n,step,1,step,1)
term1= one/term1

call dcopy(i3n,gplus,1,fvec,1)
call daxpy(i3n,-1.0d0,gold,1,fvec,1)
call dcopy(i3n,gplus,1,gold,1) ! save gradient for next iteration
call dgemv( 'N',i3n,i3n,-1.0d0,hessold,i3n,step,1,one,fvec,1) ! so fvec is V from paper
term2=ddot(i3n,fvec,1,step,1)
term2=-term2/(term1*term1)
call dcopy(i3n,fvec,1,gplus,1)
call dcopy(i3n,step,1,gminus,1)
deltah=zero
call outerp(i3n,1.0d0,gplus,gminus,deltah)
!deltah=deltah+transpose(deltah)
call outerp(i3n,1.0d0,gminus,gplus,deltah)
call outerp(i3n,term2,gminus,gminus,deltah)
deltah=term1*deltah
hessian=hessold+deltah


else ! compute new hessian by finite difference


   if(ncycle.eq.1 .and. FORWARD)then
      print*,'Fast evaluation of first hessian'
      do i=1,numint
      stat1=.false.
      stat2=.true.
      icol=icol+1
      xold=q(intadd(i,1),intadd(i,2))
      q(intadd(i,1),intadd(i,2))=q(intadd(i,1),intadd(i,2))+delta
      call scfq(q,ninput,ff,stat1,stat2,gminus,numint)
      gminus=(gminus-gplus)/delta
      q(intadd(i,1),intadd(i,2))=xold
      hessian(1:i3n,icol:icol)=gminus(1:i3n,1:1)
   end do


   else






   do i=1,numint 
      stat1=.false.
      stat2=.true.
      icol=icol+1
      xold=q(intadd(i,1),intadd(i,2))
      q(intadd(i,1),intadd(i,2))=q(intadd(i,1),intadd(i,2))+delta
      call scfq(q,ninput,ff,stat1,stat2,gplus,numint)
      q(intadd(i,1),intadd(i,2))=q(intadd(i,1),intadd(i,2))-2.0d0*delta
      stat1=.false.
      stat2=.true.

      call scfq(q,ninput,ff,stat1,stat2,gminus,numint)

      gplus=(gplus-gminus)/two/delta
      q(intadd(i,1),intadd(i,2))=xold
      hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)
   end do
   end if
end if





! symmetrize hessian
hessian=hessian+transpose(hessian)
hessian=hessian/two
hessold=hessian
hessian=hessian/627.51d0 ! convert to atomic units
do i=1,i3n
if(inttype(i).eq.'Bond')then
hessian(1:i3n,i:i)=hessian(1:i3n,i:i)*.529177
hessian(i:i,1:i3n)=hessian(i:i,1:i3n)*.529177
end if
end do


!print*,'fcm with fd of gint'
!call matprt(hessian,i3n,i3n,i3n,i3n)


!diagonalize hessian using tred3/tql3
if(allocated(scrmat))deallocate(scrmat)
if(allocated(scrvec))deallocate(scrvec)
allocate(scrmat(i3n,i3n))
allocate(scrvec(i3n))
allocate(vec2(i3n))
call tred3(i3n,i3n,hessian,scrvec,vec2,scrmat)
call tql3(i3n,i3n,scrvec,vec2,scrmat,iout) ! eigenvals are in scrvec
!call matprt(scrmat,i3n,i3n,i3n,i3n)

! check eigenvals and print warning messages
ineg=0
do i=1,i3n
if(scrvec(i).lt.0.0d0)ineg=ineg+1
end do
PRINT*,'Force constant matrix has ',ineg,' negative eigenvalues'
print*,''
PRINT*,'Modes corresponding to negative eigenvalues are:'
print*,'-------------------------------------------------'
do i=1,i3n
if(scrvec(i).lt.0.0d0)then
   write(*,*)'Eigenvalue:',scrvec(i)
   print*,'Atom    Coord. Type    Value   Vector Component'
   do j=1,i3n
      fac=1.0d0
      if(inttype(j).ne.'Bond')fac=180.0d0/pi
      write(*,40)intadd(j,1),'          ',inttype(j),'  ',q(intadd(j,1),intadd(j,2))*fac,scrmat(j,i)
   end do
   print*,''
end if
end do
40 format(i3,A,A,A,f10.4,f10.4,f10.4)
if(MODES)return

! compute F=U^T |gradient> 
!allocate(fvec(i3n))
!allocate(scale(i3n))
!allocate(step(i3n))
!convert gradient to hartrees
greal=greal/627.51d0
do j=1,i3n
if(inttype(j).eq.'Bond')greal(j,1)=greal(j,1)*.529177d0
end do
call dgemv( 'T', i3n,i3n,one,scrmat,i3n,greal,1,zero,fvec,1)


! determine lambda(p)
imode=1
if(ncycle.gt.1)overlap= ddot(i3n,scrmat(1,imode),1,oldvec,1)
print*,'New mode being followed has ',overlap,' overlap with previous mode'
call dcopy(i3n,scrmat(1,imode),1,oldvec,1)

root= scrvec(imode)*scrvec(imode)  + four * fvec(imode) * fvec(imode)
rlamp=half*scrvec(imode) + half * dsqrt(root)



allocate(scr3(i3n,i3n))
allocate(scr4(i3n,i3n))
allocate(vec3(i3n))
scr3=zero
do i=2,i3n
scr3(i-1,i-1)=scrvec(i)
scr3(i-1,i3n)=fvec(i)
scr3(i3n,i-1)=fvec(i)
end do

call tred3(i3n,i3n,scr3,vec3,vec2,hessian)
call tql3(i3n,i3n,vec3,vec2,scr3,iout) ! eigenvals are in scrvec

print*,'RFO step: Lambda0= ',rlamp,' LambdaN= ',vec3(1) 
open(unit=1,file='lambda')
write(1,*)rlamp
write(1,*)vec3(1)
close(1)

! scale factors
scale=zero
scale(imode)=fvec(imode)/(scrvec(imode)-rlamp)

do i=2,i3n
scale(i)=fvec(i)/(scrvec(i)-vec3(1))
end do
scale=-scale

step=zero
do i=1,i3n
CALL DAXPY(i3n,scale(i),scrmat(1,i),1,step,1)
end do

!compute length of step
rlength=ddot(i3n,step,1,step,1)
rlength=dsqrt(rlength)
open(unit=1,file='step')
write(1,*)rlength
close(1)



if(rlength.gt.stpsz)then
factor=stpsz/rlength
print*,'Step size exceeds maximum',stpsz,' Computed step scaled by ',factor
step=factor*step
!rlength=ddot(i3n,step,1,step,1)
!rlength=dsqrt(rlength)
!print*,'new length',rlength
end if

!convert bond step lengths to angstroms
do i=1,i3n
if(inttype(i).eq.'Bond')step(i)=step(i)*autoang
end do

do i=1,nopt
     do j=1,i3n
       if(intadd(j,1).eq.internal(1,i).and. intadd(j,2).eq.internal(2,i))goto 55
     end do
55   q(  internal(1,i), internal(2,i)  ) = q(internal(1,i), internal(2,i) )+step(j)
end do

deallocate(scr3)
deallocate(scr4)
deallocate(vec3)
deallocate(vec2)
goto 10









!write(*,90)i,scrvec(i)
!90 format(i4,f20.4)


! eigenvals are in scrvec(i) and vecs are in scrmat








call cpusec(time2)


end subroutine make_hess2
