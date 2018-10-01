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
end interface
parameter(stpsz=0.3d0)
double precision,allocatable,dimension(:,:)::hessian,gplus,gminus,hessold3 &
,scr3,scr4
double precision,allocatable,dimension(:)::vec2,amass,fvec,scale,step,vec3
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




numint=3*ninput-6
delta=1.0D-5
i3n=numint

allocate(hessian(i3n,i3n))
allocate(gplus(i3n,1))
allocate(gminus(i3n,1))
allocate(fvec(i3n))
allocate(scale(i3n))
allocate(step(i3n))

optimize=.true. ! set this to fool program into not printing xyz all the 
call cpusec(time1)
hessian=zero
allocate(amass(i3n))

SAVE_TREE=.true.
stat1=.false.
stat2=.true.



! compute gradient and see if below tolerance

10 call scfq(q,ninput,ff,stat1,stat2,gplus,i3n)
write(*,*)'Internal coordinate gradients Hartree/Bohr for gaussian'
write(*,*)'Type     Value       Gradient'
write(*,*)'-----------------------------'
rnorm2=zero
do i=1,numint
!write(*,*)inttype(i),q(intadd(i,1),intadd(i,2)),gint(i)
fac2=.529177/627.51d0
if(inttype(i).ne.'Bond')fac=180.0d0/pi
if(inttype(i).ne.'Bond')fac2=1.0d0/627.51d0
write(*,70)inttype(i),q(intadd(i,1),intadd(i,2)),gplus(i,1)*fac2
rnorm2=rnorm2+gplus(i,1)*gplus(i,1)*fac2*fac2
end do
rnorm2=dsqrt(rnorm2)
print*,'Internal Gradient Norm =',rnorm2
70 format(A,f12.5,f15.5)
if(rnorm2.lt.OPTTOL*.52917/627.51)STOP




icol=0
do i=1,numint 
stat1=.false.
stat2=.true.
icol=icol+1
xold=q(intadd(i,1),intadd(i,2))
q(intadd(i,1),intadd(i,2))=q(intadd(i,1),intadd(i,2))+delta
call scfq(q,ninput,ff,stat1,stat2,gplus,numint)
q(intadd(i,1),intadd(i,2))=q(intadd(i,1),intadd(i,2))-2.0d0*delta
call scfq(q,ninput,ff,stat1,stat2,gminus,numint)
gplus=(gplus-gminus)/two/delta
q(intadd(i,1),intadd(i,2))=xold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)
!amass(icol)=atmass(zeff(i))
end do

! symmetrize hessian
hessian=hessian+transpose(hessian)
hessian=hessian/two
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

call eig(hessian,scrmat,i3n,i3n,0)
call matprt(scrmat,i3n,i3n,i3n,i3n)
do i=1,i3n
scrvec(i)=hessian(i,i)
print*,i,scrvec(i)
end do


!call tred3(i3n,i3n,hessian,scrvec,vec2,scrmat)
!call tql3(i3n,i3n,scrvec,vec2,scrmat,iout) ! eigenvals are in scrvec
!call matprt(scrmat,i3n,i3n,i3n,i3n)

! check eigenvals and print warning messages
ineg=0
do i=1,i3n
!if( dabs(scrvec(i)).lt.1D-4)scrvec(i)=zero
!print*,scrvec(i)
if(scrvec(i).lt.0.0d0)ineg=ineg+1
end do
PRINT*,'Force constant matrix has ',ineg,' negative eigenvalues'
PRINT*,'Modes corresponding to negative eigenvalues are'
do i=1,i3n
if(scrvec(i).lt.0.0d0)then
   write(*,*)'Eigenvalue:',scrvec(i)
   print*,'---------------------------'
   do j=1,i3n
      fac=1.0d0
      if(inttype(j).ne.'Bond')fac=180.0d0/pi
      write(*,40)intadd(j,1),' ',inttype(j),q(intadd(j,1),intadd(j,2))*fac,scrmat(j,i)
   end do
   print*,''
end if
end do
40 format(i3,A,A,f10.4,f10.4,f10.4,f10.4)


! compute F=U^T |gradient> 
!allocate(fvec(i3n))
!allocate(scale(i3n))
!allocate(step(i3n))
!convert gradient to hartrees
gint=gint/627.51d0
do j=1,i3n
if(inttype(j).eq.'Bond')gint(j)=gint(j)*.529177d0
end do
call dgemv( 'T', i3n,i3n,one,scrmat,i3n,gint,1,zero,fvec,1)


! determine lambda(p)
imode=i3n
root= scrvec(imode)*scrvec(imode)  + four * fvec(imode) * fvec(imode)
rlamp=half*scrvec(imode) + half * dsqrt(root)



allocate(scr3(i3n,i3n))
allocate(scr4(i3n,i3n))
allocate(vec3(i3n))
scr3=zero
do i=1,i3n-1
scr3(i,i)=scrvec(i)
scr3(i,i3n)=fvec(i)
scr3(i3n,i)=fvec(i)
end do

!call tred3(i3n,i3n,scr3,vec3,vec2,scrmat)
!call tql3(i3n,i3n,vec3,vec2,scr3,iout) ! eigenvals are in scrvec
call eig(scr3,scrmat
print*,'RFO step: Lambda0= ',rlamp,' LambdaN= ',vec3(1) 

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
print*,'initila length',rlength
if(rlength.gt.stpsz)then
factor=stpsz/rlength
print*,'Step size exceeds maximum (0.3). Computed step scaled by ',factor
!factor=dsqrt(factor)
step=factor*step
rlength=ddot(i3n,step,1,step,1)
rlength=dsqrt(rlength)
print*,'new length',rlength

end if

print*,'STEP VECTOR'
do i=1,i3n
if(inttype(i).eq.'Bond')step(i)=step(i)
print*,i,step(i)
end do

do i=1,nopt
     do j=1,i3n
       if(intadd(j,1).eq.internal(1,i).and. intadd(j,2).eq.internal(2,i))goto 55
     end do
55   q(  internal(1,i), internal(2,i)  ) = q(internal(1,i), internal(2,i) )+step(j)
end do

write(*,*)'GEOMETRY AT CURRENT ITERATE'
     do i=1,ninput
write(*,301)zstore(i),q(i,1),opt(i,1),q(i,2)*180.0d0/pi,opt(i,2),q(i,3)*180.0d0/pi, &
opt(i,3),ref(i,1),ref(i,2),ref(i,3)
     end do
301 format(i3,f10.5,i3,f10.5,i3,f10.5,i3,i3,i3,i3)

goto 10









!write(*,90)i,scrvec(i)
!90 format(i4,f20.4)


! eigenvals are in scrvec(i) and vecs are in scrmat








call cpusec(time2)


end subroutine make_hess2
