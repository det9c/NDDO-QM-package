subroutine uhf
!
!  this is a modification of the rhf subroutine
!
!
use constants
use indices
use scratch_array
use tables
use control
implicit double precision(a-h,o-z)

interface
subroutine guess(density)
double precision,dimension(:),intent(inout)::density
end subroutine guess

subroutine f1uhf(fock,density,abdens)
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density,abdens
end subroutine f1uhf

subroutine trace(S,N,out)
double precision,dimension(:,:),intent(in)::S
integer::N
double precision,intent(out)::out
end subroutine trace

subroutine f2uhf(fock,density,abdens)
double precision,dimension(:),intent(inout)::fock
double precision,dimension(:),intent(in)::density,abdens
end subroutine f2uhf

              SUBROUTINE mmult(r,t,z,NBASIS,icol)
    DOUBLE PRECISION,dimension(:,:),intent(inout)::r,t,z
     double precision:: total
     integer::nbasis,icol
end subroutine mmult

subroutine square(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::matrix
double precision,dimension(:),intent(in)::vector
integer,intent(in)::idim1
end subroutine square

subroutine mkvector(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::matrix
double precision,dimension(:),intent(inout)::vector
integer,intent(in)::idim1
end subroutine mkvector

subroutine pulay(oldfock,error,fock,iter,ldim)
double precision,dimension(:,:),intent(in)::oldfock,error
double precision,dimension(:),intent(inout)::fock
integer,intent(in)::iter,ldim
end subroutine pulay


end interface


double precision,dimension(:,:),allocatable::vecs,scrmat2,oldfock,error
double precision,dimension(:,:),allocatable::oldfockb,errorb
double precision,dimension(:),allocatable::vec2,work,focka,density,scrvec3
! extra stuff needed for uhf
double precision,dimension(:),allocatable::fockb
integer,dimension(:),allocatable::iwork,ifail
maxiter=100

2 if(.not.save_tree)then
print*,'-----------------------------------------------'
print*,'|                                             |'
print*,'|    HARTREE FOCK SELF CONSISTENT FIELD       |'
print*,'-----------------------------------------------'
print*,'Unrestricted Hartree Fock Field to be determined...'
!ihigh=num_basis-(nelectrons/2)+1
end if

 allocate(density(ndim1))
if(TRIAL)then
density=s3
else
call guess(density)
end if



! get some memory
allocate(focka(ndim1))
allocate(fockb(ndim1))
if(allocated(adens))deallocate(adens)
if(allocated(bdens))deallocate(bdens)
allocate(adens(ndim1))
allocate(bdens(ndim1))
if(allocated(scrvec))deallocate(scrvec)
if(allocated(scrmat))deallocate(scrmat)
allocate(scrvec(num_basis))
allocate(scrmat(num_basis,num_basis))
allocate(scrmat2(num_basis,num_basis))
allocate(vecs(num_basis,num_basis))
allocate(vec2(num_basis))
allocate(scrvec3(ndim1))
if(DIIS)then
allocate(oldfock(ndim1,maxiter))
allocate(oldfockb(ndim1,maxiter))
nlen=num_basis**2
allocate(error(nlen,maxiter))
allocate(errorb(nlen,maxiter))
!allocate(scrmat2(num_basis,num_basis))
end if



!icol=nelectrons/2

! determine # of paired and unpaired electrons and max summation
! indices for alpha and beta density matrices
!print*,'Spin Multiplicity = ',int(mult)
iunpair=int(mult)-1
nclose=nelectrons-iunpair
numalpha=nclose/2+iunpair
numbeta=nclose/2
indexa=num_basis-numalpha+1
indexb=num_basis-numbeta+1
!print*,'# Alpha electrons =',numalpha
!print*,'# Beta  electrons =',numbeta

if(lapack)then 
nroota=numalpha+1 !by default compute occupied subspace and 1 virtual unless user wants more (or less)! to be implemented -> note this can be done simply by letting user define nroot from the input deck.  
nrootb=numbeta+1
if(.not.save_tree)then
write(*,*)'SCF will use incomplete diagonalizations'
write(*,*)'Program will determine',nroota,'lowest energy orbitals of alpha spin'
write(*,*)'Program will determine',nrootb,'lowest energy orbitals of beta spin'
end if
!lapack specific stuff
EVTOL = 2.0D0*DLAMCH('S')
lwork=8*num_basis
jwork=5*num_basis
else
if(.not.save_tree)then
write(*,*)'SCF will determine ALL eigenvectors of occupied/virtual spaces.'
write(*,*)'Note that the virtual space in not needed in SCF level calculations.'
write(*,*)'If matrices are large it is useful to set Lapack keyword in input'
write(*,*)'file which restricts diagonalization to the occupied subspace'
end if
end if


adens=density/2
bdens=density/2
! put some arbitrary noise into the alpha and beta density matrices
! i got this from thiel's code.  i have no idea why it is done this
! way.
if(numalpha.ne.numbeta)then
adens=adens*two*float(numalpha)/(numalpha+numbeta)
bdens=bdens*two*float(numbeta)/(numalpha+numbeta)
!if(numalpha==numbeta)then
!da=0.98d0
!db=two-da
!do i=1,num_basis
!dc=da
!da=db
!db=dc
!adens(i+offset1(i))=adens(i+offset1(i))*da
!bdens(i+offset1(i))=bdens(i+offset1(i))*db
!end do
!end if
end if
if(TRIAL)then
density=s3
adens=s3a
bdens=s3b
end if




! do the scf iterations

! put this in a more sensible spot one day
if(densityin)then
if(.not.save_tree)then
write(*,*)'Reading initial density from disk'
end if
open(unit=1,file='ALPHADENS')
read(1,*)adens
close(1)
open(unit=1,file='BETADENS')
read(1,*)bdens
close(1)
density=adens+bdens
end if


energy=zero
iter=0
idiis=0
if(.not.save_tree)then
write(*,*)'Beginning SCF iterations. Convergence tolerance = ',scftol
write(*,*)'CYCLE          ENERGY                    DELTA E                   DELTA P'
end if
! scf iteration
1 iter=iter+1
if(iter>100)then
   if(.not. DIIS)then
   print*,'********************************************************************************'
   print*,'SCF FAILED TO CONVERGE IN ',ITER,' CYCLES. SCF WILL REPEAT WITH DIIS TURNED ON'
    print*,'********************************************************************************'
   DIIS=.true.
   deallocate(focka)
   deallocate(fockb)
   deallocate(scrmat2)
   deallocate(vecs)
   deallocate(vec2)
   deallocate(scrvec3)
   deallocate(density)
   goto 2
   else
   print*,'SCF FAILED TO CONVERGE AFTER DIIS RESTART'
   STOP
   end if
!   goto 20
end if
call f1uhf(focka,density,adens)
call f2uhf(focka,density,adens)
call f1uhf(fockb,density,bdens)
call f2uhf(fockb,density,bdens)

old=energy
scrvec3=density*h
scrvec3=scrvec3+adens*focka
scrvec3=scrvec3+bdens*fockb
call square(scrmat,scrvec3,num_basis)

energy=half*sum(scrmat)
delta=abs(energy-old)
if(.not.save_tree)then
write(*,40)iter,energy,delta,pdelta
end if
40 format(i3,f25.10,f25.10,f25.10)
if(pdelta.lt.scftol.and.delta.lt.scftol)goto 20

scrvec3=density

if(DIIS .and. iter.gt.2)then
idiis=idiis+1
! compute FP-PF
call square(scrmat,focka,num_basis)
call square(vecs,adens,num_basis)
call dgemm( 'N', 'N', num_basis, num_basis, num_basis,one,vecs &
 , num_basis,scrmat , num_basis,zero,scrmat2,num_basis )
call dgemm( 'N', 'N', num_basis, num_basis, num_basis,one,scrmat &
 , num_basis,vecs , num_basis,-1.0d0,scrmat2,num_basis )
call dcopy(nlen,scrmat2(1,1),1,error(1,idiis),1)
call dcopy(ndim1,focka,1,oldfock(1,idiis),1)
if(iter.gt.4)call pulay(oldfock,error,focka,idiis,idiis+1)
! now the beta extrapolation
call square(scrmat,fockb,num_basis)
call square(vecs,bdens,num_basis)
call dgemm( 'N', 'N', num_basis, num_basis, num_basis,one,vecs &
 , num_basis,scrmat , num_basis,zero,scrmat2,num_basis )
call dgemm( 'N', 'N', num_basis, num_basis, num_basis,one,scrmat &
 , num_basis,vecs , num_basis,-1.0d0,scrmat2,num_basis )
call dcopy(nlen,scrmat2(1,1),1,errorb(1,idiis),1)
call dcopy(ndim1,fockb,1,oldfockb(1,idiis),1)
if(iter.gt.4)call pulay(oldfockb,errorb,fockb,idiis,idiis+1)
end if



! diagonalize fock matrix with lapack or tqli/tred
!note that beta must be diagonalized before alpha b/c alpha homo energy is written last
if(.not.Lapack)then 
call square(scrmat,fockb,num_basis)
call tred3(num_basis,num_basis,scrmat,scrvec,vec2,vecs)
call tql3(num_basis,num_basis,scrvec,vec2,vecs,iout)
call eigsrt(scrvec,vecs,num_basis,num_basis)
call dgemm( 'N', 'T', num_basis, num_basis, numbeta, one,vecs(1:num_basis,indexb:num_basis) &
 , num_basis,vecs(1:num_basis,indexb:num_basis) , num_basis,zero,scrmat,num_basis )
call mkvector(scrmat,bdens,num_basis)
call square(scrmat,focka,num_basis)
call tred3(num_basis,num_basis,scrmat,scrvec,vec2,vecs)
call tql3(num_basis,num_basis,scrvec,vec2,vecs,iout)
call eigsrt(scrvec,vecs,num_basis,num_basis)
call dgemm( 'N', 'T', num_basis, num_basis, numalpha, one,vecs(1:num_basis,indexa:num_basis) &
 , num_basis,vecs(1:num_basis,indexa:num_basis) , num_basis,zero,scrmat,num_basis )
call mkvector(scrmat,adens,num_basis)

else



allocate(work(lwork))
allocate(iwork(jwork))
allocate(ifail(num_basis))
call square(scrmat,fockb,num_basis)
call dsyevx('V','I','U',num_basis,scrmat,num_basis,0.,0.,1,nrootb,EVTOL, &
nevs,scrvec,vecs,num_basis,work,lwork,iwork,ifail,info)
deallocate(work)
deallocate(iwork)
deallocate(ifail)
do i=nrootb+1,num_basis
scrvec(i)=1D20
end do
call eigsrt(scrvec,vecs,num_basis,num_basis)
call dgemm( 'N', 'T', num_basis, num_basis, numbeta, one,vecs(1:num_basis,indexb:num_basis) &
 , num_basis,vecs(1:num_basis,indexb:num_basis) , num_basis,zero,scrmat,num_basis )
call mkvector(scrmat,bdens,num_basis)
allocate(work(lwork))
allocate(iwork(jwork))
allocate(ifail(num_basis))
call square(scrmat,focka,num_basis)
call dsyevx('V','I','U',num_basis,scrmat,num_basis,0.,0.,1,nroota,EVTOL, &
nevs,scrvec,vecs,num_basis,work,lwork,iwork,ifail,info)
deallocate(work)
deallocate(iwork)
deallocate(ifail)
do i=nroota+1,num_basis
scrvec(i)=1D20
end do


call eigsrt(scrvec,vecs,num_basis,num_basis)
call dgemm( 'N', 'T', num_basis, num_basis, numalpha, one,vecs(1:num_basis,indexa:num_basis) &
 , num_basis,vecs(1:num_basis,indexa:num_basis) , num_basis,zero,scrmat,num_basis )
call mkvector(scrmat,adens,num_basis)


end if





density=adens+bdens




pdelta=zero
do i=1,ndim1
pdelta=pdelta+(density(i)-scrvec3(i))**2
end do
! compute energy
go to 1

20 etotal=energy+enuc
if(.not.save_tree)then
write(*,*)'***********************************'
write(*,*)'            SCF RESULTS               '
print*,'Nuclear repulsion energy (eV) = ',enuc
print*,'Electronic energy (eV) = ',energy
print*,'SCF Total energy (eV) = ',etotal
print*,'homo energy (eV) = ',scrvec(indexa)
!print*,'alpha lumo energy (eV) = ',scrvec(ihigh-1)
call square(scrmat,density,num_basis)
call trace(scrmat,num_basis,out)
call square(scrmat,adens,num_basis)
call trace(scrmat,num_basis,out1)
call square(scrmat,bdens,num_basis)
call trace(scrmat,num_basis,out2)
print*,'Total # electrons at convergence (Tr P)= ',nint(out)
print*,'Total # alpha electrons at convergence [Tr P(alpha)] = ',nint(out1)
print*,'Total # beta electrons at convergence [Tr P(beta)] = ',nint(out2)
write(*,*)'***********************************'
end if
!call matprt2(fock,num_basis,num_basis,num_basis,num_basis,scrmat)
if(keep)then
  if(allocated(s3))deallocate(s3)
  if(allocated(s3a))deallocate(s3a)
  if(allocated(s3b))deallocate(s3b)  
  allocate(s3(ndim1))
  s3=density
  allocate(s3a(ndim1))
  s3a=adens
  allocate(s3b(ndim1))
  s3b=bdens
end if

if(densityout)then ! if true then write density to file

open(unit=1,file='ALPHADENS')
write(1,*)adens
close(1)
open(unit=1,file='BETADENS')
write(1,*)bdens
close(1)
if(.not.save_tree)then
write(*,*)'UHF alpha/beta densities have been written to disk'
end if
end if


deallocate(focka)
deallocate(fockb)
deallocate(density)
deallocate(vecs)
deallocate(vec2)
deallocate(scrvec)
deallocate(scrmat)
deallocate(scrvec3)
if(allocated(oldfock))deallocate(oldfock)
if(allocated(oldfockb))deallocate(oldfockb)
if(allocated(error))deallocate(error)
if(allocated(errorb))deallocate(errorb)


!10 format(i4,f15.8,f15.8,f15.8)

end subroutine uhf

 


