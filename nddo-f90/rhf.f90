subroutine rhf
use constants
use indices
use scratch_array
use tables
use control
implicit double precision(a-h,o-z)
include 'mpif.h'
interface
subroutine mkvector(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::matrix
double precision,dimension(:),intent(inout)::vector
integer,intent(in)::idim1
end subroutine mkvector

subroutine guess(density)
double precision,dimension(:,:),intent(inout)::density
end subroutine guess

subroutine trace(S,N,out)
double precision,dimension(:,:),intent(in)::S
integer::N
double precision,intent(out)::out
end subroutine trace

subroutine square(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::matrix
double precision,dimension(:),intent(in)::vector
integer,intent(in)::idim1
end subroutine square

              SUBROUTINE mmult(r,t,z,NBASIS,icol)
    DOUBLE PRECISION,dimension(:,:),intent(inout)::r,t,z
     double precision:: total
     integer::nbasis,icol
end subroutine mmult

subroutine pulay(oldfock,error,fock,iter,ldim)
double precision,dimension(:,:),intent(in)::oldfock,error
double precision,dimension(:),intent(inout)::fock
integer,intent(in)::iter,ldim
end subroutine pulay


subroutine build_g_se(density,gmat,gtmp,gmat2,gtmp2) !(density,g)
double precision,dimension(:,:),intent(in)::density
double precision,dimension(:,:),intent(inout)::gmat,gtmp,gmat2,gtmp2
end subroutine build_g_se


end interface


double precision,dimension(:,:),allocatable::vecs,test,oldfock,error,scrmat2,density,gmat,fock,gtmp,olddens,gmat2,gtmp2,buf,scrmat3,buf2
double precision,dimension(:),allocatable::vec2,work,scrvec3
integer,dimension(:),allocatable::iwork,ifail
integer,dimension(mpi_status_size)::status,request,request2
integer,dimension(2)::ispan,ispan2
maxiter=60



if(.not.save_tree)then
if(myrank.eq.0)print*,'-----------------------------------------------'
if(myrank.eq.0)print*,'|                                             |'
if(myrank.eq.0)print*,'|    HARTREE FOCK SELF CONSISTENT FIELD       |'
if(myrank.eq.0)print*,'-----------------------------------------------'
if(myrank.eq.0)print*,''
end if
ihigh=num_basis-(nelectrons/2)+1


icols_cpu=ilastbf(ilast_atom_on_cpu(myrank+1)) - ifirstbf(ifirst_atom_on_cpu(myrank+1)) + 1


2 allocate(density(num_basis,icols_cpu))
allocate(gmat2(num_basis,icols_cpu))
allocate(gtmp(num_basis,icols_cpu))
allocate(olddens(num_basis,icols_cpu))
allocate(fock(num_basis,icols_cpu))
allocate(gtmp2(num_basis,icols_cpu))
allocate(scrmat(num_basis,icols_cpu))
allocate(pair_dens2(denslist2_cpu))
allocate(pair_dens1(denslist1_cpu))
allocate(scrmat2(num_basis,num_basis))
allocate(scrmat3(num_basis,num_basis))
allocate(denslist_map2(ineighbors_total))
allocate(denslist_map1(ineighbors_total))
!2 allocate(density(ndim1))

!THIS WILL HAVE TO BE FIXED LATER ONCE WE DECIDE HOW TO HANDLE FINITE DIFFERENCE GRADIENTS

if(TRIAL)then
!density=s3
else
call guess(density)
end if


! put this in the above if construct later
!open(unit=19,file='density')
!read(19,*)density
!close(19)



! build the diagonal fock matrix elements
if(.not.save_tree)then
if(myrank.eq.0)print*,'Restricted Hartree Fock Field to be determined...'
end if
if(densityin)then
if(.not.save_tree)then
if(myrank.eq.0)write(*,*)'Reading initial density from disk'
end if
open(unit=1,file='DENSITY')
read(1,*)density
close(1)
end if



! do the scf iterations
if(allocated(scrvec))deallocate(scrvec)
!if(allocated(scrmat))deallocate(scrmat)
allocate(scrvec(num_basis))
!allocate(scrmat(num_basis,num_basis))
allocate(vecs(num_basis,num_basis))
allocate(vec2(num_basis))
allocate(scrvec3(ndim1))
!allocate(fock(ndim1))
if(DIIS)then
allocate(oldfock(ndim1,maxiter))
nlen=num_basis**2
allocate(error(nlen,maxiter))
!allocate(scrmat2(num_basis,num_basis))
end if


icol=nelectrons/2



iter=0
idiis=0
pdelta=zero
energy=0.
if(.not.save_tree)then
if(myrank.eq.0)write(*,*)'Beginning SCF iterations. Convergence tolerance = ',scftol
if(myrank.eq.0)write(*,*)'CYCLE          ENERGY                    DELTA E                    ELECTRONS'
end if





1 iter=iter+1
call mpi_barrier(mpi_comm_world,ierr)

!if(myrank.eq.1)call matprt(density,num_basis,icols_cpu,num_basis,icols_cpu)





!///////////////////////////distribute density/////////////////////////////////////////

iprinter=-1
particles=0
if(myrank.eq.iprinter)call matprt(density,num_basis,icols_cpu,num_basis,icols_cpu)
if(nprocs.gt.1)then
 call cpusec(time2)
 if(icols_cpu.gt.100)then
  print*,'Too many matrix columns per core. Increase number of cores!',icols_cpu
  stop
 end if 
 allocate(buf(num_basis,100))
 allocate(buf2(num_basis,100))
 jdim=num_basis*100
 buf=0.
 do i=1,icols_cpu
 do j=1,num_basis
 buf(j,i)=density(j,i)
 end do
 end do 
 call mpi_barrier(mpi_comm_world,ierr)
 ispan(1)=ifirst_atom_on_cpu(myrank+1)
 ispan(2)=ilast_atom_on_cpu(myrank+1)
!////////////////////////////////////////////////////////////////
 irow=0 
 irow1=0
do i=1,ineighbors_total
 iatom=neighbors(1,i)
 jatom=neighbors(2,i)
 istart=ifirstbf(iatom)
 iend=ilastbf(iatom)
 if(jatom .ge. ispan(1) .and. jatom .le. ispan(2))then 
 jstart=ifirstbf(jatom)
 jend=ilastbf(jatom)
 ishift=ifirstbf(ispan(1))-1

 if(myrank.eq.iprinter)print*,'-------------'
 if(myrank.eq.iprinter)print*,'atoms:',iatom,jatom,jstart,jend
 denslist_map2(i)=irow+1
  do jj=jstart,jend
  do ii=istart,iend
  irow=irow+1
  pair_dens2(irow)=buf(ii,jj-ishift)
  if(ii.eq.jj)particles=particles+buf(ii,jj-ishift)
!if(myrank.eq.iprinter)  print*,irow,pair_dens2(irow),ii,jj,jj-ishift
   end do
   end do

  denslist_map1(i)=irow1+1 
  do jj=jstart,jend
  do ii=jstart,jend
  irow1=irow1+1
  pair_dens1(irow1)=buf(ii,jj-ishift)
  if(myrank.eq.iprinter)  print*,irow1,pair_dens1(irow1),ii,jj,jj-ishift
  end do
  end do

  end if
end do
!/////////////////////////////////////////
call mpi_barrier(mpi_comm_world,ierr)

if(myrank.eq.iprinter)print*,'span',ispan
 

do j=1,nprocs-1 !//// cycle the density columns   
! if(myrank.eq.0)print*,'jth shift',j

!jreq=myrank*10
if( myrank+1 .ne. nprocs)then
 call mpi_isend(buf,jdim,mpi_double_precision,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(buf,jdim,mpi_double_precision,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(buf2,jdim,mpi_double_precision,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(buf2,jdim,mpi_double_precision,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
!if(myrank.eq.0)print*,'receiving from',j
 call mpi_wait(request,status,ierror)
 buf=buf2
 call mpi_barrier(mpi_comm_world,ierr)

 if( myrank+1 .ne. nprocs)then
 call mpi_isend(ispan,2,mpi_integer,myrank+1,myrank,mpi_comm_world,request2,ierr)
 else
 call mpi_isend(ispan,2,mpi_integer,0,myrank,mpi_comm_world,request2,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(ispan2,2,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(ispan2,2,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 ispan=ispan2
 call mpi_wait(request2,status,ierror)
 !if(myrank.eq.iprinter)print*,'span',ispan,j
 call mpi_barrier(mpi_comm_world,ierr)
!////////////////////////////////////////////////////////////////
 do i=1,ineighbors_total
 iatom=neighbors(1,i)
 jatom=neighbors(2,i)
 istart=ifirstbf(iatom)
 iend=ilastbf(iatom)
 if(jatom.ge.ispan(1) .and. jatom.le.ispan(2) )then
!  do iii=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)
!  if(jatom.eq.iii )goto 80
!  end do
 jstart=ifirstbf(jatom)
 jend=ilastbf(jatom)
 ishift=ifirstbf(ispan(1))-1
 if(myrank.eq.iprinter)print*,'-------------'
 if(myrank.eq.iprinter)print*,'atoms pumper:',iatom,jatom
 denslist_map2(i)=irow+1
  do jj=jstart,jend
  do ii=istart,iend
  irow=irow+1
  pair_dens2(irow)=buf(ii,jj-ishift)
  if(ii.eq.jj)particles=particles+buf(ii,jj-ishift)
!if(myrank.eq.iprinter)  print*,irow,pair_dens2(irow),ii,jj,jj-ishift
  end do
  end do

  denslist_map1(i)=irow1+1
  do jj=jstart,jend
  do ii=jstart,jend
  irow1=irow1+1
  pair_dens1(irow1)=buf(ii,jj-ishift)
  !if(myrank.eq.iprinter)  print*,irow1,pair_dens1(irow1),ii,jj,jj-ishift
  end do
  end do
  end if
! 80 continue
 end do
!/////////////////////////////////////////

 call mpi_barrier(mpi_comm_world,ierr)
 !print*,'density after',myrank,j,buf(1,1)
end do

else  ! if only one core

irow=0
irow1=0
!print*,ineighbors_total,'damn'
do i=1,ineighbors_total
 iatom=neighbors(1,i)
 jatom=neighbors(2,i)
 istart=ifirstbf(iatom)
 iend=ilastbf(iatom)
 jstart=ifirstbf(jatom)
 jend=ilastbf(jatom)
!print*,'-------------'
!print*,'atoms:',iatom,jatom
 denslist_map2(i)=irow+1
  do jj=jstart,jend
  do ii=istart,iend
  irow=irow+1
  pair_dens2(irow)=density(ii,jj)
 if(ii.eq.jj)particles=particles+density(ii,jj)
!  print*,irow,pair_dens2(irow)
  end do
  end do

denslist_map1(i)=irow1+1
  do jj=jstart,jend
  do ii=jstart,jend
  irow1=irow1+1
  pair_dens1(irow1)=density(ii,jj)
if(myrank.eq.iprinter)  print*,irow1,pair_dens1(irow1),ii,jj,jj-ishift
  end do
  end do
end do
end if
!////////////////////////////// end cycle density


!do i=1,ineighbors_total
!if(myrank.eq.iprinter)print*,'map',neighbors(1,i),neighbors(2,i),denslist_map1(i)
!end do





total_elec=0
call mpi_reduce(particles,total_elec,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'Trace of distributed density matrix: ',total_elec,' electrons'
call mpi_barrier(mpi_comm_world,ierr)




!if(myrank.eq.0)call matprt(density,num_basis,icols_cpu,num_basis,icols_cpu)

call build_g_se(density,fock,gtmp,scrmat,gtmp2)
call mpi_barrier(mpi_comm_world,ierr)



!if(myrank.eq.0)call matprt(scrmat,num_basis,icols_cpu,num_basis,icols_cpu)

!rsum=0.
!do j=1,icols_cpu
!do i=1,num_basis
!rsum=rsum+scrmat(i,j)
!end do
!end do
!print*,'rsum',rsum,myrank

!rt=0.
!call mpi_reduce(rsum,rt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'total sum',rt




energylocal=half*sum(density*scrmat)
old=energy
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(energylocal,energy,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'is this it?',energy

call mpi_barrier(mpi_comm_world,ierr)



!//////////////////////////////////////////////////////////////////////////
! collect fock matrix at root
if(allocated(buf))deallocate(buf)

!if(myrank.gt.0)deallocate(buf)
!deallocate(scrmat)
!allocate(scrmat(num_basis,num_basis))
jdim=icols_cpu*num_basis
jreq=myrank
!if(myrank.gt.0)call mpi_isend(scrmat,jdim,mpi_double_precision,0,myrank,mpi_comm_world,jreq,ierr)
if(myrank.gt.0)then
call mpi_isend(fock,jdim,mpi_double_precision,0,myrank,mpi_comm_world,request,ierr)
call mpi_wait(request,status,ierror)
end if
if(myrank.eq.0)then
scrmat2(1:num_basis,1:icols_cpu)=fock(1:num_basis,1:icols_cpu)
if(nprocs.gt.1)then
do i=1,nprocs-1
icols=ilastbf(ilast_atom_on_cpu(i+1)) - ifirstbf(ifirst_atom_on_cpu(i+1)) + 1
allocate(buf(num_basis,icols))
jdim=num_basis*icols
call mpi_recv(buf,jdim,mpi_double_precision,i,i,mpi_comm_world,mpi_status_ignore,ierr)
scrmat2(1:num_basis,ifirstbf(ifirst_atom_on_cpu(i+1)):ilastbf(ilast_atom_on_cpu(i+1)))=buf
deallocate(buf)
end do
end if
end if

call mpi_barrier(mpi_comm_world,ierr)

!/////////////////////////////////////////////////////////////////////////



delta=abs(energy-old)
if(.not.save_tree)then
if(myrank.eq.0)write(*,40)iter,energy,delta,total_elec
end if

call mpi_barrier(mpi_comm_world,ierr)

40 format(i3,f25.10,f25.10,f25.10)
call mpi_bcast(delta,1,mpi_double_precision,0,mpi_comm_world,ierr)
if(delta.lt.scftol)goto 20

if(myrank .eq. 0)then
call tred3(num_basis,num_basis,scrmat2,scrvec,vec2,vecs)
call tql3(num_basis,num_basis,scrvec,vec2,vecs,iout)
call eigsrt(scrvec,vecs,num_basis,num_basis)
44 continue
call dgemm( 'N', 'T', num_basis, num_basis, icol, two,vecs(1:num_basis,ihigh:num_basis) &
 , num_basis,vecs(1:num_basis,ihigh:num_basis) , num_basis,zero,scrmat2,num_basis )
!open(unit=1,file='vecs')
!write(1,*)vecs
!close(1)
end if
call mpi_barrier(mpi_comm_world,ierr)


! distribute density
if(myrank.eq.0)then

 do i=1,icols_cpu
 do j=1,num_basis
 density(j,i)=scrmat2(j,i)
 end do
 end do
end if
call mpi_barrier(mpi_comm_world,ierr)
if(myrank.eq.0 .and. nprocs .gt. 1)then
do i=1,nprocs-1
icols=ilastbf(ilast_atom_on_cpu(i+1)) - ifirstbf(ifirst_atom_on_cpu(i+1)) + 1
jdim=num_basis*icols
call mpi_isend(scrmat2(1,ifirstbf(ifirst_atom_on_cpu(i+1))),jdim,mpi_double_precision,i,i,mpi_comm_world,request,ierr)
call mpi_wait(request,status,ierror)
end do
end if

if( myrank .ne. 0)then
jdim=num_basis*icols_cpu
call mpi_recv(density,jdim,mpi_double_precision,0,myrank,mpi_comm_world,mpi_status_ignore,ierr)
end if




call mpi_barrier(mpi_comm_world,ierr)





if(allocated(buf))deallocate(buf)
if(allocated(buf2))deallocate(buf2)

go to 1

20 etotal=energy+enuc
if(.not.save_tree)then
if(myrank.eq.0) write(*,*)'***********************************'
if(myrank.eq.0)write(*,*)'            SCF RESULTS              '
if(myrank.eq.0)print*,'Nuclear repulsion energy (eV) = ',enuc
if(myrank.eq.0)print*,'Electronic energy (eV) = ',energy
if(myrank.eq.0)print*,'SCF Total energy (eV) = ',etotal
if(myrank.eq.0)print*,'homo energy (eV) = ',scrvec(ihigh)
if(myrank.eq.0)print*,'lumo energy (eV) = ',scrvec(ihigh-1)
!xcall square(scrmat,density,num_basis)
!call trace(density,num_basis,out1)
if(myrank.eq.0)print*,'Total # electrons at convergence (Tr P)= ',total_elec
if(myrank.eq.0)write(*,*)'***********************************'
end if


if(keep)then
    if(allocated(s3))then
    deallocate(s3)
    end if
  allocate(s3(ndim1))
!x  s3=density
end if

if(densityout)then ! if true then write density to file

open(unit=1,file='DENSITY')
write(1,*)density
close(1)
if(.not.save_tree)then
if(myrank.eq.0)write(*,*)'RHF density has been written to disk'
end if
end if





deallocate(fock)
deallocate(density)
deallocate(vecs)
deallocate(vec2)
deallocate(scrvec)
deallocate(scrmat)
deallocate(scrvec3)
if(allocated(oldfock))deallocate(oldfock)
if(allocated(error))deallocate(error)

!10 format(i4,f15.8,f15.8,f15.8)

end subroutine rhf

     SUBROUTINE mmult(r,t,z,NBASIS,icol)
implicit double precision(a-h,o-z)
      DOUBLE PRECISION,dimension(:,:),intent(inout)::r,t,z
     double precision:: total
     integer::nbasis,icol
      do  i=1,NBASIS
         do  j=1,NBASIS
            total=0.0D0
            do  k=1,icol
               total=total+r(i,k)*t(k,j)
   end do
               z(i,j)=total
 end do
 end do
         
         end subroutine mmult



