! broadcast the denominator in the schmidt loop ****************

! can reduce message size of buff in the transformation of fock to MO basis and also C*fmo contribution to orbital gradient


!
!
! may not need to allocate density or pass it to build_g anymore since density is in pair_dens1 and 2
!
!  
! fix F*C loop to optimize. think about moving loop over occupieds into neighbor loop. may reduce flop count

subroutine rhf_conjugate
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


subroutine initguess_cgscf(vecs)
double precision,dimension(:,:),intent(inout)::vecs
end subroutine initguess_cgscf



end interface


double precision,dimension(:,:),allocatable::vecs,test,oldfock,error,scrmat2,density,gmat,fock,gtmp,olddens,gmat2,gtmp2,buf,scrmat3,buf2,vecs_grad,fxd,fmo,fxdmo
double precision,dimension(:,:),allocatable::step_old,grad_old
double precision,dimension(:),allocatable::vec2,work,scrvec3,occnum_on_cpu2
integer,dimension(:),allocatable::iwork,ifail,num_neighbors_cpu,ifirst_neighbor_cpu
integer,dimension(:,:),allocatable::neighbors_cpu
integer,dimension(mpi_status_size)::status,request,request2
integer,dimension(2)::ispan,ispan2
double precision,dimension(9,9)::dmat1,dmat2
integer::ineighbors_total_cpu,ineighbors_total_max
character(2)::update1
maxiter=60


update1='SD'


allocate(num_neighbors_cpu(numat))
allocate(ifirst_neighbor_cpu(numat))




if(.not.save_tree)then
if(myrank.eq.0)print*,'-----------------------------------------------'
if(myrank.eq.0)print*,'|                                             |'
if(myrank.eq.0)print*,'|    HARTREE FOCK SELF CONSISTENT FIELD       |'
if(myrank.eq.0)print*,'     Steepest Descent / Conjugate Gradient    |'
if(myrank.eq.0)print*,'-----------------------------------------------'
if(myrank.eq.0)print*,''
end if
ihigh=num_basis-(nelectrons/2)+1

!////////////////////////////////////////////////////////////////////
icols_cpu=ilastbf(ilast_atom_on_cpu(myrank+1)) - ifirstbf(ifirst_atom_on_cpu(myrank+1)) + 1
icols_cpu_max=icols_cpu
list_length=0
imax=1
call mpi_gather(icols_cpu,1,mpi_integer,list_length,1,mpi_integer,0,mpi_comm_world,ierr)
if(myrank.eq.0)then
do i=1,nprocs
 if(list_length(i).gt.icols_cpu_max)then
 icols_cpu_max=list_length(i)
 imax=i
 end if
end do
print*,'Maximum number of rows and columns is',icols_cpu_max,'on task ',imax-1
print*,'Message size for Fock matrix elements is ',icols_cpu_max*num_basis*8/1d6,' MB'
end if
call mpi_bcast(icols_cpu_max,1,mpi_integer,0,mpi_comm_world,ierr)


if(icols_cpu .gt. icols_cpu_max)then ! This can never happen, but just to be sure
   print*,'Number of Fock matrix columns exceeds max buffer length on task ',myrank
  call error_message
  stop
end if
!//////////////////////////////////////////////////////////////////







!////////////////////////////////////////////////////////////////////
list_length=0
imax=1
ineighbors_total_max=ineighbors_total
call mpi_gather(ineighbors_total,1,mpi_integer,list_length,1,mpi_integer,0,mpi_comm_world,ierr)
if(myrank.eq.0)then
do i=1,nprocs
 if(list_length(i).gt.ineighbors_total_max)then
 ineighbors_total_max=list_length(i)
 imax=i
 end if
end do
print*,'Maximum length of neighbor list is ',ineighbors_total_max,'on task ',imax-1
print*,'Message size for neighbor list is ',ineighbors_total_max*4/1d6,' MB'
end if
call mpi_bcast(ineighbors_total_max,1,mpi_integer,0,mpi_comm_world,ierr)

if(ineighbors_total .gt. ineighbors_total_max)then ! This can never happen, but just to be sure
   print*,'Number of neighbors per core exceeds max buffer length on task',myrank
  call error_message
  stop
end if
allocate(neighbors_cpu(2,ineighbors_total_max))

!//////////////////////////////////////////////////////////////////









2 allocate(density(num_basis,icols_cpu))
allocate(gmat2(num_basis,icols_cpu))
allocate(gtmp(num_basis,icols_cpu))
allocate(olddens(num_basis,icols_cpu))
allocate(fock(num_basis,icols_cpu))
allocate(gtmp2(num_basis,icols_cpu))
allocate(scrmat(num_basis,icols_cpu))
allocate(pair_dens2(denslist2_cpu))
allocate(pair_dens1(denslist1_cpu))
!allocate(scrmat2(num_basis,num_basis))
!allocate(scrmat3(num_basis,num_basis))
allocate(denslist_map2(ineighbors_total))
allocate(denslist_map1(ineighbors_total))
allocate(occnum_on_cpu2(icols_cpu_max))
!ioccupy=occupied_on_cpu(myrank+1)
! if(ioccupy.gt.100)then
!  print*,'Too many occupied orbitals per core. Increase number of cores!'
!  dj=float((nelectrons/2)/100)
!  write(*,90)'Job requires at least ',floor(dj)+1, ' tasks with current settings in rhf_conjugate.f90.'
!  call error_message
!  stop
! end if
!90 format(A,i5,A)
! if(icols_cpu.gt.100)then
!  print*,'Too many Fock-matrix columns per core. Increase number of cores!'
! dj=float(num_basis/100)
! write(*,90)'Job requires at least ',floor(dj)+1, ' tasks with current settings in rhf_conjugate.f90.'
! call error_message 
!  stop
! end if

!if(ioccupy .gt. icols_cpu_max)then ! i don't think this can ever happen either, but just to be sure
!  print*,'Number of occupied orbitals per core exceeds max buffer length '
!  call error_message
!  stop
!end if


icol=nelectrons/2
pair_dens2=0.
pair_dens1=0.
call initguess_cgscf_atomic
ioccupy=occupied_on_cpu(myrank+1)
if(ioccupy .gt. icols_cpu_max)then ! i don't think this can ever happen either, but just to be sure
  print*,'Number of occupied orbitals per core exceeds max buffer length '
  call error_message
  stop
end if


allocate(vecs(num_basis,ioccupy))
vecs=vecstt
allocate(vecs_grad(num_basis,ioccupy))
allocate(step_old(num_basis,ioccupy))
allocate(grad_old(num_basis,ioccupy))
allocate(fxd(num_basis,ioccupy))
allocate(fmo(nvecs_global,ioccupy))
!allocate(fxdmo(num_basis,ioccupy))
!2 allocate(density(ndim1))

!THIS WILL HAVE TO BE FIXED LATER ONCE WE DECIDE HOW TO HANDLE FINITE DIFFERENCE GRADIENTS

if(TRIAL)then
!density=s3
else
!call guess(density)
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
!if(allocated(scrvec))deallocate(scrvec)
!if(allocated(scrmat))deallocate(scrmat)
!allocate(scrvec(num_basis))
!allocate(scrmat(num_basis,num_basis))
!allocate(vecs(num_basis,num_basis))
!allocate(vec2(num_basis))
!allocate(scrvec3(ndim1))
!allocate(fock(ndim1))
if(DIIS)then
allocate(oldfock(ndim1,maxiter))
nlen=num_basis**2
allocate(error(nlen,maxiter))
!allocate(scrmat2(num_basis,num_basis))
end if


icol=nelectrons/2




call mpi_barrier(mpi_comm_world,ierr)
!call initguess_cgscf_atomic(vecs)

!-------------------------------
!use naive initial guess for now
!vecs=0.
!icolstart=ifirst_occ_on_cpu(myrank+1)
!icolend=ilast_occ_on_cpu(myrank+1)
!ishift=icolstart-1
!do i=icolstart,icolend
!vecs(i,i-ishift)=1.0d0
!end do
!allocate(oldfock(num_basis,num_basis))
!open(unit=1,file='vecs')
!read(1,*)oldfock
!close(1)
!icolstart=ifirst_occ_on_cpu(myrank+1)+ihigh-1
!icolend=icolstart+ioccupy-1
!ii=0
!do i=icolstart,icolend
!ii=ii+1
!do  j=1,num_basis
!vecs(j,ii)=oldfock(j,i)
!end do
!end do


iter=0
idiis=0
pdelta=zero
energy=0.
if(.not.save_tree)then
if(myrank.eq.0)write(*,*)'Beginning SCF iterations. Convergence tolerance = ',scftol
if(myrank.eq.0)write(*,*)'CG iterations begin when orbital gradient ||F|| becomes less than: ',sdtocg
if(myrank.eq.0)write(*,*)'STEP          ENERGY                    DEL. E      ELECTRONS   Time(s)   ||F||'
end if


time_begin=mpi_wtime()
timeend=time_begin
scftol=3d-5
1 iter=iter+1
timea=mpi_wtime()
call mpi_barrier(mpi_comm_world,ierr)

!if(myrank.eq.1)call matprt(density,num_basis,icols_cpu,num_basis,icols_cpu)

! build density on each core. 
!print*,'density from vecs'
density=0.0
pair_dens1=0.
pair_dens2=0.
if(nprocs.gt.1)then
 allocate(buf(num_basis,icols_cpu_max))
 allocate(buf2(num_basis,icols_cpu_max))
 buf=0
 buf2=0.
 do i=1,ioccupy
 do j=1,num_basis
 buf(j,i)=vecs(j,i)
 end do
 end do
 call mpi_barrier(mpi_comm_world,ierr)
!////////////////////////////////////////////////////////////////


 irow=0
 irow1=0

do i=1,ineighbors_total
 iatom=neighbors(1,i)
 jatom=neighbors(2,i)
 istart=ifirstbf(iatom)
 iend=ilastbf(iatom)
 jstart=ifirstbf(jatom)
 jend=ilastbf(jatom)

 denslist_map2(i)=irow+1
 irow=irow+nbas(species(iatom))*nbas(species(jatom))
 ibegin=denslist_map2(i)-1

 denslist_map1(i)=irow1+1
 irow1=irow1+nbas(species(jatom))*nbas(species(jatom))
 ibegin2=denslist_map1(i)-1


 do k=1,ioccupy ! cycles with the list
 inum=0
 inum2=0
  do jj=jstart,jend
    cjk=buf(jj,k)*occnum_on_cpu(k)
       do ii=istart,iend
       inum=inum+1
       pair_dens2(ibegin+inum)=pair_dens2(ibegin+inum)+buf(ii,k)*cjk
       end do
       do ii=jstart,jend
       inum2=inum2+1
       pair_dens1(ibegin2+inum2)=pair_dens1(ibegin2+inum2)+buf(ii,k)*cjk
       end do
  end do
end do

end do ! end neighbor loop. now cycle the vecs

call mpi_barrier(mpi_comm_world,ierr)
jdim=num_basis*icols_cpu_max

do j=1,nprocs-1 !//// cycle the vector columns
!if(myrank.eq.0)print*,'cycle',j
!time1=mpi_wtime()
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
 call mpi_wait(request,status,ierror)
 buf=buf2
 call mpi_barrier(mpi_comm_world,ierr)
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(ioccupy,1,mpi_integer,myrank+1,myrank,mpi_comm_world,request2,ierr)
 else
 call mpi_isend(ioccupy,1,mpi_integer,0,myrank,mpi_comm_world,request2,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(ioccupy2,1,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(ioccupy2,1,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
call mpi_wait(request2,status,ierror)
 ioccupy=ioccupy2
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(occnum_on_cpu,icols_cpu_max,mpi_double_precision,myrank+1,myrank,mpi_comm_world,request3,ierr)
 else
 call mpi_isend(occnum_on_cpu,icols_cpu_max,mpi_double_precision,0,myrank,mpi_comm_world,request3,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(occnum_on_cpu2,icols_cpu_max,mpi_double_precision,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(occnum_on_cpu2,icols_cpu_max,mpi_double_precision,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request3,status,ierror)
 occnum_on_cpu=occnum_on_cpu2



 !if(myrank.eq.iprinter)print*,'span',ispan,j
 call mpi_barrier(mpi_comm_world,ierr)
!time2=mpi_wtime()
do i=1,ineighbors_total
 iatom=neighbors(1,i)
 jatom=neighbors(2,i)
 istart=ifirstbf(iatom)
 iend=ilastbf(iatom)
 jstart=ifirstbf(jatom)
 jend=ilastbf(jatom)

 ibegin=denslist_map2(i)-1
 ibegin2=denslist_map1(i)-1
 do k=1,ioccupy ! cycles with the list
 inum=0
 inum2=0
  do jj=jstart,jend
    cjk=buf(jj,k)*occnum_on_cpu(k)
       do ii=istart,iend
       inum=inum+1
       pair_dens2(ibegin+inum)=pair_dens2(ibegin+inum)+buf(ii,k)*cjk
       end do
       do ii=jstart,jend
       inum2=inum2+1
       pair_dens1(ibegin2+inum2)=pair_dens1(ibegin2+inum2)+buf(ii,k)*cjk
       end do
  end do
end do

end do ! end neighbor loop


!time3=mpi_wtime()
!totaltime=time3-time1
!gx1=100.0d0*(time2-time1)/totaltime
!gx2=100.0d0*(time3-time2)/totaltime
!if(myrank.eq.0)print*,'density build com cal',gx1,gx2



  
end do

 call mpi_barrier(mpi_comm_world,ierr)

! put occupation numbers back 
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(occnum_on_cpu,icols_cpu_max,mpi_double_precision,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(occnum_on_cpu,icols_cpu_max,mpi_double_precision,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(occnum_on_cpu2,icols_cpu_max,mpi_double_precision,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(occnum_on_cpu2,icols_cpu_max,mpi_double_precision,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 occnum_on_cpu=occnum_on_cpu2
 call mpi_barrier(mpi_comm_world,ierr)





  else


 allocate(buf(num_basis,icols_cpu_max))
 allocate(buf2(num_basis,icols_cpu_max))
 buf=0
 buf2=0.
 do i=1,ioccupy
 do j=1,num_basis
 buf(j,i)=vecs(j,i)
 end do
 end do
 call mpi_barrier(mpi_comm_world,ierr)
!////////////////////////////////////////////////////////////////
 irow=0
 irow1=0
do i=1,ineighbors_total
 iatom=neighbors(1,i)
 jatom=neighbors(2,i)
 istart=ifirstbf(iatom)
iend=ilastbf(iatom)
 jstart=ifirstbf(jatom)
 jend=ilastbf(jatom)

 denslist_map2(i)=irow+1
 irow=irow+nbas(species(iatom))*nbas(species(jatom))
ibegin=denslist_map2(i)-1

 denslist_map1(i)=irow1+1
 irow1=irow1+nbas(species(jatom))*nbas(species(jatom))
ibegin2=denslist_map1(i)-1


do k=1,ioccupy ! cycles with the list
 inum=0
 inum2=0.
  do jj=jstart,jend
    cjk=buf(jj,k)*occnum_on_cpu(k)
       do ii=istart,iend
       inum=inum+1
       pair_dens2(ibegin+inum)=pair_dens2(ibegin+inum)+buf(ii,k)*cjk
       end do
       do ii=jstart,jend
       inum2=inum2+1
       pair_dens1(ibegin2+inum2)=pair_dens1(ibegin2+inum2)+buf(ii,k)*cjk
       end do
  end do
end do

end do

end if


!pair_dens2= !2.*pair_dens2
!pair_dens1= !2.*pair_dens1

timeb=mpi_wtime()
!if(myrank.eq.0)print*,'density time is',timeb-timea
!do i=1,ineighbors_total
!if(myrank.eq.iprinter)print*,'map',neighbors(1,i),neighbors(2,i),denslist_map1(i)
!end do

total_elec=0
particles=0.
do i=1,ineighbors_total
 iatom=neighbors(1,i)
 jatom=neighbors(2,i)
 if(iatom.eq.jatom)then
 istart=ifirstbf(iatom)
 iend=ilastbf(iatom)
 ibegin=denslist_map2(i)-1
 inum=0
 onhere=0
       do ii=istart,iend
       do jj=istart,iend
       inum=inum+1
       if(ii.eq.jj)onhere=onhere+pair_dens2(ibegin+inum)
       if(ii.eq.jj)particles=particles+pair_dens2(ibegin+inum)
       end do
       end do
!       print*,'atom',iatom,'has ',onhere,' electrons'
  end if
end do ! end neighbor loop

call mpi_reduce(particles,total_elec,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'Trace of distributed density matrix: ',total_elec,' electrons'

timec=mpi_wtime()
!if(myrank.eq.0)print*,'trace time is',timec-timeb


call mpi_barrier(mpi_comm_world,ierr)
call build_g_se(density,fock,gtmp,scrmat,gtmp2)
call mpi_barrier(mpi_comm_world,ierr)



timed=mpi_wtime()
!if(myrank.eq.0)print*,'fock build',timed-timec


! compute energy
icounter=0
ishift2=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
energylocal=0.0
do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)
 do kk=1,num_neighbors(i)
 jspot=ifirst_neighbor(i)+kk-1
 j=neighbors(2,jspot)
 icounter=icounter+1
 istart=ifirstbf(i)
 iend=ilastbf(i)
 jstart=ifirstbf(j)
 jend=ilastbf(j)
 idens_row2=denslist_map2(icounter)-1
 ni=nbas(species(i))
 nj=nbas(species(j))
 ll=0
 do i1=1,nj
 do i2=1,ni
 ll=ll+1
 dmat2(i1,i2)=pair_dens2(idens_row2+ll)
 end do
 end do
 ishift=ifirstbf(i)-1
 jshift=ifirstbf(j)-1
 inum=0
 do ii=istart,iend
 do jj=jstart,jend
   energylocal=energylocal+dmat2(jj-jshift,ii-ishift )*scrmat(jj,ii-ishift2)
 end do
 end do
 end do
end do
energylocal=energylocal*half
old=energy
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(energylocal,energy,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'is this it?',energy
delta=abs(energy-old)
timestart=timeend
if(.not.save_tree)then
timeend=mpi_wtime()
if(myrank.eq.0)write(*,40)iter,energy,delta,total_elec,timeend-timestart,rnorm_init,update1
end if
40 format(i5,f25.10,f25.10,f10.2,f10.2,f15.6,A3)
call mpi_bcast(delta,1,mpi_double_precision,0,mpi_comm_world,ierr)
if(delta.lt.scftol)goto 20
timee=mpi_wtime()
!if(myrank.eq.0)print*,'energy time',timee-timed



call mpi_barrier(mpi_comm_world,ierr)




!//////////////////////////////////////////////////////////////////////////
! compute orbital gradient F*C by cycling fock
vecs_grad=0.
ioccupy=occupied_on_cpu(myrank+1)
 ispan(1)=ifirst_atom_on_cpu(myrank+1)
 ispan(2)=ilast_atom_on_cpu(myrank+1)
if(nprocs.gt.1)then
 buf=0.
 buf2=0.
 do i=1,icols_cpu
 do j=1,num_basis
 buf(j,i)=fock(j,i)
 end do
 end do
icolstart=ifirstbf(ifirst_atom_on_cpu(myrank+1))
icolend=ilastbf(ilast_atom_on_cpu(myrank+1))
ishift=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
ioccupy=occupied_on_cpu(myrank+1)
!----------------------computes gradient of own vecs
do k=1,ioccupy
  do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)
    do irow_basis=ifirstbf(i),ilastbf(i)
       vecs_grad(irow_basis,k)=0.
       do kk=1,num_neighbors(i)
          jspot=ifirst_neighbor(i)+kk-1
          j=neighbors(2,jspot)
          do jrow_basis=ifirstbf(j),ilastbf(j)
          vecs_grad(irow_basis,k)=vecs_grad(irow_basis,k)+fock(jrow_basis,irow_basis-ishift)*vecs(jrow_basis,k)
          end do
       end do
     end do
   end do
end do


!--------------------done computing gradient of own vecs. now pass the fock vecs and continue
call mpi_barrier(mpi_comm_world,ierr)
do j=1,nprocs-1 !//// cycle the fock columns
!time1=mpi_wtime()
jdim=num_basis*icols_cpu_max
!if(myrank.eq.0)then
! ispeedy=0
! do i=1,icols_cpu
! do jh=1,num_basis
! if(dabs(buf(jh,i)) .gt. 1d-12)ispeedy=ispeedy+1
! end do
! end do
! print*,'out of',jdim,'numbers only',ispeedy,'matter',myrank
! print*,'that is only',100.0*ispeedy/jdim,'% bro',ispeedy*8/1d6
!end if
 




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
 call mpi_wait(request2,status,ierror)
 ispan=ispan2
 call mpi_barrier(mpi_comm_world,ierr)


!!!!!!!!!!!!!!!!!  FIX THIS HUGE SEND AFTER WE GET EVERYTHING IN ORDER

 jdim=2*ineighbors_total_max
 if( myrank+1 .ne. nprocs)then
! call mpi_isend(neighbors(1:2,1:ineighbors_total_max),jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
! else
! call mpi_isend(neighbors(1:2,1:ineighbors_total_max),jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 call mpi_isend(neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)


 end if
 if( myrank .ne. 0)then
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
neighbors(1:2,1:ineighbors_total_max)=neighbors_cpu

 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(num_neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(num_neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 num_neighbors=num_neighbors_cpu
 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 ifirst_neighbor=ifirst_neighbor_cpu
!time2=mpi_wtime()
 !if(myrank.eq.iprinter)print*,'span',ispan,j
 call mpi_barrier(mpi_comm_world,ierr)
!////////////////////////////////////////////////////////////////
icolstart=ifirstbf(ispan(1))
ishift=icolstart-1
do k=1,ioccupy
  do i=ispan(1),ispan(2)
    do irow_basis=ifirstbf(i),ilastbf(i)
       vecs_grad(irow_basis,k)=0.
       do kk=1,num_neighbors(i)
          jspot=ifirst_neighbor(i)+kk-1
          jj=neighbors(2,jspot)
          do jrow_basis=ifirstbf(jj),ilastbf(jj)
          vecs_grad(irow_basis,k)=vecs_grad(irow_basis,k)+buf(jrow_basis,irow_basis-ishift)*vecs(jrow_basis,k)
          end do
       end do
     end do
   end do
end do
 call mpi_barrier(mpi_comm_world,ierr)


!time3=mpi_wtime()
!totaltime=time3-time1
!gx1=100.0d0*(time2-time1)/totaltime
!gx2=100.0d0*(time3-time2)/totaltime
!if(myrank.eq.0)print*,'FC grad com cal',gx1,gx2

end do ! end loop over cpu




! put everything back on correct task before proceeding
 jdim=2*ineighbors_total_max
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 neighbors(1:2,1:ineighbors_total_max)=neighbors_cpu
 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(num_neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(num_neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 num_neighbors=num_neighbors_cpu
 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 ifirst_neighbor=ifirst_neighbor_cpu



else

ishift=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
ioccupy=occupied_on_cpu(myrank+1)
!----------------------computes density of own vecs

do k=1,ioccupy

  do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)

    do irow_basis=ifirstbf(i),ilastbf(i)
       vecs_grad(irow_basis,k)=0. 

       do kk=1,num_neighbors(i)
          jspot=ifirst_neighbor(i)+kk-1
          j=neighbors(2,jspot)
          do jrow_basis=ifirstbf(j),ilastbf(j)
          vecs_grad(irow_basis,k)=vecs_grad(irow_basis,k)+fock(jrow_basis,irow_basis-ishift)*vecs(jrow_basis,k)
          end do
       end do   

     end do

   end do

end do




end if

timef=mpi_wtime()
!if(myrank.eq.0)print*,'FC time',timef-timee


!sumt=0.
!energylocal=sum(vecs_grad)
!call mpi_barrier(mpi_comm_world,ierr)
!call mpi_reduce(energylocal,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'debug FC sum is',-2.0*sumt 



! transform fock matrix to MO basis

!
!//////////////////////////////////////////////////////////////////////////
fmo=0.
ioccupy=occupied_on_cpu(myrank+1)
 ispan(1)=ifirst_occ_on_cpu(myrank+1)
 ispan(2)=ilast_occ_on_cpu(myrank+1)
if(nprocs.gt.1)then
 buf=0.
 buf2=0.
 do i=1,ioccupy
 do j=1,num_basis
 buf(j,i)=vecs(j,i)
 end do
 end do
icolstart=ifirst_occ_on_cpu(myrank+1)
icolend=ilast_occ_on_cpu(myrank+1)
ishift=icolstart-1

do j=1,ioccupy
   do i=icolstart,icolend
   fmo(i,j)=ddot(num_basis,buf(1,i-ishift),1,vecs_grad(1,j),1)
   end do
end do


!--------------------done computing gradient of own vecs. now pass the fock vecs and continue
call mpi_barrier(mpi_comm_world,ierr)
do j=1,nprocs-1 !//// cycle the vector columns
jdim=num_basis*icols_cpu_max  ! can reduce this message size!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
 call mpi_wait(request2,status,ierror)
 ispan=ispan2
 call mpi_barrier(mpi_comm_world,ierr)

icolstart=ispan(1)
icolend=ispan(2)
ishift=icolstart-1

 do jj=1,ioccupy
   do i=icolstart,icolend
   fmo(i,jj)=ddot(num_basis,buf(1,i-ishift),1,vecs_grad(1,jj),1)
   end do
end do



 call mpi_barrier(mpi_comm_world,ierr)
end do ! end loop over cpu


else


ioccupy=occupied_on_cpu(myrank+1)

do j=1,ioccupy
   do i=1,ioccupy
   fmo(i,j)=ddot(num_basis,vecs(1,i),1,vecs_grad(1,j),1)
   end do
end do


end if

timeg=mpi_wtime()
!if(myrank.eq.0)print*,'FMo time',timeg-timef




sumt=0.
energylocal=sum(fmo)
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(energylocal,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'debug Fmo sum is',sumt



! compute overlap contribution to orbital gradient by cycling vecs C*Fmo 
!//////////////////////////////////////////////////////////////////////////

ioccupy=occupied_on_cpu(myrank+1)
 ispan(1)=ifirst_occ_on_cpu(myrank+1)
 ispan(2)=ilast_occ_on_cpu(myrank+1)
if(nprocs.gt.1)then
 buf=0.
 buf2=0.
 do i=1,ioccupy
 do j=1,num_basis
 buf(j,i)=vecs(j,i)
 end do
 end do
jcolstart=ifirst_occ_on_cpu(myrank+1)
jcolend=ilast_occ_on_cpu(myrank+1)
jshift=jcolstart-1
do k=1,ioccupy
 do ilam=1,num_basis
   total=0.0
    do j=jcolstart,jcolend
    total=total+buf(ilam,j-jshift)*fmo(j,k)
    end do
  vecs_grad(ilam,k)=vecs_grad(ilam,k)-total
 end do
end do
!--------------------done computing gradient of own vecs. now pass the fock vecs and continue
call mpi_barrier(mpi_comm_world,ierr)
do j=1,nprocs-1 !//// cycle the vector columns
 jdim=num_basis*icols_cpu_max  ! can reduce this message size!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
 call mpi_wait(request2,status,ierror)
 ispan=ispan2
 call mpi_barrier(mpi_comm_world,ierr)

jcolstart=ispan(1)
jcolend=ispan(2)
jshift=jcolstart-1

do k=1,ioccupy
 do ilam=1,num_basis
   total=0.0
    do jj=jcolstart,jcolend
    total=total+ buf(ilam,jj-jshift)*fmo(jj,k)
    end do
  vecs_grad(ilam,k)=vecs_grad(ilam,k)-total
 end do
end do





 call mpi_barrier(mpi_comm_world,ierr)
end do ! end loop over cpu

else

ioccupy=occupied_on_cpu(myrank+1)

do k=1,ioccupy
 do ilam=1,num_basis
   total=0.0
    do j=1,ioccupy
    total=total+vecs(ilam,j)*fmo(j,k)
    end do
  vecs_grad(ilam,k)=vecs_grad(ilam,k)-total
 end do
end do



end if




vecs_grad=-2.0d0*vecs_grad


jdim=num_basis*occupied_on_cpu(myrank+1)
rnorm=ddot(jdim,vecs_grad(1,1),1,vecs_grad(1,1),1)
sumt=0.
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(rnorm,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
rnorm_init=dsqrt(sumt)

!if(myrank.eq.0 )print*,'Initial gradient norm of E=Tr(PF)',rnorm_init

grad_old=vecs_grad


timeh=mpi_wtime()
!if(myrank.eq.0)print*,'overla piece to grad',timeh-timeg
imax=1
if(rnorm_init.lt. sdtocg .and. myrank.eq.0 .and. DEBUG)then
print*,'-------------------------------------------------------'
print*,'            Conjugate Gradient Microiterations    '
print*,'CG Step             Energy             ||F||        Delta E  '
print*,'-------------------------------------------------------'
end if
update1='SD'
if(rnorm_init.lt. sdtocg)then
imax=20
update1='CG'
end if
call mpi_bcast(imax,1,mpi_integer,0,mpi_comm_world,ierr)


do ii=1,imax

!//////////////////////////////////////////////////////////////////////////
! compute  F*D by cycling fock
fxd=0.
ioccupy=occupied_on_cpu(myrank+1)
 ispan(1)=ifirst_atom_on_cpu(myrank+1)
 ispan(2)=ilast_atom_on_cpu(myrank+1)
if(nprocs.gt.1)then
 buf=0.
 buf2=0.
 do i=1,icols_cpu
 do j=1,num_basis
 buf(j,i)=fock(j,i)
 end do
 end do
icolstart=ifirstbf(ifirst_atom_on_cpu(myrank+1))
icolend=ilastbf(ilast_atom_on_cpu(myrank+1))
ishift=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
ioccupy=occupied_on_cpu(myrank+1)
!----------------------computes gradient of own vecs
do k=1,ioccupy
  do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)
    do irow_basis=ifirstbf(i),ilastbf(i)
       fxd(irow_basis,k)=0.
       do kk=1,num_neighbors(i)
          jspot=ifirst_neighbor(i)+kk-1
          j=neighbors(2,jspot)
          do jrow_basis=ifirstbf(j),ilastbf(j)
          fxd(irow_basis,k)=fxd(irow_basis,k)+fock(jrow_basis,irow_basis-ishift)*vecs_grad(jrow_basis,k)
          end do
       end do
     end do
   end do
end do
!--------------------done computing gradient of own vecs. now pass the fock vecs and continue
call mpi_barrier(mpi_comm_world,ierr)
do j=1,nprocs-1 !//// cycle the fock columns
jdim=num_basis*icols_cpu_max
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
 call mpi_wait(request2,status,ierror)
 ispan=ispan2
 call mpi_barrier(mpi_comm_world,ierr)


!!!!!!!!!!!!!!!!!  FIX THIS HUGE SEND AFTER WE GET EVERYTHING IN ORDER

 jdim=2*ineighbors_total_max
 if( myrank+1 .ne. nprocs)then
! call mpi_isend(neighbors(1:2,1:ineighbors_total_max),jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
! else
! call mpi_isend(neighbors(1:2,1:ineighbors_total_max),jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 call mpi_isend(neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
neighbors(1:2,1:ineighbors_total_max)=neighbors_cpu

 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(num_neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(num_neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 num_neighbors=num_neighbors_cpu
 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 ifirst_neighbor=ifirst_neighbor_cpu

 !if(myrank.eq.iprinter)print*,'span',ispan,j
 call mpi_barrier(mpi_comm_world,ierr)
!////////////////////////////////////////////////////////////////
icolstart=ifirstbf(ispan(1))
ishift=icolstart-1
do k=1,ioccupy
  do i=ispan(1),ispan(2)
    do irow_basis=ifirstbf(i),ilastbf(i)
       fxd(irow_basis,k)=0.
       do kk=1,num_neighbors(i)
          jspot=ifirst_neighbor(i)+kk-1
          jj=neighbors(2,jspot)
          do jrow_basis=ifirstbf(jj),ilastbf(jj)
          fxd(irow_basis,k)=fxd(irow_basis,k)+buf(jrow_basis,irow_basis-ishift)*vecs_grad(jrow_basis,k)
          end do
       end do
     end do
   end do
end do
 call mpi_barrier(mpi_comm_world,ierr)
end do ! end loop over cpu

! put everything back on correct task before proceeding
 jdim=2*ineighbors_total_max
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 neighbors(1:2,1:ineighbors_total_max)=neighbors_cpu
 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(num_neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(num_neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 num_neighbors=num_neighbors_cpu
 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 ifirst_neighbor=ifirst_neighbor_cpu




else

ishift=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
ioccupy=occupied_on_cpu(myrank+1)
!----------------------computes density of own vecs

do k=1,ioccupy

  do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)

    do irow_basis=ifirstbf(i),ilastbf(i)
       fxd(irow_basis,k)=0.

       do kk=1,num_neighbors(i)
          jspot=ifirst_neighbor(i)+kk-1
          j=neighbors(2,jspot)
          do jrow_basis=ifirstbf(j),ilastbf(j)
          fxd(irow_basis,k)=fxd(irow_basis,k)+fock(jrow_basis,irow_basis-ishift)*vecs_grad(jrow_basis,k)
          end do
       end do

     end do

   end do

end do


end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!sumt=0.
!energylocal=sum(fxd)
!call mpi_barrier(mpi_comm_world,ierr)
!call mpi_reduce(energylocal,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'debug FD sum is',sumt

timei=mpi_wtime()
!if(myrank.eq.0)print*,'FD part',timei-timeh



! compute trace C+FD and D+FD

trace1=0.
trace2=0.
do k=1,ioccupy
   total1=0.0
   total2=0.0
   do i=1,num_basis
   total1=total1+vecs(i,k)*fxd(i,k)
   total2=total2+vecs_grad(i,k)*fxd(i,k)
   end do
 trace1=trace1+total1
 trace2=trace2+total2
end do



sumt=0.
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(trace1,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
if(myrank.eq.0 .and. DEBUG)print*,'aterm',sumt
aterm=sumt
sumt=0.
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(trace2,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
if(myrank.eq.0 .and. DEBUG)print*,'bterm',sumt
bterm=sumt

timej=mpi_wtime()
!if(myrank.eq.0)print*,'trace to get aterm bterm',timej-timei


! get c and d
!//////////////////////////////////////////////////////////////////////////

ioccupy=occupied_on_cpu(myrank+1)
 ispan(1)=ifirst_occ_on_cpu(myrank+1)
 ispan(2)=ilast_occ_on_cpu(myrank+1)
fxd=0.
if(nprocs.gt.1)then
 buf=0.
 buf2=0.
 do i=1,ioccupy
 do j=1,num_basis
 buf(j,i)=vecs_grad(j,i)
 end do
 end do
jcolstart=ifirst_occ_on_cpu(myrank+1)
jcolend=ilast_occ_on_cpu(myrank+1)
jshift=jcolstart-1
do k=1,ioccupy
 do ilam=1,num_basis
   total=0.0
    do j=jcolstart,jcolend
    total=total+buf(ilam,j-jshift)*fmo(j,k)
    end do
  fxd(ilam,k)=fxd(ilam,k)-total
 end do
end do
!--------------------done computing gradient of own vecs. now pass the fock vecs and continue
call mpi_barrier(mpi_comm_world,ierr)
do j=1,nprocs-1 !//// cycle the vector columns
 jdim=num_basis*icols_cpu_max  ! can reduce this message size!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
 call mpi_wait(request2,status,ierror)
 ispan=ispan2
 call mpi_barrier(mpi_comm_world,ierr)

jcolstart=ispan(1)
jcolend=ispan(2)
jshift=jcolstart-1

do k=1,ioccupy
 do ilam=1,num_basis
   total=0.0
    do jj=jcolstart,jcolend
    total=total+buf(ilam,jj-jshift)*fmo(jj,k)
    end do
  fxd(ilam,k)=fxd(ilam,k)-total
 end do
end do

 call mpi_barrier(mpi_comm_world,ierr)
end do ! end loop over cpu

else

ioccupy=occupied_on_cpu(myrank+1)

do k=1,ioccupy
 do ilam=1,num_basis
   total=0.0
    do j=1,ioccupy
    total=total+vecs_grad(ilam,j)*fmo(j,k)
    end do
  fxd(ilam,k)=fxd(ilam,k)-total
 end do
end do


end if


timek=mpi_wtime()
!if(myrank.eq.0)print*,'FD+mo',timek-timej



trace1=0.
trace2=0.
do k=1,ioccupy
   total1=0.0
   total2=0.0
   do i=1,num_basis
   total1=total1+vecs(i,k)*fxd(i,k)
   total2=total2+vecs_grad(i,k)*fxd(i,k)
   end do
 trace1=trace1+total1
 trace2=trace2+total2
end do


sumt=0.
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(trace1,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
if(myrank.eq.0 .and. DEBUG)print*,'cterm',sumt
cterm=sumt
sumt=0.
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(trace2,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
if(myrank.eq.0 .and. DEBUG)print*,'dterm',sumt
dterm=sumt

timel=mpi_wtime()
!if(myrank.eq.0)print*,'trace to get cterm dterm',timel-timek


if(myrank.eq.0)then
aterm=aterm*2.
bterm=bterm*2.
cterm=cterm*2.
dterm=dterm*2.
aterm=aterm+cterm
bterm=bterm+dterm
rlambda=-aterm/bterm
if(DEBUG)print*,'lambda',rlambda
end if
call mpi_bcast(rlambda,1,mpi_double_precision,0,mpi_comm_world,ierr)

!vecs=vecs+rlambda*vecs_grad

jdim=num_basis*occupied_on_cpu(myrank+1)
CALL DAXPY(jdim,rlambda,vecs_grad(1,1),1,vecs(1,1),1)


call mpi_barrier(mpi_comm_world,ierr)
if(imax.eq.1)exit



step_old=vecs_grad ! copy old step


!//////////////////////////////////////////////////////////////////////////
! compute orbital gradient F*C by cycling fock
vecs_grad=0.
ioccupy=occupied_on_cpu(myrank+1)
 ispan(1)=ifirst_atom_on_cpu(myrank+1)
 ispan(2)=ilast_atom_on_cpu(myrank+1)
if(nprocs.gt.1)then
 buf=0.
 buf2=0.
 do i=1,icols_cpu
 do j=1,num_basis
 buf(j,i)=fock(j,i)
 end do
 end do
icolstart=ifirstbf(ifirst_atom_on_cpu(myrank+1))
icolend=ilastbf(ilast_atom_on_cpu(myrank+1))
ishift=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
ioccupy=occupied_on_cpu(myrank+1)
!----------------------computes gradient of own vecs
do k=1,ioccupy
  do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)
    do irow_basis=ifirstbf(i),ilastbf(i)
       vecs_grad(irow_basis,k)=0.
       do kk=1,num_neighbors(i)
          jspot=ifirst_neighbor(i)+kk-1
          j=neighbors(2,jspot)
          do jrow_basis=ifirstbf(j),ilastbf(j)
          vecs_grad(irow_basis,k)=vecs_grad(irow_basis,k)+fock(jrow_basis,irow_basis-ishift)*vecs(jrow_basis,k)
          end do
       end do
     end do
   end do
end do


!--------------------done computing gradient of own vecs. now pass the fock vecs and continue
call mpi_barrier(mpi_comm_world,ierr)
do j=1,nprocs-1 !//// cycle the fock columns
!time1=mpi_wtime()
jdim=num_basis*icols_cpu_max
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
 call mpi_wait(request2,status,ierror)
 ispan=ispan2
 call mpi_barrier(mpi_comm_world,ierr)


!!!!!!!!!!!!!!!!!  FIX THIS HUGE SEND AFTER WE GET EVERYTHING IN ORDER

 jdim=2*ineighbors_total_max
 if( myrank+1 .ne. nprocs)then
! call mpi_isend(neighbors(1:2,1:ineighbors_total_max),jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
! else
! call mpi_isend(neighbors(1:2,1:ineighbors_total_max),jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 call mpi_isend(neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)


 end if
 if( myrank .ne. 0)then
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
neighbors(1:2,1:ineighbors_total_max)=neighbors_cpu

 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(num_neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(num_neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 num_neighbors=num_neighbors_cpu
 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 ifirst_neighbor=ifirst_neighbor_cpu
!time2=mpi_wtime()
 !if(myrank.eq.iprinter)print*,'span',ispan,j
 call mpi_barrier(mpi_comm_world,ierr)
!////////////////////////////////////////////////////////////////
icolstart=ifirstbf(ispan(1))
ishift=icolstart-1
do k=1,ioccupy
  do i=ispan(1),ispan(2)
    do irow_basis=ifirstbf(i),ilastbf(i)
       vecs_grad(irow_basis,k)=0.
       do kk=1,num_neighbors(i)
          jspot=ifirst_neighbor(i)+kk-1
          jj=neighbors(2,jspot)
          do jrow_basis=ifirstbf(jj),ilastbf(jj)
          vecs_grad(irow_basis,k)=vecs_grad(irow_basis,k)+buf(jrow_basis,irow_basis-ishift)*vecs(jrow_basis,k)
          end do
       end do
     end do
   end do
end do
 call mpi_barrier(mpi_comm_world,ierr)


!time3=mpi_wtime()
!totaltime=time3-time1
!gx1=100.0d0*(time2-time1)/totaltime
!gx2=100.0d0*(time3-time2)/totaltime
!if(myrank.eq.0)print*,'FC grad com cal',gx1,gx2

end do ! end loop over cpu




! put everything back on correct task before proceeding
 jdim=2*ineighbors_total_max
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 neighbors(1:2,1:ineighbors_total_max)=neighbors_cpu
 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(num_neighbors,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(num_neighbors,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(num_neighbors_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 num_neighbors=num_neighbors_cpu
 call mpi_barrier(mpi_comm_world,ierr)
 jdim=numat
 if( myrank+1 .ne. nprocs)then
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,myrank+1,myrank,mpi_comm_world,request,ierr)
 else
 call mpi_isend(ifirst_neighbor,jdim,mpi_integer,0,myrank,mpi_comm_world,request,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(ifirst_neighbor_cpu,jdim,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request,status,ierror)
 ifirst_neighbor=ifirst_neighbor_cpu



else

ishift=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
ioccupy=occupied_on_cpu(myrank+1)
!----------------------computes density of own vecs

do k=1,ioccupy

  do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)

    do irow_basis=ifirstbf(i),ilastbf(i)
       vecs_grad(irow_basis,k)=0.

       do kk=1,num_neighbors(i)
          jspot=ifirst_neighbor(i)+kk-1
          j=neighbors(2,jspot)
          do jrow_basis=ifirstbf(j),ilastbf(j)
          vecs_grad(irow_basis,k)=vecs_grad(irow_basis,k)+fock(jrow_basis,irow_basis-ishift)*vecs(jrow_basis,k)
          end do
       end do

     end do

   end do

end do

end if



ioccupy=occupied_on_cpu(myrank+1)
 ispan(1)=ifirst_occ_on_cpu(myrank+1)
 ispan(2)=ilast_occ_on_cpu(myrank+1)
if(nprocs.gt.1)then
 buf=0.
 buf2=0.
 do i=1,ioccupy
 do j=1,num_basis
 buf(j,i)=vecs(j,i)
 end do
 end do
jcolstart=ifirst_occ_on_cpu(myrank+1)
jcolend=ilast_occ_on_cpu(myrank+1)
jshift=jcolstart-1
do k=1,ioccupy
 do ilam=1,num_basis
   total=0.0
    do j=jcolstart,jcolend
    total=total+buf(ilam,j-jshift)*fmo(j,k)
    end do
  vecs_grad(ilam,k)=vecs_grad(ilam,k)-total
 end do
end do
!--------------------done computing gradient of own vecs. now pass the fock vecs and continue
call mpi_barrier(mpi_comm_world,ierr)
do j=1,nprocs-1 !//// cycle the vector columns
 jdim=num_basis*icols_cpu_max  ! can reduce this message size!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
 call mpi_wait(request2,status,ierror)
 ispan=ispan2
 call mpi_barrier(mpi_comm_world,ierr)

jcolstart=ispan(1)
jcolend=ispan(2)
jshift=jcolstart-1

do k=1,ioccupy
 do ilam=1,num_basis
   total=0.0
    do jj=jcolstart,jcolend
    total=total+ buf(ilam,jj-jshift)*fmo(jj,k)
    end do
  vecs_grad(ilam,k)=vecs_grad(ilam,k)-total
 end do
end do





 call mpi_barrier(mpi_comm_world,ierr)
end do ! end loop over cpu

else

ioccupy=occupied_on_cpu(myrank+1)

do k=1,ioccupy
 do ilam=1,num_basis
   total=0.0
    do j=1,ioccupy
    total=total+vecs(ilam,j)*fmo(j,k)
    end do
  vecs_grad(ilam,k)=vecs_grad(ilam,k)-total
 end do
end do



end if


vecs_grad=-2.0d0*vecs_grad !new H

jdim=num_basis*occupied_on_cpu(myrank+1)
rnorm=ddot(jdim,vecs_grad(1,1),1,vecs_grad(1,1),1)
sumt=0.
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(rnorm,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
rnorm=dsqrt(sumt)

if(myrank.eq.0 .and. DEBUG)write(*,33)ii,'      ',0.,'     ',rnorm,'      ',0.
33 format(i3,A,f20.8,A,f20.8,A,f20.8)
call mpi_bcast(rnorm,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_barrier(mpi_comm_world,ierr)
if(rnorm.lt.1d-3)exit

fxd=vecs_grad ! use fxd as scratch3

!vecs_grad=vecs_grad-grad_old
CALL DAXPY(jdim,-1.0d0,grad_old(1,1),1,vecs_grad(1,1),1)
top=ddot(jdim,vecs_grad(1,1),1,fxd(1,1),1)
bottom=ddot(jdim,grad_old(1,1),1,grad_old(1,1),1)
grad_old=fxd
sumt=0.
call mpi_barrier(mpi_comm_world,ierr)
call mpi_reduce(top,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
top=sumt
sumt=0.
call mpi_reduce(bottom,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
bottom=sumt
if(myrank.eq.0)then
rgamma=top/bottom
!print*,'gamma',rgamma
end if
call mpi_barrier(mpi_comm_world,ierr)
call mpi_bcast(rgamma,1,mpi_double_precision,0,mpi_comm_world,ierr)
!vecs_grad=fxd+rgamma*step_old
CALL DAXPY(jdim,rgamma,step_old(1,1),1,fxd(1,1),1)
call dcopy(jdim,fxd(1,1),1,vecs_grad(1,1),1)
call mpi_barrier(mpi_comm_world,ierr)


end do ! end loop over conjugate gradient iterations


if(myrank.eq.0)then
!print*,'schmidt'
vecs_grad(1,1)=ddot(num_basis,vecs(1,1),1,vecs(1,1),1)
!schmidt
do i=2,ioccupy
do j=1,i-1
top=ddot(num_basis,vecs(1,j),1,vecs(1,i),1)
bottom=vecs_grad(j,1)!ddot(num_basis,vecs(1,j),1,vecs(1,j),1)
top=-top/bottom
CALL DAXPY(num_basis,top,vecs(1,j),1,vecs(1,i),1)
end do
vecs_grad(i,1)=ddot(num_basis,vecs(1,i),1,vecs(1,i),1)
end do
!print*,'normalize'
!normalize
!do i=1,ioccupy
!factor=1.0d0/dsqrt(vecs_grad(i,1))
!call dscal(num_basis,factor,vecs(1,i),1)
!end do
!print*,'orthogonal'
buf(1:num_basis,1:ioccupy)=vecs
!do i=1,ioccupy
!do j=1,ioccupy
!top=ddot(num_basis,vecs(1,j),1,vecs(1,i),1)
!print*,j,i,top
!end do
!end do
end if 
call mpi_barrier(mpi_comm_world,ierr)
jdim=num_basis*icols_cpu_max
call mpi_bcast(buf,jdim,mpi_double_precision,0,mpi_comm_world,ierr)
ioccupy_node=ioccupy
call mpi_bcast(ioccupy_node,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_barrier(mpi_comm_world,ierr)
do jj=1,nprocs-1

if(myrank.ge.jj)then
 do i=1,ioccupy
 do j=1,ioccupy_node
 top=ddot(num_basis,buf(1,j),1,vecs(1,i),1)
 bottom=ddot(num_basis,buf(1,j),1,buf(1,j),1) ! can broadcast this too
 top=-top/bottom
 CALL DAXPY(num_basis,top,buf(1,j),1,vecs(1,i),1)
 end do
 end do
end if

if(myrank.eq.jj)then
vecs_grad(1,1)=ddot(num_basis,vecs(1,1),1,vecs(1,1),1)
do i=2,ioccupy
do j=1,i-1
top=ddot(num_basis,vecs(1,j),1,vecs(1,i),1)
bottom=vecs_grad(j,1)!ddot(num_basis,vecs(1,j),1,vecs(1,j),1)
top=-top/bottom
CALL DAXPY(num_basis,top,vecs(1,j),1,vecs(1,i),1)
end do
vecs_grad(i,1)=ddot(num_basis,vecs(1,i),1,vecs(1,i),1)
end do
end if
call mpi_barrier(mpi_comm_world,ierr)
buf(1:num_basis,1:ioccupy)=vecs
jdim=num_basis*icols_cpu_max
call mpi_bcast(buf,jdim,mpi_double_precision,jj,mpi_comm_world,ierr)
ioccupy_node=ioccupy
call mpi_bcast(ioccupy_node,1,mpi_integer,jj,mpi_comm_world,ierr)
call mpi_barrier(mpi_comm_world,ierr)

end do

!print*,'normalize'
!normalize
do i=1,ioccupy
factor=1.0d0/dsqrt(vecs_grad(i,1))
call dscal(num_basis,factor,vecs(1,i),1)
end do



!total=0.
!do k=1,ioccupy
!do i=1,num_basis
!total=total+vecs(i,k)
!end do
!end do
!sumt=0.
!call mpi_barrier(mpi_comm_world,ierr)
!call mpi_reduce(total,sumt,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'vex sum is scmidt',sumt











timem=mpi_wtime()
!if(myrank.eq.0)print*,'schmidt',timem-timel



!if(myrank.eq.0)print*,'total iteration time',timem-timea





call mpi_barrier(mpi_comm_world,ierr)


!/////////////////////////////////////////////////////////////////////////




if(allocated(buf))deallocate(buf)
if(allocated(buf2))deallocate(buf2)

go to 1

20 call mpi_barrier(mpi_comm_world,ierr)
etotal=energy+enuc
if(.not.save_tree)then
if(myrank.eq.0) write(*,*)'***********************************'
if(myrank.eq.0)write(*,*)'            SCF RESULTS              '
if(myrank.eq.0)print*,'Nuclear repulsion energy (eV) = ',enuc
if(myrank.eq.0)print*,'Electronic energy (eV) = ',energy
if(myrank.eq.0)print*,'SCF Total energy (eV) = ',etotal
!if(myrank.eq.0)print*,'homo energy (eV) = ',scrvec(ihigh)
!if(myrank.eq.0)print*,'lumo energy (eV) = ',scrvec(ihigh-1)
!xcall square(scrmat,density,num_basis)
!call trace(density,num_basis,out1)
!if(myrank.eq.0)print*,'Total # electrons at convergence (Tr P)= ',nint(out1)
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
call mpi_barrier(mpi_comm_world,ierr)

end subroutine rhf_conjugate





subroutine error_message
print*,''
print*,''
print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
print*,''
print*,''
print*,'####### ######  ######  ####### ######      '
print*,'#       #     # #     # #     # #     #     '
print*,'#       #     # #     # #     # #     #     '
print*,'#####   ######  ######  #     # ######      '
print*,'#       #   #   #   #   #     # #   #       '
print*,'#       #    #  #    #  #     # #    #      '
print*,'####### #     # #     # ####### #     #     '
print*,''
print*,''
print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
print*,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
end subroutine error_message




