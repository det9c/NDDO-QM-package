subroutine initguess_cgscf(vecs)
use constants
use indices
use scratch_array
use tables
use control
implicit double precision(a-h,o-z)
include 'mpif.h'
interface
subroutine guess(density)
double precision,dimension(:,:),intent(inout)::density
end subroutine guess
subroutine build_g_se(density,gmat,gtmp,gmat2,gtmp2) !(density,g)
double precision,dimension(:,:),intent(in)::density
double precision,dimension(:,:),intent(inout)::gmat,gtmp,gmat2,gtmp2
end subroutine build_g_se
end interface
double precision,dimension(:,:),allocatable::density,fock,gtmp,gmat2,gtmp2,buf,buf2,vecs_local
double precision,dimension(:,:),allocatable::scrmati,scrmat2i
double precision,dimension(:),allocatable::vec2i,scrveci
double precision,dimension(:,:),intent(inout)::vecs
integer,dimension(mpi_status_size)::status,request,request2
integer,dimension(2)::ispan,ispan2
ihigh=num_basis-(nelectrons/2)+1
icols_cpu=ilastbf(ilast_atom_on_cpu(myrank+1)) - ifirstbf(ifirst_atom_on_cpu(myrank+1)) + 1
allocate(density(num_basis,icols_cpu))
allocate(gmat2(num_basis,icols_cpu))
allocate(gtmp(num_basis,icols_cpu))
allocate(fock(num_basis,icols_cpu))
allocate(gtmp2(num_basis,icols_cpu))
allocate(scrmati(num_basis,icols_cpu))
allocate(scrmat2i(num_basis,num_basis))
allocate(vecs_local(num_basis,num_basis))

call guess(density)

allocate(scrveci(num_basis))
allocate(vec2i(num_basis))
icol=nelectrons/2

!///////////////////////////distribute density/////////////////////////////////////////
particles=0
if(nprocs.gt.1)then
 call cpusec(time2)
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
 denslist_map2(i)=irow+1
  do jj=jstart,jend
  do ii=istart,iend
  irow=irow+1
  pair_dens2(irow)=buf(ii,jj-ishift)
  if(ii.eq.jj)particles=particles+buf(ii,jj-ishift)
   end do
   end do
  denslist_map1(i)=irow1+1 
  do jj=jstart,jend
  do ii=jstart,jend
  irow1=irow1+1
  pair_dens1(irow1)=buf(ii,jj-ishift)
  end do
  end do

  end if
end do
call mpi_barrier(mpi_comm_world,ierr)

do j=1,nprocs-1 !//// cycle the density columns   
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
 ispan=ispan2
 call mpi_wait(request2,status,ierror)
 call mpi_barrier(mpi_comm_world,ierr)
 do i=1,ineighbors_total
 iatom=neighbors(1,i)
 jatom=neighbors(2,i)
 istart=ifirstbf(iatom)
 iend=ilastbf(iatom)
 if(jatom.ge.ispan(1) .and. jatom.le.ispan(2) )then
 jstart=ifirstbf(jatom)
 jend=ilastbf(jatom)
 ishift=ifirstbf(ispan(1))-1
 denslist_map2(i)=irow+1
  do jj=jstart,jend
  do ii=istart,iend
  irow=irow+1
  pair_dens2(irow)=buf(ii,jj-ishift)
  if(ii.eq.jj)particles=particles+buf(ii,jj-ishift)
  end do
  end do
  denslist_map1(i)=irow1+1
  do jj=jstart,jend
  do ii=jstart,jend
  irow1=irow1+1
  pair_dens1(irow1)=buf(ii,jj-ishift)
  end do
  end do
  end if
 end do
 call mpi_barrier(mpi_comm_world,ierr)
end do

else  ! if only one core
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
  do jj=jstart,jend
  do ii=istart,iend
  irow=irow+1
  pair_dens2(irow)=density(ii,jj)
 if(ii.eq.jj)particles=particles+density(ii,jj)
  end do
  end do
denslist_map1(i)=irow1+1
  do jj=jstart,jend
  do ii=jstart,jend
  irow1=irow1+1
  pair_dens1(irow1)=density(ii,jj)
  end do
  end do
end do
end if
!////////////////////////////// end cycle density

call mpi_barrier(mpi_comm_world,ierr)
call build_g_se(density,fock,gtmp,scrmat,gtmp2)
call mpi_barrier(mpi_comm_world,ierr)


! collect fock matrix at root
if(allocated(buf))deallocate(buf)
jdim=icols_cpu*num_basis
jreq=myrank
if(myrank.gt.0)then
call mpi_isend(fock,jdim,mpi_double_precision,0,myrank,mpi_comm_world,request,ierr)
call mpi_wait(request,status,ierror)
end if
if(myrank.eq.0)then
scrmat2i(1:num_basis,1:icols_cpu)=fock(1:num_basis,1:icols_cpu)
if(nprocs.gt.1)then
do i=1,nprocs-1
icols=ilastbf(ilast_atom_on_cpu(i+1)) - ifirstbf(ifirst_atom_on_cpu(i+1)) + 1
allocate(buf(num_basis,icols))
jdim=num_basis*icols
call mpi_recv(buf,jdim,mpi_double_precision,i,i,mpi_comm_world,mpi_status_ignore,ierr)
scrmat2i(1:num_basis,ifirstbf(ifirst_atom_on_cpu(i+1)):ilastbf(ilast_atom_on_cpu(i+1)))=buf
deallocate(buf)
end do
end if
end if
if(myrank .eq. 0)then
call tred3(num_basis,num_basis,scrmat2i,scrveci,vec2i,vecs_local)
call tql3(num_basis,num_basis,scrveci,vec2i,vecs_local,iout)
call eigsrt(scrveci,vecs_local,num_basis,num_basis)
end if
call mpi_barrier(mpi_comm_world,ierr)

! distribute vectors
vecs=0.
jdim=num_basis*num_basis
call mpi_bcast(vecs_local,jdim,mpi_double_precision,0,mpi_comm_world,ierr)
ioccupy=occupied_on_cpu(myrank+1)
icolstart=ifirst_occ_on_cpu(myrank+1)+ihigh-1
icolend=icolstart+ioccupy-1
ii=0
do i=icolstart,icolend
ii=ii+1
do  j=1,num_basis
vecs(j,ii)=vecs_local(j,i)
end do
end do

deallocate(density)
deallocate(gmat2)
deallocate(gtmp)
deallocate(fock)
deallocate(gtmp2)
deallocate(scrmati)
deallocate(scrmat2i)
deallocate(vecs_local)




end subroutine initguess_cgscf
