 !
!
! can reduce size of vecs_guess/valsguess to 9 by icols_cpu since only need one center blocks 
!



subroutine initguess_cgscf_atomic
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
double precision,dimension(:,:),allocatable::vecs
integer,dimension(mpi_status_size)::status,request,request2
integer,dimension(2)::ispan,ispan2

double precision,dimension(:,:),allocatable::fock_small,vecs_small,vecs_guess
double precision,dimension(:),allocatable::vec_small,vec2_small,vals_guess,vals_all,vals_all2,occnum,occnum_local
integer,dimension(:),allocatable::iorder1,iorder2,ifirst_occ_guess,ilast_occ_guess

if(myrank.eq.0)print*,'Atomic initial guess'

ihigh=num_basis-(nelectrons/2)+1
icols_cpu=ilastbf(ilast_atom_on_cpu(myrank+1)) - ifirstbf(ifirst_atom_on_cpu(myrank+1)) + 1
allocate(density(num_basis,icols_cpu))
allocate(gmat2(num_basis,icols_cpu))
allocate(gtmp(num_basis,icols_cpu))
allocate(fock(num_basis,icols_cpu))
allocate(gtmp2(num_basis,icols_cpu))
allocate(vecs_local(num_basis,icols_cpu))
allocate(vecs_guess(num_basis,icols_cpu))
allocate(vals_guess(icols_cpu))
allocate(iorder1(num_basis))
allocate(iorder2(num_basis))
allocate(occnum(num_basis))
allocate(occnum_local(icols_cpu_max))
allocate(occnum_on_cpu(icols_cpu_max))

call guess(density)

icol=nelectrons/2

!///////////////////////////distribute density/////////////////////////////////////////
particles=0
if(nprocs.gt.1)then
 call cpusec(time2)
 allocate(buf(num_basis,icols_cpu_max))
 allocate(buf2(num_basis,icols_cpu_max))
 jdim=num_basis*icols_cpu_max
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
 call mpi_isend(ispan,2,mpi_integer,myrank+1,myrank,mpi_comm_world,request3,ierr)
 else
 call mpi_isend(ispan,2,mpi_integer,0,myrank,mpi_comm_world,request3,ierr)
 end if
 if( myrank .ne. 0)then
 call mpi_recv(ispan2,2,mpi_integer,myrank-1,myrank-1,mpi_comm_world,mpi_status_ignore,ierr)
 else
 call mpi_recv(ispan2,2,mpi_integer,nprocs-1,nprocs-1,mpi_comm_world,mpi_status_ignore,ierr)
 end if
 call mpi_wait(request3,status,ierror)
 ispan=ispan2
 

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

!deallocate(density)





icolumn=0
ishift=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
irow1=0!ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
vecs_guess=0.0
do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)

 ni=nbas(species(i))
 allocate(fock_small(ni,ni))
 allocate(vecs_small(ni,ni))
 allocate(vec_small(ni))
 allocate(vec2_small(ni))
 fock_small=0.0
 vecs_small=0.0
 vec_small=0.0
 vec2_small=0.0
 fock_small(1:ni,1:ni)=fock(ifirstbf(i):ilastbf(i),ifirstbf(i)-ishift:ilastbf(i)-ishift)
 call tred3(ni,ni,fock_small,vec_small,vec2_small,vecs_small)
 call tql3(ni,ni,vec_small,vec2_small,vecs_small,iout)
 call eigsrt(vec_small,vecs_small,ni,ni)

 do kk=1,ni
  icolumn=icolumn+1
  irow1=irow1+1
  vals_guess(irow1)=vec_small(kk)
  iii=ifirstbf(i)-1
  do jj=1,ni
     iii=iii+1
    vecs_guess(iii,icolumn)=vecs_small(jj,kk)
  end do
  end do

 deallocate(fock_small)
 deallocate(vecs_small)
 deallocate(vec_small)
 deallocate(vec2_small)

end do












if(allocated(buf))deallocate(buf)
allocate(vals_all(num_basis))
allocate(vals_all2(num_basis))

jdim=icols_cpu
jreq=myrank
if(myrank.gt.0)then
call mpi_isend(vals_guess,jdim,mpi_double_precision,0,myrank,mpi_comm_world,request,ierr)
call mpi_wait(request,status,ierror)
end if
if(myrank.eq.0)then
vals_all(1:icols_cpu)=vals_guess(1:icols_cpu)
if(nprocs.gt.1)then
icount=icols_cpu
do i=1,nprocs-1
icols=ilastbf(ilast_atom_on_cpu(i+1)) - ifirstbf(ifirst_atom_on_cpu(i+1)) + 1
allocate(buf(icols,1))
jdim=icols
call mpi_recv(buf,jdim,mpi_double_precision,i,i,mpi_comm_world,mpi_status_ignore,ierr)
do j=1,jdim
icount=icount+1
vals_all(icount)=buf(j,1)
end do
deallocate(buf)
end do
end if
end if










if(myrank.eq.0)then
!print*,'vals_all',vals_all
call indexx(num_basis,vals_all,iorder1)
call rank(num_basis,iorder1,iorder2)
!print*,iorder2
vals_all2=vals_all
call sort(num_basis,vals_all)
!print*,vals_all


ihigh=nelectrons/2
degen_tol=1d-3
efermi=vals_all(ihigh)
if(debug)print*,'efermi is',efermi,ihigh
! loop over eigenvalues and find window of degeneracy
iequal=0
ibase=10000000
itop=-100000000

do i=1,num_basis
ediff=(vals_all(i) - efermi)
ediff=abs(ediff)
if(ediff.le.degen_tol)then
iequal=iequal+1
if(i.lt.ibase)ibase=i
if(i.gt.itop)itop=i
end if
end do
if(debug)print*,'degeneracy level is',iequal,ibase,itop
!set occupation numbers
occnum=0.0
raufbau=0


do i=1,ibase-1
!occnum(i)=2.0d0
raufbau=raufbau+2.0d0
end do
dsmear=float(nelectrons)-raufbau
if(debug)print*,'there are',dsmear,'electrons left'
dslots=itop-ibase+1
if(debug)print*,'putting them in',dslots,' orbitals'
dsmear=dsmear/dslots
if(debug)print*,'each fermi orbital to hold ',dsmear,'particles'

vals_all=vals_all2
do i=1,num_basis
if(vals_all(i).lt.efermi) occnum(i)=2.0d0
 ediff=(vals_all(i) - efermi)
 ediff=abs(ediff)
 if(ediff.le.degen_tol)occnum(i)=dsmear
end do



! do i=ibase,itop
 !occnum(i)=dsmear
! end do
if(debug)then
 do i=1,num_basis
 print*,i,' has occupancy',occnum(i)
 end do
end if
end if

call mpi_bcast(occnum,num_basis,mpi_double_precision,0,mpi_comm_world,ierr)

ibegin = ifirstbf(ifirst_atom_on_cpu(myrank+1))
iend   =   ilastbf(ilast_atom_on_cpu(myrank+1))

num_occ=0
do i=ibegin,iend
if(occnum(i).ne.0.0d0)num_occ=num_occ+1
end do

list_length=0
call mpi_gather(num_occ,1,mpi_integer,list_length,1,mpi_integer,0,mpi_comm_world,ierr)

if(myrank.eq.0)then
nvecs=0
do i=1,nprocs
if(debug)print*,'there are',list_length(i),'occupieds on cpu',i-1
nvecs=nvecs+list_length(i)
end do
if(debug)print*,'there are ',nvecs,'total with nonzero occupancy'
end if
call mpi_bcast(nvecs,1,mpi_integer,0,mpi_comm_world,ierr)
nvecs_global=nvecs
a=float(nvecs)/float(nprocs)
ieach=floor(a)
occupied_on_cpu=ieach
j=ieach
k=mod(nvecs,nprocs)
if(k.ne.0)then
j=nvecs-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
occupied_on_cpu(icount)=occupied_on_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
ifirst_occ_on_cpu(1)=1
ilast_occ_on_cpu(1)=ifirst_occ_on_cpu(1)+occupied_on_cpu(1)-1
do i=2,nprocs
ifirst_occ_on_cpu(i)=ilast_occ_on_cpu(i-1)+1
ilast_occ_on_cpu(i)=ifirst_occ_on_cpu(i)+occupied_on_cpu(i)-1
end do

if(myrank.eq.0)then
if(debug)then
print*,'Occupied Orbital Distribution Among Processors'
print*,'Processor          First Orbital        Last Orbital       Total'
do i=1,nprocs
print*,i,'       ',ifirst_occ_on_cpu(i),'       ',ilast_occ_on_cpu(i),occupied_on_cpu(i)
end do
end if
end if



call mpi_bcast(list_length,nprocs,mpi_integer,0,mpi_comm_world,ierr)
nbelow=0
allocate(ifirst_occ_guess(nprocs))
allocate(ilast_occ_guess(nprocs))
ifirst_occ_guess=0
ilast_occ_guess=0
if(list_length(1).ne.0)then
ifirst_occ_guess(1)=nbelow+1
ilast_occ_guess(1)=ifirst_occ_guess(1)+list_length(1)-1
nbelow=ilast_occ_guess(1)
end if
do i=2,nprocs
if(list_length(i).ne.0)then
ifirst_occ_guess(i)=nbelow+1
ilast_occ_guess(i)=ifirst_occ_guess(i)+list_length(i)-1
nbelow=ilast_occ_guess(i)
end if
end do

if(myrank.eq.0)then
if(debug)then
print*,'Occupied Orbital Distribution Among Processors before distribution'
print*,'Processor          First Orbital        Last Orbital       Total'
do i=1,nprocs
print*,i,'       ',ifirst_occ_guess(i),'       ',ilast_occ_guess(i),list_length(i)
end do
end if
end if


icol=0
ispot=0
if(allocated(buf))deallocate(buf)
if(allocated(buf2))deallocate(buf2)
allocate(vecs(num_basis,occupied_on_cpu(myrank+1)))
do jj=0,nprocs-1
if(list_length(jj+1).eq.0)goto 2345

allocate(buf(num_basis,list_length(jj+1)))
allocate(buf2(num_basis,list_length(jj+1)))
if(myrank.eq.jj)then
  icol=0
  ibegin = ifirstbf(ifirst_atom_on_cpu(myrank+1)) - 1
  occnum_local=0
!  call matprt(vecs_guess,num_basis,icols_cpu,num_basis,icols_cpu)
  do i=1,icols_cpu
   ibegin=ibegin+1
!   print*,'occ',occnum(ibegin)
   if(occnum(ibegin).ne.0.0d0)then
   icol=icol+1
   call dcopy(num_basis,vecs_guess(1,i),1,buf(1,icol),1)
   occnum_local(icol)=occnum(ibegin)
   end if
  end do
!call matprt(buf,num_basis,list_length(jj+1),num_basis,list_length(jj+1)) 
end if
call mpi_barrier(mpi_comm_world,ierr)
jdim=num_basis*list_length(jj+1)
call mpi_bcast(buf,jdim,mpi_double_precision,jj,mpi_comm_world,ierr)
call mpi_bcast(occnum_local,list_length(jj+1),mpi_double_precision,jj,mpi_comm_world,ierr)

ilow=ifirst_occ_guess(jj+1)
ihigh=ilast_occ_guess(jj+1)
ihere=0
do k=ilow,ihigh
 ihere=ihere+1
 if(k .ge. ifirst_occ_on_cpu(myrank+1) .and. k .le. ilast_occ_on_cpu(myrank+1))then
 ispot=ispot+1
 occnum_on_cpu(ispot)=occnum_local(ihere)
 do mm=1,num_basis
 vecs(mm,ispot)=buf(mm,ihere)
 end do
 end if
end do

call mpi_barrier(mpi_comm_world,ierr)



deallocate(buf)
deallocate(buf2)
2345 continue
end do



!do jj=0,nprocs-1
!if(myrank.eq.jj)then
!print*,'vecs on',jj
!call matprt(vecs,num_basis,occupied_on_cpu(jj+1),num_basis,occupied_on_cpu(jj+1))
!end if
!call mpi_barrier(mpi_comm_world,ierr)
!end do



call mpi_barrier(mpi_comm_world,ierr)



allocate(vecstt(num_basis,occupied_on_cpu(myrank+1)))
vecstt=vecs
















deallocate(gmat2)
deallocate(gtmp)
deallocate(fock)
deallocate(gtmp2)
deallocate(vecs_local)




end subroutine initguess_cgscf_atomic
