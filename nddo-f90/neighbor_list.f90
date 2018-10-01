subroutine neighbor_list
use control
use tables
use constants

implicit double precision (a-h,o-z)
 include 'mpif.h'
integer,allocatable,dimension(:,:)::old_list
double precision,dimension(3,3)::cell_copy,sinv
 integer,dimension(3)::indx

! allocate num_neighbors,ifirst_neighbor in input.f90


!if(allocated(neighbors))deallocate(neighbors)
!allocate(neighbors(1,2))
!allocate(old_list(1,2))

irow=0


 if(pbc)then ! wrap the coordinates

         cell_copy=cell
         call migs(cell_copy,3,sinv,indx)
         do i=1,numat
          ssx=sinv(1,1)*x(i)+sinv(1,2)*y(i)+sinv(1,3)*z(i)
          ssy=sinv(2,1)*x(i)+sinv(2,2)*y(i)+sinv(2,3)*z(i)
          ssz=sinv(3,1)*x(i)+sinv(3,2)*y(i)+sinv(3,3)*z(i)
          xss=ssx-floor(ssx)
          yss=ssy-floor(ssy)
          zss=ssz-floor(ssz)
          fracs(1,i)=xss
          fracs(2,i)=yss
          fracs(3,i)=zss
!          print*,i,fracs(1,i),fracs(2,i),fracs(3,i)
         end do
!         stop
  end if






do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)
  near=0
  ifirst_neighbor(i)=irow+1
   do j=1,numat

             rxij=x(i)-x(j)
             ryij=y(i)-y(j)
             rzij=z(i)-z(j)
             
             if(pbc)then
                 sxij=fracs(1,i)-fracs(1,j)
                 syij=fracs(2,i)-fracs(2,j)
                 szij=fracs(3,i)-fracs(3,j)
                 sxij=sxij-anint(sxij)
                 syij=syij-anint(syij)
                 szij=szij-anint(szij)
                 rxij=cell(1,1)*sxij+cell(1,2)*syij+cell(1,3)*szij
                 ryij=cell(2,1)*sxij+cell(2,2)*syij+cell(2,3)*szij
                 rzij=cell(3,1)*sxij+cell(3,2)*syij+cell(3,3)*szij
             end if
         rijsq=rxij*rxij + ryij*ryij + rzij*rzij

         if(rijsq.lt.cutoff_sq)then
 !           deallocate(neighbors)
 !           allocate(neighbors(irow+1,2))
 !           neighbors(1:irow,1:2)=old_list(1:irow,1:2)
 !           neighbors(irow+1,1)=j
 !           neighbors(irow+1,2)=i
 !           deallocate(old_list)
 !           allocate(old_list(irow+1,2))
 !           old_list=neighbors
            near=near+1
            irow=irow+1
           neighbors(1,irow)=i
           neighbors(2,irow)=j

         end if
   end do
num_neighbors(i)=near
end do
ineighbors_total=irow

!if(debug)then
!if(myrank.eq.0)print*,'there are ',ineighbors_total,'pairs on the list'
!end if
if(ineighbors_total  .gt. 5000000)then
print*,'Number of neighbors per core exceeds max buffer length on task',myrank
call error_message
stop
end if

mem=0

!do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)
!print*,i,num_neighbors(i)
!end do
do i=1,ineighbors_total
mem=mem+nbas(species(neighbors(1,i)))*nbas(species(neighbors(2,i)))
end do
denslist2_cpu=mem
!dmem=dmem*8/1d6

list_length=0
call mpi_gather(mem,1,mpi_integer,list_length,1,mpi_integer,0,mpi_comm_world,ierr)

if(myrank.eq.0)then
if(debug)then
print*,''
print*,''
itotalints_all_cpu=0
totalmem_all_cpu=0
print*,'       Memory for Replicated (i|j) density on each core'
print*,'     Processor                   Total Mem (MB)'
do i=1,nprocs
print*,i,'       ',list_length(i)*8/1d6
end do
print*,'---------------------------------------------------------------'
end if
end if

mem2=0
do i=1,ineighbors_total
mem2=mem2+nbas(species(neighbors(2,i)))*nbas(species(neighbors(2,i)))
end do
denslist1_cpu=mem2
!dmem=dmem*8/1d6

list_length=0
call mpi_gather(mem2,1,mpi_integer,list_length,1,mpi_integer,0,mpi_comm_world,ierr)

if(myrank.eq.0)then
if(debug)then
print*,''
print*,''
itotalints_all_cpu=0
totalmem_all_cpu=0
print*,'       Memory for Replicated (i|i) density on each core'
print*,'     Processor                   Total Mem (MB)'
do i=1,nprocs
print*,i,'       ',list_length(i)*8/1d6
end do
print*,'---------------------------------------------------------------'
end if
end if




return
!deallocate(old_list)

!do i=1,irow
!if(myrank.eq.0)print*,i,neighbors(i,1),neighbors(i,2)
!end do


if(myrank.eq.0)then
a=float(ineighbors_total)/float(nprocs)
ieach=floor(a)
iatom_pairs_cpu=ieach
j=ieach
k=mod(ineighbors_total,nprocs)
if(k.ne.0)then
j=ineighbors_total-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
iatom_pairs_cpu(icount)=iatom_pairs_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
!         print*,'Number of CPUs  =',nprocs
!         print*,'Each cpu will handle',ieach,' ps pairs'
!         if(k.ne.0)print*,'Runaway slave will handle',j, 'ps pairs'
ifirst_pair_on_cpu(1)=1
ilast_pair_on_cpu(1)=ifirst_pair_on_cpu(1)+iatom_pairs_cpu(1)-1
do i=2,nprocs
ifirst_pair_on_cpu(i)=ilast_pair_on_cpu(i-1)+1
ilast_pair_on_cpu(i)=ifirst_pair_on_cpu(i)+iatom_pairs_cpu(i)-1
end do
print*,'pair Distribution Among Processors'
print*,'Processor          First pair            Last pair         Total'
do i=1,nprocs
print*,i,'       ',ifirst_pair_on_cpu(i),'       ',ilast_pair_on_cpu(i),ilast_pair_on_cpu(i)-ifirst_pair_on_cpu(i)+1
end do
end if
call mpi_bcast(ifirst_pair_on_cpu,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ilast_pair_on_cpu,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(iatom_pairs_cpu,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_barrier(mpi_comm_world,ierr)


!print*,myrank,ifirst_pair_on_cpu,ilast_pair_on_cpu


!do j=1,natoms
!print*,'atom j',j,num_neighbors(j)
!end do
!stop
!print*,'-----------------'
!do i=1,num_neighbors(j)
!print*,ipair_list(ifirst_neighbor(j)+i-1)
!end do
!end do

end subroutine neighbor_list
