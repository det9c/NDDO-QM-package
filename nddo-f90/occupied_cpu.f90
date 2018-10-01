subroutine occupied_cpu
use control
use tables 
use constants
implicit double precision (a-h,o-z)
       include 'mpif.h'

nhole=nelectrons/2
a=float(nhole)/float(nprocs)
ieach=floor(a)

!allocate(atoms_on_cpu(nprocs))
occupied_on_cpu=ieach

j=ieach
k=mod(nhole,nprocs)
if(k.ne.0)then
j=nhole-(nprocs)*ieach
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






end subroutine occupied_cpu
