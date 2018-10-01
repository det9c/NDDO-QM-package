subroutine atoms_cpu(natoms)
use control
use tables 
implicit double precision (a-h,o-z)
integer,intent(in)::natoms
       include 'mpif.h'


a=float(natoms)/float(nprocs)
ieach=floor(a)

!allocate(atoms_on_cpu(nprocs))
atoms_on_cpu=ieach

j=ieach
k=mod(natoms,nprocs)
if(k.ne.0)then
j=natoms-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
atoms_on_cpu(icount)=atoms_on_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if


ifirst_atom_on_cpu(1)=1
ilast_atom_on_cpu(1)=ifirst_atom_on_cpu(1)+atoms_on_cpu(1)-1
do i=2,nprocs
ifirst_atom_on_cpu(i)=ilast_atom_on_cpu(i-1)+1
ilast_atom_on_cpu(i)=ifirst_atom_on_cpu(i)+atoms_on_cpu(i)-1
end do

if(myrank.eq.0)then
if(debug)then
print*,'atom Distribution Among Processors'
print*,'Processor          First atom            Last atom         Total'
do i=1,nprocs
print*,i,'       ',ifirst_atom_on_cpu(i),'       ',ilast_atom_on_cpu(i),atoms_on_cpu(i)
end do
end if
end if






end subroutine atoms_cpu
