subroutine hcore
   use constants
use tables
use indices
use control
implicit double precision (a-h,o-z)
include 'mpif.h'
double precision,dimension(:),allocatable::htemp
double precision,dimension(:,:),allocatable::wout
interface
 subroutine twoe_hcore_ints(i,j,wout)
  double precision,dimension(:,:),intent(inout)::wout
integer,intent(in)::i,j
end subroutine twoe_hcore_ints
end interface
!integer,dimension(:),allocatable::map
ir=0
jc=0
ir2=0
jc2=  0
nelectrons=0
thresh=1d-12
! loopf 

allocate(ifirsthc1c(numat))
allocate(ilasthc1c(numat))


call mpi_barrier(mpi_comm_world,ierr)


do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)
ifirsthc1c(i)=itotalh+1
!compute # of electron
nelectrons=nelectrons+eff_core(zeff(i))
   ni=nbas(species(i))
lmax=10
if(ni>4)lmax=45
!      istart=ifirst(i)
!      iend=ilast(i)
      istart=ifirstbf(i)
      iend=ilastbf(i)
!      istart2=ifirst2(i)
!      iend2=ilast2(i)           
if(ni==1)then
!   H(istart+offset1(istart))=uss(species(i))
     t1=uss(species(i))
     ! add a neighborlist here 

      do j=1,numat
         if(j==i)goto 1
          allocate(wout(pairs(species(i)),pairs(species(j))))
     call twoe_hcore_ints(i,j,wout)
     t1=t1-eff_core(zeff(j))*wout(1,1)
     deallocate(wout)
1        continue
      end do
!      H(istart+offset1(istart))=H(istart+offset1(istart))-t1
      if(dabs(t1).gt.thresh)then
         itotalh=itotalh+1
       if(itotalh.gt.hdim_max)call hcore_over
      H(itotalh)=t1
           hpair(1,itotalh)=istart
           hpair(2,itotalh)=istart
!      print*,istart,istart,H(itotalh)
      end if
!if(sparkles)then
!   allocate(wout(pairs(species(i)),1))
!   do j=1,nsparkle
!   call twoe_sparkle(i,j,wout)
!   H(istart+offset1(istart))=H(istart+offset1(istart))-spcharge(zeffsp(j))*wout(1,1)
!   end do
!   deallocate(wout)

!end if





   else
!      H(istart+offset1(istart))=uss(species(i))
allocate(htemp(lmax))
   htemp=zero
   htemp(1)=uss(species(i))
   htemp(5)=upp(species(i))
   htemp(8)=upp(species(i))
   htemp(10)=upp(species(i))
   if(ni>4)then
   htemp(15)=udd(species(i))
   htemp(21)=udd(species(i))
   htemp(28)=udd(species(i))
   htemp(36)=udd(species(i))
   htemp(45)=udd(species(i))
   end if

    
      do j=1,numat
         if(j==i)goto 2
        allocate(wout(pairs(species(i)),pairs(species(j))))
         charge=eff_core(zeff(j))
       call twoe_hcore_ints(i,j,wout)
           do k=1,lmax
         htemp(k)=htemp(k)-charge*wout(k,1)
          end do
          deallocate(wout)
2      continue
      end do


!if(sparkles)then
!   allocate(wout(pairs(species(i)),1))
!   do j=1,nsparkle
!   call twoe_sparkle(i,j,wout)

!   do k=0,lmax-1
! htemp(k+1)=htemp(k+1)+spcharge(zeffsp(j))*wout(k+1,1)
!          end do
!   end do
!deallocate(wout)
!end if


!H(istart+offset1(istart))=H(istart+offset1(istart))-htemp(1)
!H(offset1(istart+1)+istart)=-htemp(2)
!H(offset1(istart+2)+istart)=-htemp(3)
!H(offset1(istart+3)+istart)=-htemp(4)
!print*,istart,istart,htemp(1)
!print*,istart+1,istart,htemp(2)
!print*,istart+2,istart,htemp(3)
!print*,istart+3,istart,htemp(4)
if(dabs(htemp(1)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(1)
hpair(1,itotalh)=istart
hpair(2,itotalh)=istart
end if
if(dabs(htemp(2)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(2)
hpair(1,itotalh)=istart+1
hpair(2,itotalh)=istart
end if
if(dabs(htemp(3)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(3)
hpair(1,itotalh)=istart+2
hpair(2,itotalh)=istart
end if
if(dabs(htemp(4)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(4)
hpair(1,itotalh)=istart+3
hpair(2,itotalh)=istart
end if




if(ni>4)then
!H(offset1(istart+4)+istart)=-htemp(11)
!H(offset1(istart+5)+istart)=-htemp(22)
!H(offset1(istart+6)+istart)=-htemp(37)
!H(offset1(istart+7)+istart)=-htemp(16)
!H(offset1(istart+8)+istart)=-htemp(29)
!print*,istart+4,istart,htemp(11)
!print*,istart+5,istart,htemp(22)
!print*,istart+6,istart,htemp(37)
!print*,istart+7,istart,htemp(16)
!print*,istart+8,istart,htemp(29)
if(dabs(htemp(11)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(11)
hpair(1,itotalh)=istart+4
hpair(2,itotalh)=istart
end if
if(dabs(htemp(22)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(22)
hpair(1,itotalh)=istart+5
hpair(2,itotalh)=istart
end if
if(dabs(htemp(37)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(37)
hpair(1,itotalh)=istart+6
hpair(2,itotalh)=istart
end if
if(dabs(htemp(16)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(16)
hpair(1,itotalh)=istart+7
hpair(2,itotalh)=istart
end if
if(dabs(htemp(29)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(29)
hpair(1,itotalh)=istart+8
hpair(2,itotalh)=istart
end if

end if


i1=istart+1
!H(i1+offset1(i1))=H(i1+offset1(i1))-htemp(5)
!H(offset1(istart+2)+i1)=-htemp(6)
!H(offset1(istart+3)+i1)=-htemp(7)
!print*,i1,i1,htemp(5)
!print*,i1,istart+2,htemp(6)
!print*,i1,istart+3,htemp(7)
if(dabs(htemp(5)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(5)
hpair(1,itotalh)=i1
hpair(2,itotalh)=i1
end if
if(dabs(htemp(6)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(6)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+2
end if
if(dabs(htemp(7)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(7)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+3
end if


if(ni>4)then

!H(offset1(istart+4)+i1)=-htemp(12)
!H(offset1(istart+5)+i1)=-htemp(23)
!H(offset1(istart+6)+i1)=-htemp(38)
!H(offset1(istart+7)+i1)=-htemp(17)
!H(offset1(istart+8)+i1)=-htemp(30)
!print*,istart+4,i1,htemp(12)
!print*,istart+5,i1,htemp(23)
!print*,istart+6,i1,htemp(38)
!print*,istart+7,i1,htemp(17)
!print*,istart+8,i1,htemp(30)
if(dabs(htemp(12)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(12)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+4
end if
if(dabs(htemp(23)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(23)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+5
end if
if(dabs(htemp(38)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(38)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+6
end if
if(dabs(htemp(17)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(17)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+7
end if
if(dabs(htemp(30)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(30)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+8
end if
end if


i1=istart+2
!H(offset1(i1)+i1)=H(offset1(i1)+i1)-htemp(8)
!H(offset1(istart+3)+i1)=-htemp(9)
!print*,i1,i1,htemp(8)
!print*,istart+3,i1,htemp(9)
if(dabs(htemp(8)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(8)
hpair(1,itotalh)=i1
hpair(2,itotalh)=i1
end if
if(dabs(htemp(9)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(9)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+3
end if






if(ni>4)then
!H(offset1(istart+4)+i1)=-htemp(13)
!H(offset1(istart+5)+i1)=-htemp(24)
!H(offset1(istart+6)+i1)=-htemp(39)
!H(offset1(istart+7)+i1)=-htemp(18)
!H(offset1(istart+8)+i1)=-htemp(31)
!print*,istart+4,i1,htemp(13)
!print*,istart+5,i1,htemp(24)
!print*,istart+6,i1,htemp(39)
!print*,istart+7,i1,htemp(18)
!print*,istart+8,i1,htemp(31)
if(dabs(htemp(13)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(13)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+4
end if
if(dabs(htemp(24)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(24)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+5
end if
if(dabs(htemp(39)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(39)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+6
end if
if(dabs(htemp(18)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(18)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+7
end if
if(dabs(htemp(31)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(31)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+8
end if

end if




i1=istart+3
!H(i1+offset1(i1))=H(i1+offset1(i1))-htemp(10)
!print*,i1,i1,htemp(10)
if(dabs(htemp(10)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(10)
hpair(1,itotalh)=i1
hpair(2,itotalh)=i1
end if


if(ni>4)then
!H(offset1(istart+4)+i1)=-htemp(14)
!H(offset1(istart+5)+i1)=-htemp(25)
!H(offset1(istart+6)+i1)=-htemp(40)
!H(offset1(istart+7)+i1)=-htemp(19)
!H(offset1(istart+8)+i1)=-htemp(32)
!print*,istart+4,i1,htemp(14)
!print*,istart+5,i1,htemp(25)
!print*,istart+6,i1,htemp(40)
!print*,istart+7,i1,htemp(19)
!print*,istart+8,i1,htemp(32)
if(dabs(htemp(14)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(14)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+4
end if
if(dabs(htemp(25)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(25)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+5
end if
if(dabs(htemp(40)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(40)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+6
end if
if(dabs(htemp(19)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(19)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+7
end if
if(dabs(htemp(32)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(32)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+8
end if



! dump d-orbitals also
i1=istart+4
!H(i1+offset1(i1))=H(i1+offset1(i1))-htemp(15)
!H(offset1(istart+5)+i1)=-htemp(26)
!H(offset1(istart+6)+i1)=-htemp(41)
!H(offset1(istart+7)+i1)=-htemp(20)
!H(offset1(istart+8)+i1)=-htemp(33)
if(dabs(htemp(15)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(15)
hpair(1,itotalh)=i1
hpair(2,itotalh)=i1
end if
if(dabs(htemp(26)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(26)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+5
end if
if(dabs(htemp(41)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(41)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+6
end if
if(dabs(htemp(20)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(20)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+7
end if
if(dabs(htemp(33)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(33)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+8
end if




i1=istart+5
!H(i1+offset1(i1))=H(i1+offset1(i1))-htemp(28)
!H(offset1(istart+6)+i1)=-htemp(43)
!H(offset1(istart+7)+i1)=-htemp(27)
!H(offset1(istart+8)+i1)=-htemp(35)
if(dabs(htemp(28)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(28)
hpair(1,itotalh)=i1
hpair(2,itotalh)=i1
end if
if(dabs(htemp(43)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(43)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+6
end if
if(dabs(htemp(27)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(27)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+7
end if
if(dabs(htemp(35)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(35)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+8
end if



i1=istart+6

!H(i1+offset1(i1))=H(i1+offset1(i1))-htemp(45)
!H(offset1(istart+7)+i1)=-htemp(42)
!H(offset1(istart+8)+i1)=-htemp(44)
if(dabs(htemp(45)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(45)
hpair(1,itotalh)=i1
hpair(2,itotalh)=i1
end if
if(dabs(htemp(42)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(42)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+7
end if
if(dabs(htemp(44)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(44)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+8
end if




i1=istart+7

!H(i1+offset1(i1))=H(i1+offset1(i1))-htemp(21)
!H(offset1(istart+8)+i1)=-htemp(34)
if(dabs(htemp(21)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(21)
hpair(1,itotalh)=i1
hpair(2,itotalh)=i1
end if
if(dabs(htemp(34)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(34)
hpair(1,itotalh)=i1
hpair(2,itotalh)=istart+8
end if



i1=istart+8
!H(i1+offset1(i1))=H(i1+offset1(i1))-htemp(36)
if(dabs(htemp(36)).gt.thresh)then
itotalh=itotalh+1
if(itotalh.gt.hdim_max)call hcore_over
H(itotalh)=htemp(36)
hpair(1,itotalh)=i1
hpair(2,itotalh)=i1
end if

end if








deallocate(htemp)
end if

ilasthc1c(i)=itotalh
end do
nelectrons=nelectrons-sys_charge

!call mpi_barrier(mpi_comm_world,ierr)
!itotalints_all_cpuh=0
!call mpi_reduce(itotalh,itotalints_all_cpuh,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'There are ',itotalints_all_cpuh,'hcore integrals above',thresh
!call mpi_barrier(mpi_comm_world,ierr)
list_length=0
call mpi_gather(itotalh,1,mpi_integer,list_length,1,mpi_integer,0,mpi_comm_world,ierr)

if(myrank.eq.0)then
if(debug)then
print*,''
print*,''
itotalints_all_cpu=0
totalmem_all_cpu=0
print*,'Length and Memory for Hcore buffer on each core'
print*,'     Processor        Number Ints           Total Mem (MB)'
do i=1,nprocs
print*,i,'       ',list_length(i),'       ',list_length(i)*8/1d6
itotalints_all_cpu=itotalints_all_cpu+list_length(i)
totalmem_all_cpu=totalmem_all_cpu+float(list_length(i))*8/1d6
end do
print*,'---------------------------------------------------------------'
print*,'      Total:      ',itotalints_all_cpu,'       ',totalmem_all_cpu
end if
end if


call mpi_barrier(mpi_comm_world,ierr)

!do i=1,nsparkle
!nelectrons=nelectrons+spcharge(zeffsp(i))
!end do

!do i=1,itotalh
!print*,'hcoreint',hpair(1,i),hpair(2,i),h(i)
!end do

ipart=0
call mpi_reduce(nelectrons,ipart,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
nelectrons=ipart
call mpi_bcast(nelectrons,1,mpi_double_precision,0,mpi_comm_world,ierr)



end subroutine hcore


subroutine hcore_over
use constants
use tables
use indices
use control
implicit double precision (a-h,o-z)
include 'mpif.h'
print*,'Core Hamiltonian buffer overrun on process',myrank
stop
end subroutine hcore_over
