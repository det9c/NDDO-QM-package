program CMOL
use control
implicit double precision (a-h,o-z)

character(50) process
      include 'mpif.h'
     call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myrank,ierr)
      call mpi_comm_size(mpi_comm_world,nprocs,ierr)

      call mpi_get_processor_name(process,length,ierr)
if(myrank.eq.0)print*,'Total number of CPUs',nprocs
!if(myrank.eq.0)print*,'Processor names and ranks follow'
!call mpi_barrier(mpi_comm_world,ierr)
!      print*,process,myrank
!call mpi_barrier(mpi_comm_world,ierr)

master=.false.
if(myrank.eq.0)master=.true.
if(master)then
print*,'-----------------------------------------------'
print*,'|          Transfer Hamiltonian               |'
print*,'|       WRITTEN BY DeCARLOS TAYLOR (2003)     |'
print*,'|   Revision 2003.A (Current thru  12/31/03)  |'
print*,'|   Revision 2005.A (Current thru  10/11/05)  |'
print*,'-----------------------------------------------'
print*,''
end if
call cpusec(time1)
ty=mpi_wtime()
call readin
call cpusec(time2)
if(master)then
write(*,*)''
write(*,*)''
write(*,200)time2-time1
ty2=mpi_wtime()
print*,'timer in driver',ty,ty2
200 format(' Total job time (sec) =  ',F10.3)

write(*,*)''
write(*,*)'Calculation completed. Program will stop'
end if
call mpi_barrier(mpi_comm_world,ierr)
call mpi_finalize(ierr)
end program cmol

