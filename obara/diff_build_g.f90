subroutine diff_build_g(ediff)
use gaussian_basis
implicit double precision (a-h,o-z)
include 'mpif.h'
interface
subroutine pack(indi,indj,ij)
integer,intent(in)::indi,indj
integer,intent(inout)::ij
end subroutine pack
end interface

double precision,intent(inout)::ediff
double precision,dimension(600)::rbuffer
integer,dimension(4,600)::inds
half=.5d0

!g=0.
ediff=0.
 do irow=1,itotal

imu=bfuncs_diff(1,irow)
inu=bfuncs_diff(2,irow)
ila=bfuncs_diff(3,irow)
isig=bfuncs_diff(4,irow)
value=diff_shell(irow)

!!!if(io<0)exit
!coulomb contribution
double=1.0d0
scale1=1.
scale2=1.
call pack(imu,inu,imu_nu)
call pack(ila,isig,ila_sig)
if(imu.ne.inu)scale1=2.0d0
if(ila.ne.isig)scale2=2.0d0
if(imu_nu .eq. ila_sig)double=.5

wp1=2.
wp2=2.
if(imu.eq.inu)wp1=1.0
if(ila.eq.isig)wp2=1.0
!g(imu_nu)=g(imu_nu)+density_global(ila_sig)*scale2*value*double
!g(ila_sig)=g(ila_sig)+density_global(imu_nu)*scale1*value*double
ediff=ediff+density_global(ila,isig)*scale2*value*double*wp1*density_global(imu,inu)
ediff=ediff+density_global(imu,inu)*scale1*value*double*wp2*density_global(ila,isig)

wp1=2.
wp2=2.
wp3=2.
wp4=2.
if(imu.eq.ila)wp1=1.0
if(imu.eq.isig)wp2=1.0
if(inu.eq.ila)wp3=1.
if(inu.eq.isig)wp4=1.


!exchange contribution
call pack(imu,ila,imu_la)
call pack(imu,isig,imu_sig)
call pack(inu,ila,inu_la)
call pack(inu,isig,inu_sig)
factor=1.0
if(ila .eq. isig)factor=half
if(imu .eq. inu )factor=half*factor
scale3=1.
if(imu .eq. ila .and. inu.ne.isig)scale3=2.
!g(imu_la)=g(imu_la)- half*density_global(inu_sig)*value*scale3*factor
ediff=ediff- half*density_global(inu,isig)*value*scale3*factor*wp1*density_global(imu,ila)
scale3=1.
factor2=1.0d0
if(imu_nu .eq. ila_sig .and. imu .ne. inu)factor2=.5d0 !(uv|uv) ->uu uv vu vv so need half since uv vu is twice
if(imu .eq. isig .and.  inu.ne.ila)scale3=2.
!g(imu_sig)=g(imu_sig)- half*density_global(inu_la)*value*scale3*factor*factor2
ediff=ediff - half*density_global(inu,ila)*value*scale3*factor*factor2*wp2*density_global(imu,isig)
scale3=1.
if(inu .eq. ila .and.  imu.ne.isig)scale3=2.
!g(inu_la)=g(inu_la)- half*density_global(imu_sig)*value*scale3*factor*factor2
ediff=ediff- half*density_global(imu,isig)*value*scale3*factor*factor2*wp3*density_global(inu,ila)
scale3=1.
if(inu.eq.isig  .and.  imu.ne.ila)scale3=2.
!g(inu_sig)=g(inu_sig)- half*density_global(imu_la)*value*scale3*factor
ediff=ediff- half*density_global(imu,ila)*value*scale3*factor*wp4*density_global(inu,isig)


end do

!call mpi_reduce(g,glocal,iveclength,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!g=glocal

!call mpi_barrier(mpi_comm_world,ierr)

ediff=ediff*0.5
return
end
