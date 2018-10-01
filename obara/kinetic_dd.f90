subroutine kinetic_dd(distsq,dexpa,dexpb,ishella,ishellb,ddints,pvector,dkddints)
! (p|p) becomes (p|s) by VRR (P_j-B_j)
use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine overlap_dp(distsq,pexp,dexpb,ishellp,ishelld,dpints,pvector)
double precision,dimension(:,:,:),intent(inout)::dpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpb
integer,intent(in)::ishellp,ishelld
end subroutine overlap_dp

subroutine kinetic_dp(distsq,pexp,dexpa,ishellp,ishelld,dpints,pvector,dkdpints)
double precision,dimension(:,:,:),intent(inout)::dkdpints
double precision,dimension(:,:,:),intent(in)::dpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpa
integer,intent(in)::ishelld,ishellp
end subroutine kinetic_dp

subroutine kinetic_pp(distsq,pexpa,pexpb,ishellpa,ishellpb,ppints,pvector,dkppints)
double precision,dimension(:,:),intent(inout)::dkppints
double precision,dimension(:,:),intent(in)::ppints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexpa,pexpb
integer,intent(in)::ishellpa,ishellpb
end subroutine kinetic_pp


end interface


double precision,dimension(:,:,:,:),intent(inout)::dkddints
double precision,dimension(:,:,:,:),intent(in)::ddints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,dexpb,dexpa
integer,intent(in)::ishella,ishellb

double precision,dimension(3)::pjjbj,psints,dkpsints
double precision,dimension(3,3)::dkdsints,dsints,dkppints
double precision,dimension(3,3,3)::dpints,dkdpints
pi=dacos(-1.0d0)


do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellb)
end do


call overlap_dp(distsq,dexpb,dexpa,ishellb,ishella,dpints,pvector) ! generates ssglobal for the shell

call kinetic_dp(distsq,dexpb,dexpa,ishellb,ishella,dpints,pvector,dkdpints) ! generates dkss_global

do l=1,3
do k=1,3
do j=1,3
do i=1,3
dkddints(i,j,k,l)=pjjbj(l)*dkdpints(i,j,k)
end do
end do
end do
end do

r=2.0d0*(dexpb+dexpa)
do k=1,3
do j=1,3
do i=1,3
dkddints(i,j,k,k)=dkddints(i,j,k,k)+dkdsints_global(i,j)/r - zeta_global*dsints_global(i,j)/dexpb
end do
end do
end do

dkddints=dkddints+2.0d0*zeta_global*ddints ! zeta_global generated in kinetic_ss 


call kinetic_pp(distsq,dexpa,dexpb,ishella,ishellb,ppints_global,pvector,dkppints)
!ppints_global from the overlap_dd call in pick2.f90 


do k=1,3
do j=1,3
do i=1,3
dkddints(i,j,k,i)=dkddints(i,j,k,i)+dkppints(j,k)/r
end do
end do
end do

do k=1,3
do j=1,3
do i=1,3
dkddints(i,j,k,j)=dkddints(i,j,k,j)+dkppints(i,k)/r
end do
end do
end do





return
end
