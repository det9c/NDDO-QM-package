subroutine kinetic_dp(distsq,pexp,dexpa,ishellp,ishelld,dpints,pvector,dkdpints)
! (p|p) becomes (p|s) by VRR (P_j-B_j)
use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine overlap_ds(distsq,sexp,dexpa,idshell,isshell,dsints,pvector)
double precision,dimension(:,:),intent(inout)::dsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexpa
integer,intent(in)::isshell,idshell
end subroutine overlap_ds

subroutine kinetic_ds(distsq,sexp,dexpa,ishelld,ishells,dsints,pvector,dkdsints)
double precision,dimension(:,:),intent(inout)::dkdsints
double precision,dimension(:,:),intent(in)::dsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexpa
integer,intent(in)::ishelld,ishells
end subroutine kinetic_ds


end interface


double precision,dimension(:,:,:),intent(inout)::dkdpints
double precision,dimension(:,:,:),intent(in)::dpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpa
integer,intent(in)::ishelld,ishellp

double precision,dimension(3)::pjjbj,psints,dkpsints
double precision,dimension(3,3)::dkdsints,dsints
pi=dacos(-1.0d0)


do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellp)
end do


call overlap_ds(distsq,pexp,dexpa,ishelld,ishellp,dsints,pvector) ! generates ssglobal for the shell
dsints_global=dsints


call kinetic_ds(distsq,pexp,dexpa,ishelld,ishellp,dsints,pvector,dkdsints) ! generates dkss_global
dkdsints_global=dkdsints

do k=1,3
do j=1,3
do i=1,3
dkdpints(i,j,k)=pjjbj(k)*dkdsints(i,j)
end do
end do
end do


r=2.0d0*(pexp+dexpa)
do j=1,3
do i=1,3
dkdpints(i,j,i)=dkdpints(i,j,i)+dkpsints_global(j)/r 
end do
end do

do j=1,3
do i=1,3
dkdpints(i,j,j)=dkdpints(i,j,j)+dkpsints_global(i)/r 
end do
end do



dkdpints=dkdpints+2.0d0*zeta_global*dpints ! zeta_global generated in kinetic_ss 


return
end
