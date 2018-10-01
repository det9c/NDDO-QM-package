subroutine nuclear_dd(distsq,dexpa,dexpb,ishella,ishellb,pvector,mindex,dnddints)
use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine nuclear_pp(distsq,pexpa,pexpb,ishellpa,ishellpb,pvector,mindex,dnppints)
double precision,dimension(:,:),intent(inout)::dnppints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexpa,pexpb
integer,intent(in)::ishellpa,ishellpb,mindex
end subroutine nuclear_pp

subroutine nuclear_dp(distsq,pexp,dexpb,ishellp,ishelld,pvector,mindex,dndpints)
double precision,dimension(:,:,:),intent(inout)::dndpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpb
integer,intent(in)::ishellp,ishelld,mindex
end subroutine nuclear_dp




end interface


double precision,dimension(:,:,:,:),intent(inout)::dnddints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,dexpa,dexpb
integer,intent(in)::ishella,ishellb,mindex
double precision,dimension(3)::pjjbj,pjjcj,dnpsints
double precision,dimension(3,3,3)::dndpints
double precision,dimension(3,3)::dnppints,dnppints1

  
pi=dacos(-1.0d0)


do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellb)
pjjcj(i)=pvector(i)-coor(i,inuc)
end do


call nuclear_dp(distsq,dexpb,dexpa,ishellb,ishella,pvector,mindex,dndpints)




do l=1,3
do k=1,3
do j=1,3
do i=1,3
dnddints(i,j,k,l)=pjjbj(l)*dndpints(i,j,k)
end do
end do
end do
end do



r=2.0d0*(dexpa+dexpb)
do k=1,3
do j=1,3
do i=1,3
dnddints(i,j,k,k)=dnddints(i,j,k,k)+(dsintsm_global(i,j)-dsintsmp1_global(i,j))/r 
end do
end do
end do


call nuclear_dp(distsq,dexpb,dexpa,ishellb,ishella,pvector,mindex+1,dndpints)

do l=1,3
do k=1,3
do j=1,3
do i=1,3
dnddints(i,j,k,l)=dnddints(i,j,k,l)-pjjcj(l)*dndpints(i,j,k)
end do
end do
end do
end do





call nuclear_pp(distsq,dexpa,dexpb,ishella,ishellb,pvector,mindex,dnppints)
call nuclear_pp(distsq,dexpa,dexpb,ishella,ishellb,pvector,mindex+1,dnppints1)

do k=1,3
do j=1,3
do i=1,3
dnddints(i,j,k,i)=dnddints(i,j,k,i)+(dnppints(j,k)-dnppints1(j,k))/r
end do
end do
end do

do k=1,3
do j=1,3
do i=1,3
dnddints(i,j,k,j)=dnddints(i,j,k,j)+(dnppints(i,k)-dnppints1(i,k))/r
end do
end do
end do






return
end
