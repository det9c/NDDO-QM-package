subroutine overlap_dd(distsq,dexpa,dexpb,ishella,ishellb,ddints,pvector)

use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine overlap_pp(distsq,pexpa,pexpb,ishellpa,ishellpb,ppints,pvector)
double precision,dimension(:,:),intent(inout)::ppints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexpa,pexpb
integer,intent(in)::ishellpa,ishellpb
end subroutine overlap_pp

subroutine overlap_dp(distsq,pexp,dexpb,ishellp,ishelld,dpints,pvector)
double precision,dimension(:,:,:),intent(inout)::dpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpb
integer,intent(in)::ishellp,ishelld
end subroutine overlap_dp

end interface


double precision,dimension(:,:,:,:),intent(inout)::ddints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,dexpa,dexpb
integer,intent(in)::ishella,ishellb
double precision,dimension(3,3,3)::dpints
double precision,dimension(3,3)::dsints,ppints
double precision,dimension(3)::pjjbj
pi=dacos(-1.0d0)

do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellb)
end do

call overlap_dp(distsq,dexpb,dexpa,ishellb,ishella,dpints,pvector)

do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddints(i,j,k,l)=pjjbj(l)*dpints(i,j,k)
end do
end do
end do
end do

r=2.0d0*(dexpa+dexpb)



do k=1,3
do j=1,3
do i=1,3
ddints(i,j,k,k)=ddints(i,j,k,k)+dsints_global(i,j)/r
end do
end do
end do

call overlap_pp(distsq,dexpa,dexpb,ishella,ishellb,ppints,pvector)
ppints_global=ppints


do k=1,3
do j=1,3
do i=1,3
ddints(i,j,k,i)=ddints(i,j,k,i)+ppints(j,k)/r
end do
end do
end do

do k=1,3
do j=1,3
do i=1,3
ddints(i,j,k,j)=ddints(i,j,k,j)+ppints(i,k)/r
end do
end do
end do



return
end
