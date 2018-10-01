subroutine overlap_dp(distsq,pexp,dexpb,ishellp,ishelld,dpints,pvector)

use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine overlap_ds(distsq,sexp,dexp,idshell,isshell,dsints,pvector)
double precision,dimension(:,:),intent(inout)::dsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexp
integer,intent(in)::isshell,idshell
end subroutine overlap_ds
end interface


double precision,dimension(:,:,:),intent(inout)::dpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpb
integer,intent(in)::ishellp,ishelld
double precision,dimension(3)::psints
double precision,dimension(3,3)::dsints
double precision,dimension(3)::pjjbj
pi=dacos(-1.0d0)

do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellp)
end do

call overlap_ds(distsq,pexp,dexpb,ishelld,ishellp,dsints,pvector)
dsints_global=dsints

do k=1,3
do j=1,3
do i=1,3
dpints(i,j,k)=pjjbj(k)*dsints(i,j)
end do
end do
end do

r=2.0d0*(pexp+dexpb)

do j=1,3
do i=1,3
dpints(i,j,i)=dpints(i,j,i)+psints_global(j)/r
end do
end do


do j=1,3
do i=1,3
dpints(i,j,j)=dpints(i,j,j)+psints_global(i)/r
end do
end do



return
end
