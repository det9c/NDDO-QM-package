subroutine nuclear_ds(distsq,sexp,dexpb,ishelld,ishells,pvector,mindex,dndsints)
use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine nuclear_ps(distsq,sexp,pexp,ishellp,ishells,pvector,mindex,dnpsints)
double precision,dimension(:),intent(inout)::dnpsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,pexp
integer,intent(in)::ishellp,ishells,mindex
end subroutine nuclear_ps




end interface


double precision,dimension(:,:),intent(inout)::dndsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexpb
integer,intent(in)::ishells,ishelld,mindex
double precision,dimension(3)::pjjbj,pjjcj,dnpsints

pi=dacos(-1.0d0)


do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishelld)
pjjcj(i)=pvector(i)-coor(i,inuc)
end do



call nuclear_ps(distsq,sexp,dexpb,ishelld,ishells,pvector,mindex,dnpsints)
psintsm_global=dnpsints
do j=1,3
do i=1,3
dndsints(i,j)=pjjbj(j)*dnpsints(i)
end do
end do


r=2.0d0*(sexp+dexpb)
do i=1,3
dndsints(i,i)=dndsints(i,i)+(ssm_global-ssmp1_global)/r 
end do

call nuclear_ps(distsq,sexp,dexpb,ishelld,ishells,pvector,mindex+1,dnpsints)
psintsmp1_global=dnpsints
do j=1,3
do i=1,3
dndsints(i,j)=dndsints(i,j)-pjjcj(j)*dnpsints(i)
end do
end do





return
end
