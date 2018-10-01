subroutine nuclear_pp(distsq,pexpa,pexpb,ishellpa,ishellpb,pvector,mindex,dnppints)
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


double precision,dimension(:,:),intent(inout)::dnppints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexpa,pexpb
integer,intent(in)::ishellpa,ishellpb,mindex
double precision,dimension(3)::pjjbj,pjjcj,dnpsints

pi=dacos(-1.0d0)


do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellpb)
pjjcj(i)=pvector(i)-coor(i,inuc)
end do



call nuclear_ps(distsq,pexpb,pexpa,ishellpa,ishellpb,pvector,mindex,dnpsints)
do j=1,3
do i=1,3
dnppints(i,j)=pjjbj(j)*dnpsints(i)
end do
end do


r=2.0d0*(pexpa+pexpb)
do i=1,3
dnppints(i,i)=dnppints(i,i)+(ssm_global-ssmp1_global)/r 
end do

call nuclear_ps(distsq,pexpb,pexpa,ishellpa,ishellpb,pvector,mindex+1,dnpsints)
do j=1,3
do i=1,3
dnppints(i,j)=dnppints(i,j)-pjjcj(j)*dnpsints(i)
end do
end do




!do i=1,3
!dnppints(i,i)=dnppints(i,i)-ssmp1_global/r 
!end do




return
end
