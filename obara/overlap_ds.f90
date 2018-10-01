subroutine overlap_ds(distsq,sexp,dexp,idshell,isshell,dsints,pvector)
! (d|s) becomes (p|s) by VRR (P_j-B_j)
use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine overlap_ps(distsq,sexp,pexp,ishellp,ishells,psints,pvector)
double precision,dimension(:),intent(inout)::psints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,pexp
integer,intent(in)::ishellp,ishells
end subroutine overlap_ps
end interface


double precision,dimension(:,:),intent(inout)::dsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexp
integer,intent(in)::isshell,idshell
double precision,dimension(3)::psints

double precision,dimension(3)::pjjbj
pi=dacos(-1.0d0)

do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,idshell)
end do

call overlap_ps(distsq,sexp,dexp,idshell,isshell,psints,pvector)
psints_global=psints

icount=0
do j=1,3
do i=1,3
dsints(i,j)=pjjbj(j)*psints(i)
end do
end do



r=2.0d0*(sexp+dexp)
do i=1,3
dsints(i,i)=dsints(i,i)+ssglobal/r
end do



return
end
