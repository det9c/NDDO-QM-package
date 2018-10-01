subroutine overlap_pp(distsq,pexpa,pexpb,ishellpa,ishellpb,ppints,pvector)
! (p|p) becomes (p|s) by VRR (P_j-B_j)
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


double precision,dimension(:,:),intent(inout)::ppints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexpa,pexpb
integer,intent(in)::ishellpa,ishellpb
double precision,dimension(3)::psints

double precision,dimension(3)::pjjbj
pi=dacos(-1.0d0)

do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellpb)
end do

call overlap_ps(distsq,pexpb,pexpa,ishellpa,ishellpb,psints,pvector)

do j=1,3
do i=1,3
ppints(i,j)=pjjbj(j)*psints(i)
end do
end do



r=2.0d0*(pexpa+pexpb)
do i=1,3
ppints(i,i)=ppints(i,i)+ssglobal/r !ssglobal is computed in overlap_ps
end do


return
end
