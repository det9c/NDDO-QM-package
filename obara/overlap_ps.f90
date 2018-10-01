subroutine overlap_ps(distsq,sexp,pexp,ishellp,ishells,psints,pvector)
! (p|s) order so need (P-A)
use gaussian_basis
implicit double precision (a-h,o-z)
double precision,dimension(:),intent(inout)::psints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,pexp
integer,intent(in)::ishellp,ishells

double precision,dimension(3)::piiai
pi=dacos(-1.0d0)

do i=1,3
piiai(i)=pvector(i)-shell_coor(i,ishellp)
end do
call overlap_ss(distsq,pexp,sexp,ssintegral)
ssglobal=ssintegral

psints=piiai*ssintegral



return
end
