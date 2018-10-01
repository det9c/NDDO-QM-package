subroutine kinetic_ds(distsq,sexp,dexp,ishelld,ishells,dsints,pvector,dkdsints)
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

subroutine kinetic_ps(distsq,sexp,pexp,ishellp,ishells,psints,pvector,dkspints)
double precision,dimension(:),intent(in)::psints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,pexp
integer,intent(in)::ishellp,ishells
double precision,intent(inout),dimension(:)::dkspints
end subroutine kinetic_ps


end interface


double precision,dimension(:,:),intent(inout)::dkdsints
double precision,dimension(:,:),intent(in)::dsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexp
integer,intent(in)::ishelld,ishells

double precision,dimension(3)::pjjbj,psints,dkpsints
pi=dacos(-1.0d0)


do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishelld)
end do


call overlap_ps(distsq,sexp,dexp,ishelld,ishells,psints,pvector) ! generates ssglobal for the shell

call kinetic_ps(distsq,sexp,dexp,ishelld,ishells,psints,pvector,dkpsints) ! generates dkss_global
dkpsints_global=dkpsints


do i=1,3
do j=1,3
dkdsints(j,i)=pjjbj(j)*dkpsints(i)
end do
end do



r=2.0d0*(sexp+dexp)
do i=1,3
dkdsints(i,i)=dkdsints(i,i)+dkss_global/r -zeta_global*ssglobal/dexp !ssglobal is computed in overlap_ps
end do

dkdsints=dkdsints+2.0d0*zeta_global*dsints ! zeta_global generated in kinetic_ss 


return
end
