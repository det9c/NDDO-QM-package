subroutine kinetic_ps(distsq,sexp,pexp,ishellp,ishells,psints,pvector,dkpsints)
use gaussian_basis
implicit double precision (a-h,o-z)

double precision,dimension(:),intent(in)::psints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,pexp
integer,intent(in)::ishellp,ishells
double precision,intent(inout),dimension(:)::dkpsints
double precision,dimension(3)::piiai

pi=dacos(-1.0d0)
do i=1,3
piiai(i)=pvector(i)-shell_coor(i,ishellp)

end do

call kinetic_ss(distsq,pexp,sexp,ssglobal,dkssout)
dkss_global=dkssout


dkpsints=piiai*dkssout + 2.0d0*zeta_global*psints



return
end
