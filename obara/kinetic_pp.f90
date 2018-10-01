subroutine kinetic_pp(distsq,pexpa,pexpb,ishellpa,ishellpb,ppints,pvector,dkppints)
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


double precision,dimension(:,:),intent(inout)::dkppints
double precision,dimension(:,:),intent(in)::ppints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexpa,pexpb
integer,intent(in)::ishellpa,ishellpb

double precision,dimension(3)::pjjbj,psints,dkpsints
pi=dacos(-1.0d0)


do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellpb)
end do


call overlap_ps(distsq,pexpb,pexpa,ishellpa,ishellpb,psints,pvector) ! generates ssglobal for the shell

call kinetic_ps(distsq,pexpb,pexpa,ishellpa,ishellpb,psints,pvector,dkpsints) ! generates dkss_global



do j=1,3
do i=1,3
dkppints(i,j)=pjjbj(j)*dkpsints(i)
end do
end do



r=2.0d0*(pexpa+pexpb)
do i=1,3
dkppints(i,i)=dkppints(i,i)+dkss_global/r !ssglobal is computed in overlap_ps
end do

dkppints=dkppints+2.0d0*zeta_global*ppints ! zeta_global generated in kinetic_ss 


return
end
