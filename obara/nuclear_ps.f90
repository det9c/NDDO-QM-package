subroutine nuclear_ps(distsq,sexp,pexp,ishellp,ishells,pvector,mindex,dnpsints)
use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine nuclear_ss(distsq,r,t,ss,pvector,mindex,snsout)
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,r,t,ss
double precision,intent(inout)::snsout
integer,intent(in)::mindex
end subroutine nuclear_ss
end interface


double precision,dimension(:),intent(inout)::dnpsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,pexp
integer,intent(in)::ishellp,ishells,mindex
double precision,dimension(3)::piiai,piici

pi=dacos(-1.0d0)


do i=1,3
piiai(i)=pvector(i)-shell_coor(i,ishellp)
piici(i)=pvector(i)-coor(i,inuc)
end do

call nuclear_ss(distsq,sexp,pexp,ssglobal,pvector,mindex,ss0)
ssm_global=ss0
call nuclear_ss(distsq,sexp,pexp,ssglobal,pvector,mindex+1,ss1)
ssmp1_global=ss1

dnpsints=piiai*ss0 - piici*ss1




return
end
