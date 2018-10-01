subroutine nuclear_dp(distsq,pexp,dexpb,ishellp,ishelld,pvector,mindex,dndpints)
use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine nuclear_ds(distsq,sexp,dexpb,ishelld,ishells,pvector,mindex,dndsints)
double precision,dimension(:,:),intent(inout)::dndsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexpb
integer,intent(in)::ishells,ishelld,mindex
end subroutine nuclear_ds


end interface


double precision,dimension(:,:,:),intent(inout)::dndpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpb
integer,intent(in)::ishellp,ishelld,mindex
double precision,dimension(3)::pjjbj,pjjcj,dnpsints
double precision,dimension(3,3)::dndsints

pi=dacos(-1.0d0)


do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellp)
pjjcj(i)=pvector(i)-coor(i,inuc)
end do



call nuclear_ds(distsq,pexp,dexpb,ishelld,ishellp,pvector,mindex,dndsints)
dsintsm_global=dndsints
do k=1,3
do j=1,3
do i=1,3
dndpints(i,j,k)=pjjbj(k)*dndsints(i,j)
end do
end do
end do


r=2.0d0*(pexp+dexpb)
do j=1,3
do i=1,3
dndpints(i,j,i)=dndpints(i,j,i)+(psintsm_global(j)-psintsmp1_global(j))/r 
end do
end do


do j=1,3
do i=1,3
dndpints(i,j,j)=dndpints(i,j,j)+(psintsm_global(i)-psintsmp1_global(i))/r
end do
end do

call nuclear_ds(distsq,pexp,dexpb,ishelld,ishellp,pvector,mindex+1,dndsints)
dsintsmp1_global=dndsints
do k=1,3
do j=1,3
do i=1,3
dndpints(i,j,k)=dndpints(i,j,k)-pjjcj(k)*dndsints(i,j)
end do
end do
end do





return
end
