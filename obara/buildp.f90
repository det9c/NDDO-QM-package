subroutine buildp(expa,expb,icola,icolb,pvector)
use gaussian_basis
double precision,intent(in)::expa,expb
integer,intent(in)::icola,icolb
double precision,dimension(3),intent(inout)::pvector

do i=1,3
pvector(i)=expa*shell_coor(i,icola)+expb*shell_coor(i,icolb)
end do
pvector=pvector/(expa+expb)
return
end
