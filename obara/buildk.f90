subroutine buildk(distsq,expa,expb,icola,icolb,dkab,za_plus_zb)
use gaussian_basis
implicit double precision (a-h,o-z)
double precision,intent(in)::distsq,expa,expb
integer,intent(in)::icola,icolb
double precision,intent(inout)::dkab,za_plus_zb


za_plus_zb=expa+expb
dkab=dkabfactor*exp(-expa*expb*distsq/za_plus_zb)
dkab=dkab/za_plus_zb

return
end
