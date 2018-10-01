subroutine repel(nuclear,r,iatom,jatom,ssss)
use constants
use tables
use control
implicit double precision (a-h,o-z)
!**************
interface
subroutine repam1(iatom,jatom,term,r)
double precision,intent(inout)::term
double precision,intent(in)::r
integer,intent(in)::iatom,jatom
end subroutine repam1
end interface
!**************


double precision,intent(inout)::nuclear,r
double precision,intent(in)::ssss
integer,intent(in)::iatom,jatom
integer,dimension(9)::isite


!nuclear=eff_core(zeff(iatom))*eff_core(zeff(jatom))*27.21d0/(r/.529167d0)
!return




nuclear=zero

if(  eff_core(zeff(iatom)).eq.6  .and. eff_core(zeff(jatom)).eq.1)then
nuclear=one+(r*dexp(-alpha(species(iatom))*r))+dexp(-alpha(species(jatom))*r)
elseif(eff_core(zeff(jatom)).eq.6 .and. eff_core(zeff(iatom)).eq.1)then
nuclear=one+(r*dexp(-alpha(species(jatom))*r))+dexp(-alpha(species(iatom))*r)
elseif(  zstore(iatom).eq.7  .and. zstore(jatom).eq.1)then
nuclear=one+(r*dexp(-alpha(species(iatom))*r))+dexp(-alpha(species(jatom))*r)
elseif(zstore(jatom).eq.7 .and. zstore(iatom).eq.1)then
nuclear=one+(r*dexp(-alpha(species(jatom))*r))+dexp(-alpha(species(iatom))*r)

!elseif(zstore(jatom).eq.79 .and. eff_core(zeff(iatom)).eq.1)then
!nuclear=one+(r*dexp(-alpha(species(jatom))*r))+dexp(-alpha(species(iatom))*r)
!elseif(zstore(iatom).eq.79 .and. eff_core(zeff(jatom)).eq.1)then
!nuclear=one+(r*dexp(-alpha(species(iatom))*r))+dexp(-alpha(species(jatom))*r)
!elseif(zstore(jatom).eq.80 .and. eff_core(zeff(iatom)).eq.1)then
!nuclear=one+(r*dexp(-alpha(species(jatom))*r))+dexp(-alpha(species(iatom))*r)
!elseif(zstore(iatom).eq.80 .and. eff_core(zeff(jatom)).eq.1)then
!nuclear=one+(r*dexp(-alpha(species(iatom))*r))+dexp(-alpha(species(jatom))*r)





else
nuclear=one+dexp(-alpha(species(iatom))*r)+dexp(-alpha(species(jatom))*r)
end if


nuclear=eff_core(zeff(iatom))*eff_core(zeff(jatom))*ssss*nuclear
!return
! if(isite(iatom).eq.isite(jatom))then

   if(method=='AM1'  .or.  core=='GAUSSIAN'  .or. method=='PM3')then
      call repam1(iatom,jatom,term,r)
      nuclear=nuclear+term
   end if
return
!else



iz=zeff(iatom)
jz=zeff(jatom)
izjz=iz+jz
rau=r/0.52918d0
if(izjz==16)then
pa=aoop
pb=boop
pc=coop
pd=doop
pe=eoop

elseif(izjz==9)then
pa=aohp
pb=bohp
pc=cohp
pd=dohp
pe=eohp

elseif(izjz==2)then
pa=ahhp
pb=bhhp
pc=chhp
pd=dhhp
pe=ehhp
end if




term=pa*exp(-pb * rau) + pc/rau**6 + pd / rau**8 + pe / rau**10
nuclear=nuclear+term*27.21d0
!end if








end subroutine repel

