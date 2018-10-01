subroutine nuclear_ss(distsq,r,t,ss,pvector,mindex,snsout)
use gaussian_basis
implicit double precision (a-h,o-z)

interface
subroutine fzero(uz,series,low)
double precision,intent(in)::uz
double precision,intent(inout)::series
integer,intent(in)::low
end subroutine fzero
end interface


double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,r,t,ss
double precision,dimension(3)::piici
double precision,intent(inout)::snsout
integer,intent(in)::mindex

pi=dacos(-1.0d0)

rplust=r+t

snsout=0.
!do inuc=1,natoms

bigU=0.
do i=1,3
piici(i)=pvector(i)-coor(i,inuc)
bigU=bigU+piici(i)*piici(i)
end do
bigU=bigU*rplust





if(bigU.gt.39.99D0.and.mindex.ne.0)then
 max=2*mindex-1
 totalt=1.0D0
 do 1155 l=max,1,-2
        totalt=totalt*l
 1155  continue
 fzzero=(totalt/(2.0D0*(2.0*bigU)**mindex))*dsqrt((pi/bigU))
elseif(bigU.gt.39.99D0.and.mindex.eq.0)then
 fzzero=0.5D0*dsqrt((pi/bigU))
else 
 call fzero(bigU,fzzero,mindex)
end if



      snsout=-zcore(inuc)*2.0D0*SQRT(rplust/pi)*ss*fzzero
!end do


return
end
