subroutine twoe_ssss(mindex,zabcd,za_plus_zb,zc_plus_zd,pvector,qvector,ssss)
use gaussian_basis
implicit double precision (a-h,o-z)
interface
subroutine fzero(uz,series,low)
double precision,intent(in)::uz
double precision,intent(inout)::series
integer,intent(in)::low
end subroutine fzero

end interface


integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
double precision,dimension(:),intent(in)::pvector,qvector
double precision,intent(inout)::ssss

pi=dacos(-1.0d0)
bigT=0.
do i=1,3
bigT=bigT+(pvector(i)-qvector(i))**2.0
end do
bigT=bigT*za_plus_zb*zc_plus_zd/zabcd

if(bigT.gt.39.99D0.and.mindex.ne.0)then
 max=2*mindex-1
 totalt=1.0D0
 do 1155 l=max,1,-2
        totalt=totalt*l
 1155  continue
 fzzero=(totalt/(2.0D0*(2.0*bigT)**mindex))*dsqrt((pi/bigT))
elseif(bigT.gt.39.99D0.and.mindex.eq.0)then
 fzzero=0.5D0*dsqrt((pi/bigT))
else
 call fzero(bigT,fzzero,mindex)
end if



ssss=   (1.0d0/dsqrt(zabcd))*kabkcd*fzzero

return
end


