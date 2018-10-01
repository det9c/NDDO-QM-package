subroutine twoe_psss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellp,pvector,qvector,wvector,psssints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface
subroutine twoe_ssss(mindex,zabcd,za_plus_zb,zc_plus_zd,pvector,qvector,ssss)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
double precision,dimension(:),intent(in)::pvector,qvector
double precision,intent(inout)::ssss
end subroutine twoe_ssss
end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:),intent(inout)::psssints
double precision,dimension(3)::piiai,wiipi

pi=dacos(-1.0d0)

do i=1,3
piiai(i)=pvector(i)-shell_coor(i,ishellp)
wiipi(i)=wvector(i)-pvector(i)
end do

call twoe_ssss(mindex,zabcd,za_plus_zb,zc_plus_zd,pvector,qvector,ssss_m)
ssss_m_global=ssss_m
call twoe_ssss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,pvector,qvector,ssss_mp1)
ssss_mp1_global=ssss_mp1
psssints=piiai*ssss_m + wiipi*ssss_mp1


return
end


