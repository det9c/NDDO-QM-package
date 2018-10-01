subroutine twoe_psps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpc,pvector,qvector,wvector,pspsints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface
subroutine twoe_psss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellp,pvector,qvector,wvector,psssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:),intent(inout)::psssints
end subroutine twoe_psss
end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::pspsints
double precision,dimension(3)::qkkck,wkkqk,psssints

pi=dacos(-1.0d0)



do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellpc)
wkkqk(i)=wvector(i)-qvector(i)
end do

call twoe_psss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,pvector,qvector,wvector,psssints)

do k=1,3
do i=1,3
pspsints(i,k)=qkkck(k)*psssints(i)
end do
end do

do i=1,3
pspsints(i,i)=pspsints(i,i)+ssss_mp1_global/2.0d0/zabcd
end do

call twoe_psss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellpa,pvector,qvector,wvector,psssints)

do k=1,3
do i=1,3
pspsints(i,k)=pspsints(i,k)+wkkqk(k)*psssints(i)
end do
end do



return
end


