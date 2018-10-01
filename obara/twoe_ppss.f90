subroutine twoe_ppss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,pvector,qvector,wvector,ppssints)
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
integer,intent(in)::ishellpa,ishellpb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::ppssints
double precision,dimension(3)::pjjbj,wjjpj,psssints

pi=dacos(-1.0d0)



do i=1,3
pjjbj(i)=pvector(i)-shell_coor(i,ishellpb)
wjjpj(i)=wvector(i)-pvector(i)
end do

call twoe_psss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,pvector,qvector,wvector,psssints)

do j=1,3
do i=1,3
ppssints(i,j)=pjjbj(j)*psssints(i)
end do
end do

factor=1.0/(2.0d0*za_plus_zb)
factor2=factor*zc_plus_zd/zabcd
do i=1,3
ppssints(i,i)=ppssints(i,i)+ factor*ssss_m_global-factor2* ssss_mp1_global
end do

call twoe_psss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellpa,pvector,qvector,wvector,psssints)
psssints_mp1_global=psssints
do j=1,3
do i=1,3
ppssints(i,j)=ppssints(i,j)+wjjpj(j)*psssints(i)
end do
end do



return
end


