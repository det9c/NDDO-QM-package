subroutine twoe_ddss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,pvector,qvector,wvector,ddssints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface

subroutine twoe_dpss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellp,pvector,qvector,wvector,dpssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:),intent(inout)::dpssints
end subroutine twoe_dpss

subroutine twoe_ppss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,pvector,qvector,wvector,ppssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::ppssints
end subroutine twoe_ppss


end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::ddssints
double precision,dimension(3)::piiai,wiipi,psssints
double precision,dimension(3,3)::ppssints,ppssintsmp1
double precision,dimension(3,3,3)::dpssints

pi=dacos(-1.0d0)

do i=1,3
piiai(i)=pvector(i)-shell_coor(i,ishelldb)
wiipi(i)=wvector(i)-pvector(i)
end do

call twoe_dpss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,pvector,qvector,wvector,dpssints)


do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddssints(i,j,k,l)=piiai(l)*dpssints(i,j,k)
end do
end do
end do
end do



factor=1.0d0/(2.0d0*za_plus_zb)
factor2=factor*(zc_plus_zd)/zabcd
do k=1,3
do j=1,3
do i=1,3
ddssints(i,j,k,k)=ddssints(i,j,k,k)+factor*dsssm_global(i,j)-factor2*dsssmp1_global(i,j)
end do
end do
end do


call twoe_ppss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,pvector,qvector,wvector,ppssints)
call twoe_ppss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,pvector,qvector,wvector,ppssintsmp1)

do k=1,3
do j=1,3
do i=1,3
ddssints(i,j,k,i)=ddssints(i,j,k,i)+factor*ppssints(j,k)-factor2*ppssintsmp1(j,k)
end do
end do
end do


do k=1,3
do j=1,3
do i=1,3
ddssints(i,j,k,j)=ddssints(i,j,k,j)+factor*ppssints(i,k)-factor2*ppssintsmp1(i,k)
end do
end do
end do




call twoe_dpss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,pvector,qvector,wvector,dpssints)
dpssmp1_global=dpssints


do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddssints(i,j,k,l)=ddssints(i,j,k,l)+wiipi(l)*dpssints(i,j,k)
end do
end do
end do
end do



return
end


