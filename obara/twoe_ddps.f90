subroutine twoe_ddps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishellpc,pvector,qvector,wvector,ddpsints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface

subroutine twoe_ddss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,pvector,qvector,wvector,ddssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::ddssints
end subroutine twoe_ddss

subroutine twoe_dpss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellp,pvector,qvector,wvector,dpssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:),intent(inout)::dpssints
end subroutine twoe_dpss


end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:),intent(inout)::ddpsints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::ppssints,ppssintsmp1
double precision,dimension(3,3,3)::dpssints
double precision,dimension(3,3,3,3)::ddssints

pi=dacos(-1.0d0)

do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellpc)
wiipi(i)=wvector(i)-qvector(i)
end do

call twoe_ddss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,pvector,qvector,wvector,ddssints)
ddssm_global=ddssints
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddpsints(i,j,k,l,m)=qkkck(m)*ddssints(i,j,k,l)
end do
end do
end do
end do
end do

factor=1.0d0/(2.0d0*zabcd)
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddpsints(i,j,k,l,k)=ddpsints(i,j,k,l,k)+factor*dpssmp1_global(i,j,l)
end do
end do
end do
end do


do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddpsints(i,j,k,l,l)=ddpsints(i,j,k,l,l)+factor*dpssmp1_global(i,j,k)
end do
end do
end do
end do


call twoe_ddss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,pvector,qvector,wvector,ddssints)
ddssmp1_global=ddssints
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddpsints(i,j,k,l,m)=ddpsints(i,j,k,l,m)+wiipi(m)*ddssints(i,j,k,l)
end do
end do
end do
end do
end do


call twoe_dpss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelldb,ishellda,pvector,qvector,wvector,dpssints)

!!! had to use klj in loop below because formula is pdss and subroutine is dpss
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddpsints(i,j,k,l,i)=ddpsints(i,j,k,l,i)+factor*dpssints(k,l,j)
end do
end do
end do
end do

do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddpsints(i,j,k,l,j)=ddpsints(i,j,k,l,j)+factor*dpssints(k,l,i)
end do
end do
end do
end do





return
end


