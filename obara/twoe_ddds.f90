subroutine twoe_ddds(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,pvector,qvector,wvector,dddsints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface

subroutine twoe_ddps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishellpc,pvector,qvector,wvector,ddpsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:),intent(inout)::ddpsints
end subroutine twoe_ddps

subroutine twoe_dpps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpc,pvector,qvector,wvector,dppsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellpb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dppsints
end subroutine twoe_dpps




end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishelldc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:),intent(inout)::dddsints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::ppssints,ppssintsmp1
double precision,dimension(3,3,3)::dpssints
double precision,dimension(3,3,3,3)::ddssints,dppsints
double precision,dimension(3,3,3,3,3)::ddpsints

do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishelldc)
wiipi(i)=wvector(i)-qvector(i)
end do

call twoe_ddps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,pvector,qvector,wvector,ddpsints)
ddpsintsm_global=ddpsints
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddsints(i,j,k,l,m,n)=qkkck(n)*ddpsints(i,j,k,l,m)
end do
end do
end do
end do
end do
end do

factor=1.0d0/(2.0d0*zc_plus_zd)
factor2=factor*(za_plus_zb)/zabcd
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddsints(i,j,k,l,m,m)=dddsints(i,j,k,l,m,m)+factor*ddssm_global(i,j,k,l)-factor2*ddssmp1_global(i,j,k,l)
end do
end do
end do
end do
end do


call twoe_ddps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,pvector,qvector,wvector,ddpsints)
ddpsintsmp1_global=ddpsints
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddsints(i,j,k,l,m,n)=dddsints(i,j,k,l,m,n)+wiipi(n)*ddpsints(i,j,k,l,m)
end do
end do
end do
end do
end do
end do



call twoe_dpps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,pvector,qvector,wvector,dppsints)


factor=1.0d0/(2.0d0*zabcd)
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddsints(i,j,k,l,m,k)=dddsints(i,j,k,l,m,k)+factor*dppsints(i,j,l,m)
end do
end do
end do
end do
end do

do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddsints(i,j,k,l,m,l)=dddsints(i,j,k,l,m,l)+factor*dppsints(i,j,k,m)
end do
end do
end do
end do
end do


! swap, sorta
call twoe_dpps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelldb,ishellda,ishelldc,pvector,qvector,wvector,dppsints)


factor=1.0d0/(2.0d0*zabcd)
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddsints(i,j,k,l,m,i)=dddsints(i,j,k,l,m,i)+factor*dppsints(k,l,j,m)
end do
end do
end do
end do
end do

do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddsints(i,j,k,l,m,j)=dddsints(i,j,k,l,m,j)+factor*dppsints(k,l,i,m)
end do
end do
end do
end do
end do








return
end


