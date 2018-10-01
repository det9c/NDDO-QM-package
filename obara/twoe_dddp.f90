subroutine twoe_dddp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishellpd,pvector,qvector,wvector,dddpints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface

subroutine twoe_ddds(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,pvector,qvector,wvector,dddsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishelldc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:),intent(inout)::dddsints
end subroutine twoe_ddds

subroutine twoe_dsdp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,ishellpd,pvector,qvector,wvector,dsdpints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:),intent(inout)::dsdpints
end subroutine twoe_dsdp


end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:,:),intent(inout)::dddpints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::ppssints,ppssintsmp1
double precision,dimension(3,3,3)::dpssints
double precision,dimension(3,3,3,3)::ddssints,dppsints
double precision,dimension(3,3,3,3,3)::ddpsints,dsdpints
double precision,dimension(3,3,3,3,3,3)::dddsints

do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellpd)
wiipi(i)=wvector(i)-qvector(i)
end do

call twoe_ddds(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,pvector,qvector,wvector,dddsints)
dddsintsm_global=dddsints
do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddpints(i,j,k,l,m,n,io)=qkkck(io)*dddsints(i,j,k,l,m,n)
end do
end do
end do
end do
end do
end do
end do

factor=1.0d0/(2.0d0*zc_plus_zd)
factor2=factor*(za_plus_zb)/zabcd
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddpints(i,j,k,l,m,n,m)=dddpints(i,j,k,l,m,n,m)+factor*ddpsintsm_global(i,j,k,l,n)-factor2*ddpsintsmp1_global(i,j,k,l,n)
end do
end do
end do
end do
end do
end do


do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddpints(i,j,k,l,m,n,n)=dddpints(i,j,k,l,m,n,n)+factor*ddpsintsm_global(i,j,k,l,m)-factor2*ddpsintsmp1_global(i,j,k,l,m)
end do
end do
end do
end do
end do
end do


call twoe_dsdp(mindex+1,zabcd,zc_plus_zd,za_plus_zb,ishelldc,ishelldb,ishellda,qvector,pvector,wvector,dsdpints)

factor=1.0d0/(2.0d0*zabcd)
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddpints(i,j,k,l,m,n,i)=dddpints(i,j,k,l,m,n,i)+factor*dsdpints(m,n,k,l,j)
end do
end do
end do
end do
end do
end do


do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddpints(i,j,k,l,m,n,j)=dddpints(i,j,k,l,m,n,j)+factor*dsdpints(m,n,k,l,i)
end do
end do
end do
end do
end do
end do


call twoe_dsdp(mindex+1,zabcd,zc_plus_zd,za_plus_zb,ishelldc,ishellda,ishelldb,qvector,pvector,wvector,dsdpints)


do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddpints(i,j,k,l,m,n,k)=dddpints(i,j,k,l,m,n,k)+factor*dsdpints(m,n,i,j,l)
end do
end do
end do
end do
end do
end do

do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddpints(i,j,k,l,m,n,l)=dddpints(i,j,k,l,m,n,l)+factor*dsdpints(m,n,i,j,k)
end do
end do
end do
end do
end do
end do





call twoe_ddds(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,pvector,qvector,wvector,dddsints)
dddsintsmp1_global=dddsints

do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dddpints(i,j,k,l,m,n,io)=dddpints(i,j,k,l,m,n,io)+wiipi(io)*dddsints(i,j,k,l,m,n)
end do
end do
end do
end do
end do
end do
end do




return
end


