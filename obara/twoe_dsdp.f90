subroutine twoe_dsdp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,ishellpd,pvector,qvector,wvector,dsdpints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface

subroutine twoe_dsds(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,pvector,qvector,wvector,dsdsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dsdsints
end subroutine twoe_dsds

subroutine twoe_dsps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellp,pvector,qvector,wvector,dspsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:),intent(inout)::dspsints
end subroutine twoe_dsps



end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:),intent(inout)::dsdpints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::dsssints,dsssintsmp1,pspsints
double precision,dimension(3,3,3)::dspsints
double precision,dimension(3,3,3,3)::dsdsints


do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellpd)
wiipi(i)=wvector(i)-qvector(i)
end do

call twoe_dsds(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,pvector,qvector,wvector,dsdsints)


do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsdpints(i,j,k,l,m)=qkkck(m)*dsdsints(i,j,k,l)
end do
end do
end do
end do
end do


factor=1.0d0/(2.0d0*zc_plus_zd)
factor2=factor*za_plus_zb/zabcd
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsdpints(i,j,k,l,k)=dsdpints(i,j,k,l,k)+factor*dspsintsm_global(i,j,l)-factor2*dspsintsmp1_global(i,j,l)
end do
end do
end do
end do


do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsdpints(i,j,k,l,l)=dsdpints(i,j,k,l,l)+factor*dspsintsm_global(i,j,k)-factor2*dspsintsmp1_global(i,j,k)
end do
end do
end do
end do

call twoe_dsds(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,pvector,qvector,wvector,dsdsints)
dsdsmp1_global=dsdsints

do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsdpints(i,j,k,l,m)=dsdpints(i,j,k,l,m)+wiipi(m)*dsdsints(i,j,k,l)
end do
end do
end do
end do
end do




!here we must switch za_plus_zb with zc_plus_zd and pvector with qvector 

call twoe_dsps(mindex+1,zabcd,zc_plus_zd,za_plus_zb,ishelldc,ishellda,qvector,pvector,wvector,dspsints)


factor=1.0d0/(2.0d0*zabcd)
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsdpints(i,j,k,l,i)=dsdpints(i,j,k,l,i)+factor*dspsints(k,l,j)
end do
end do
end do
end do 


do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsdpints(i,j,k,l,j)=dsdpints(i,j,k,l,j)+factor*dspsints(k,l,i)
end do
end do
end do
end do








return
end


