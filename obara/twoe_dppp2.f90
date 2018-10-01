subroutine twoe_dppp2(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpc,ishellpd,pvector,qvector,wvector,dpppints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface
subroutine twoe_dspp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpc,ishellpd,pvector,qvector,wvector,dsppints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellpd,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dsppints
end subroutine twoe_dspp

subroutine twoe_ppps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,ishellpc,pvector,qvector,wvector,pppsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,intent(inout),dimension(3,3,3)::pppsints
end subroutine twoe_ppps

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
integer,intent(in)::ishelld,ishellpd,ishellpc,ishellpb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:),intent(inout)::dpppints
double precision,dimension(3)::piiai,wiipi,psssints
double precision,dimension(3,3)::dsssints,pspsints
double precision,dimension(3,3,3,3)::dppsints,dsppints
double precision,dimension(3,3,3)::pppsints,dspsints


do i=1,3
piiai(i)=pvector(i)-shell_coor(i,ishellpb)
wiipi(i)=wvector(i)-pvector(i)
end do


call twoe_dspp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpc,ishellpd,pvector,qvector,wvector,dsppints)

do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,l,m)=piiai(k)*dsppints(i,j,l,m)
end do
end do
end do
end do
end do


factor=1.0/(2.0d0*zabcd)
do m=1,3
do k=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,k,m)=dpppints(i,j,k,k,m)+factor*dspsintsmp1_global(i,j,m)
end do
end do
end do
end do




call twoe_ppps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpc,ishellpd,pvector,qvector,wvector,pppsints)

factor=1.0/(2.0d0*za_plus_zb)
factor2=factor*zc_plus_zd/zabcd
do m=1,3
do l=1,3
do j=1,3
do i=1,3
dpppints(i,j,i,l,m)=dpppints(i,j,i,l,m)+factor*pppsints(j,l,m)
end do
end do
end do
end do


do m=1,3
do l=1,3
do j=1,3
do i=1,3
dpppints(i,j,j,l,m)=dpppints(i,j,j,l,m)+factor*pppsints(i,l,m)
end do
end do
end do
end do


call twoe_ppps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpc,ishellpd,pvector,qvector,wvector,pppsints)


do m=1,3
do l=1,3
do j=1,3
do i=1,3
dpppints(i,j,i,l,m)=dpppints(i,j,i,l,m)-factor2*pppsints(j,l,m)
end do
end do
end do
end do

do m=1,3
do l=1,3
do j=1,3
do i=1,3
dpppints(i,j,j,l,m)=dpppints(i,j,j,l,m)-factor2*pppsints(i,l,m)
end do
end do
end do
end do


call twoe_dsps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpc,pvector,qvector,wvector,dspsints)

factor=1.0/(2.0d0*zabcd)
do k=1,3
do l=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,l,k)=dpppints(i,j,k,l,k)+factor*dspsints(i,j,l)
end do
end do
end do
end do



call twoe_dspp(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpc,ishellpd,pvector,qvector,wvector,dsppints)

do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,l,m)=dpppints(i,j,k,l,m)+wiipi(k)*dsppints(i,j,l,m)
end do
end do
end do
end do
end do






return
end


