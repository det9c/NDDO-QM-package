subroutine twoe_dppp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpc,ishellpd,pvector,qvector,wvector,dpppints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface
subroutine twoe_dpps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpc,pvector,qvector,wvector,dppsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellpb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dppsints
end subroutine twoe_dpps

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
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::dsssints,pspsints
double precision,dimension(3,3,3,3)::dppsints
double precision,dimension(3,3,3)::pppsints,dspsints

dpppints=0.

do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellpc)
wiipi(i)=wvector(i)-qvector(i)
end do


call twoe_dpps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpd,pvector,qvector,wvector,dppsints)

do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,l,m)=qkkck(l)*dppsints(i,j,k,m)
end do
end do
end do
end do
end do


factor=1.0/(2.0d0*zc_plus_zd)
factor2=factor*za_plus_zb/zabcd
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,l,l)=dpppints(i,j,k,l,l)+factor*dpssm_global(i,j,k)-factor2*dpssmp1_global(i,j,k)
end do
end do
end do
end do


call twoe_ppps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpd,pvector,qvector,wvector,pppsints)


factor=1.0/(2.0d0*zabcd)
do m=1,3
do k=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,i,m)=dpppints(i,j,k,i,m)+factor*pppsints(j,k,m)
end do
end do
end do
end do




do m=1,3
do k=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,j,m)=dpppints(i,j,k,j,m)+factor*pppsints(i,k,m)
end do
end do
end do
end do




call twoe_dsps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpd,pvector,qvector,wvector,dspsints)





do m=1,3
do k=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,k,m)=dpppints(i,j,k,k,m)+factor*dspsints(i,j,m)
end do
end do
end do
end do








call twoe_dpps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpd,pvector,qvector,wvector,dppsints)

do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpppints(i,j,k,l,m)=dpppints(i,j,k,l,m)+wiipi(l)*dppsints(i,j,k,m)
end do
end do
end do
end do
end do




return
end


