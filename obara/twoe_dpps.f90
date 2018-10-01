subroutine twoe_dpps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpc,pvector,qvector,wvector,dppsints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface
subroutine twoe_ppss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,pvector,qvector,wvector,ppssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::ppssints
end subroutine twoe_ppss

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
integer,intent(in)::ishelld,ishellpb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dppsints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::dsssints,ppssints
double precision,dimension(3,3,3)::dpssints



do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellpc)
wiipi(i)=wvector(i)-qvector(i)
end do


call twoe_dpss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,pvector,qvector,wvector,dpssints)
dpssm_global=dpssints
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dppsints(i,j,k,l)=qkkck(l)*dpssints(i,j,k)
end do
end do
end do
end do

factor=1.0/(2.0d0*zabcd)
do k=1,3
do j=1,3
do i=1,3
dppsints(i,j,k,k)=dppsints(i,j,k,k)+factor*dsssmp1_global(i,j)
end do
end do
end do

call twoe_ppss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,pvector,qvector,wvector,ppssints)


do k=1,3
do j=1,3
do i=1,3
dppsints(i,j,k,i)=dppsints(i,j,k,i)+factor*ppssints(j,k)
end do
end do
end do

do k=1,3
do j=1,3
do i=1,3
dppsints(i,j,k,j)=dppsints(i,j,k,j)+factor*ppssints(i,k)
end do
end do
end do


call twoe_dpss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,pvector,qvector,wvector,dpssints)
dpssmp1_global=dpssints
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dppsints(i,j,k,l)=dppsints(i,j,k,l)+wiipi(l)*dpssints(i,j,k)
end do
end do
end do
end do



return
end


