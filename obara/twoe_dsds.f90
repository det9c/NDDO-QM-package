subroutine twoe_dsds(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,pvector,qvector,wvector,dsdsints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface

subroutine twoe_dsps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellp,pvector,qvector,wvector,dspsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:),intent(inout)::dspsints
end subroutine twoe_dsps

subroutine twoe_psps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,pvector,qvector,wvector,pspsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::pspsints
end subroutine twoe_psps




end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dsdsints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::dsssints,dsssintsmp1,pspsints
double precision,dimension(3,3,3)::dspsints

pi=dacos(-1.0d0)

do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishelldc)
wiipi(i)=wvector(i)-qvector(i)
end do

call twoe_dsps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,pvector,qvector,wvector,dspsints)
dspsintsm_global=dspsints

do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsdsints(i,j,k,l)=qkkck(l)*dspsints(i,j,k)
end do
end do
end do
end do



factor=1.0d0/(2.0d0*zc_plus_zd)
factor2=factor*za_plus_zb/zabcd
do k=1,3
do j=1,3
do i=1,3
dsdsints(i,j,k,k)=dsdsints(i,j,k,k)+factor*dsssm_global(i,j)-factor2*dsssmp1_global(i,j)
end do
end do
end do


call twoe_psps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,pvector,qvector,wvector,pspsints)


factor=1.0d0/(2.0d0*zabcd)
do k=1,3
do j=1,3
do i=1,3
dsdsints(i,j,k,i)=dsdsints(i,j,k,i)+factor*pspsints(j,k)
end do
end do
end do

do k=1,3
do j=1,3
do i=1,3
dsdsints(i,j,k,j)=dsdsints(i,j,k,j)+factor*pspsints(i,k)
end do
end do
end do




call twoe_dsps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,pvector,qvector,wvector,dspsints)
dspsintsmp1_global=dspsints

do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsdsints(i,j,k,l)=dsdsints(i,j,k,l)+wiipi(l)*dspsints(i,j,k)
end do
end do
end do
end do




return
end


