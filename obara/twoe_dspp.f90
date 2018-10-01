subroutine twoe_dspp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpc,ishellpd,pvector,qvector,wvector,dsppints)
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
integer,intent(in)::ishelld,ishellpd,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dsppints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::dsssints,pspsints
double precision,dimension(3,3,3)::dspsints



do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellpc)
wiipi(i)=wvector(i)-qvector(i)
end do


call twoe_dsps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpd,pvector,qvector,wvector,dspsints)


do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsppints(i,j,k,l)=qkkck(k)*dspsints(i,j,l)
end do
end do
end do
end do

factor=1.0/(2.0d0*zc_plus_zd)
factor2=factor*za_plus_zb/zabcd
do k=1,3
do j=1,3
do i=1,3
dsppints(i,j,k,k)=dsppints(i,j,k,k)+factor*dsssm_global(i,j)-factor2*dsssmp1_global(i,j)
end do
end do
end do


call twoe_psps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpd,pvector,qvector,wvector,pspsints)


factor=1.0/(2.0d0*zabcd)
do l=1,3
do j=1,3
do i=1,3
dsppints(i,j,i,l)=dsppints(i,j,i,l)+factor*pspsints(j,l)
end do
end do
end do


do l=1,3
do j=1,3
do i=1,3
dsppints(i,j,j,l)=dsppints(i,j,j,l)+factor*pspsints(i,l)
end do
end do
end do



call twoe_dsps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpd,pvector,qvector,wvector,dspsints)
dspsintsmp1_global=dspsints

do l=1,3
do k=1,3
do j=1,3
do i=1,3
dsppints(i,j,k,l)=dsppints(i,j,k,l)+wiipi(k)*dspsints(i,j,l)
end do
end do
end do
end do





return
end


