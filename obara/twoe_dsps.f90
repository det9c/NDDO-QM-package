subroutine twoe_dsps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellp,pvector,qvector,wvector,dspsints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface
subroutine twoe_dsss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,pvector,qvector,wvector,dsssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::dsssints
end subroutine twoe_dsss
end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:),intent(inout)::dspsints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::dsssints

pi=dacos(-1.0d0)

do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellp)
wiipi(i)=wvector(i)-qvector(i)
end do

call twoe_dsss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,pvector,qvector,wvector,dsssints)
dsssm_global=dsssints

do k=1,3
do j=1,3
do i=1,3
dspsints(i,j,k)=qkkck(k)*dsssints(i,j)
end do
end do
end do

call twoe_dsss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelld,pvector,qvector,wvector,dsssints)
dsssmp1_global=dsssints
do k=1,3
do j=1,3
do i=1,3
dspsints(i,j,k)=dspsints(i,j,k)+wiipi(k)*dsssints(i,j)
end do
end do
end do


factor=1.0d0/(2.0d0*zabcd)
do j=1,3
do i=1,3
dspsints(i,j,i)=dspsints(i,j,i)+factor*psssm_global(j)
end do
end do

do j=1,3
do i=1,3
dspsints(i,j,j)=dspsints(i,j,j)+factor*psssm_global(i)
end do
end do




return
end


