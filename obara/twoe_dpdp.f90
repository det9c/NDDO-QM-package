subroutine twoe_dpdp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishellpb,ishelldc,ishellpd,pvector,qvector,wvector,dpdpints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface

subroutine twoe_dsdp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,ishellpd,pvector,qvector,wvector,dsdpints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:),intent(inout)::dsdpints
end subroutine twoe_dsdp

subroutine twoe_dpps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpc,pvector,qvector,wvector,dppsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellpb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dppsints
end subroutine twoe_dpps


subroutine twoe_dspp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpc,ishellpd,pvector,qvector,wvector,dsppints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellpd,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dsppints
end subroutine twoe_dspp



end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishellpb,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:),intent(inout)::dpdpints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::ppssints,ppssintsmp1
double precision,dimension(3,3,3)::dpssints
double precision,dimension(3,3,3,3)::ddssints,dppsints,dppsintsmp1,dsppints
double precision,dimension(3,3,3,3,3)::ddpsints,dsdpints

do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellpd)
wiipi(i)=wvector(i)-qvector(i)
end do

! need to swap bra and ket 
call twoe_dsdp(mindex,zabcd,zc_plus_zd,za_plus_zb,ishelldc,ishellda,ishellpb,qvector,pvector,wvector,dsdpints)


do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpdpints(i,j,k,l,m,n)=qkkck(n)*dsdpints(l,m,i,j,k)
end do
end do
end do
end do
end do
end do



! dsds for last term is swapped b.c of reversal above

factor=1.0d0/(2.0d0*zabcd)
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpdpints(i,j,k,l,m,k)=dpdpints(i,j,k,l,m,k)+factor*dsdsmp1_global(l,m,i,j)
end do
end do
end do
end do
end do



call twoe_dsdp(mindex+1,zabcd,zc_plus_zd,za_plus_zb,ishelldc,ishellda,ishellpb,qvector,pvector,wvector,dsdpints)


do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpdpints(i,j,k,l,m,n)=dpdpints(i,j,k,l,m,n)+wiipi(n)*dsdpints(l,m,i,j,k)
end do
end do
end do
end do
end do
end do



call twoe_dpps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishellpb,ishelldc,pvector,qvector,wvector,dppsints)
call twoe_dpps(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishellpb,ishelldc,pvector,qvector,wvector,dppsintsmp1)

factor=1.0/(2.0d0*zc_plus_zd)
factor2=factor*za_plus_zb/zabcd
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpdpints(i,j,k,l,m,l)=dpdpints(i,j,k,l,m,l)+factor*dppsints(i,j,k,m)-factor2*dppsintsmp1(i,j,k,m)
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
dpdpints(i,j,k,l,m,m)=dpdpints(i,j,k,l,m,m)+factor*dppsints(i,j,k,l)-factor2*dppsintsmp1(i,j,k,l)
end do
end do
end do
end do
end do


! swap on this call

call twoe_dspp(mindex+1,zabcd,zc_plus_zd,za_plus_zb,ishelldc,ishellda,ishellpb,qvector,pvector,wvector,dsppints)

factor=1.0d0/(2.0d0*zabcd)
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
dpdpints(i,j,k,l,m,i)=dpdpints(i,j,k,l,m,i)+factor*dsppints(l,m,j,k)
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
dpdpints(i,j,k,l,m,j)=dpdpints(i,j,k,l,m,j)+factor*dsppints(l,m,i,k)
end do
end do
end do
end do
end do





return
end


