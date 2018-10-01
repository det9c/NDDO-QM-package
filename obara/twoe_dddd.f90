subroutine twoe_dddd(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishelldd,pvector,qvector,wvector,ddddints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface
subroutine twoe_dddp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishellpd,pvector,qvector,wvector,dddpints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:,:),intent(inout)::dddpints
end subroutine twoe_dddp

subroutine twoe_ddpp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishellpc,ishellpd,pvector,qvector,wvector,ddppints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishellpc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:),intent(inout)::ddppints
end subroutine twoe_ddpp

subroutine twoe_dpdp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishellpb,ishelldc,ishellpd,pvector,qvector,wvector,dpdpints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishellpb,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:),intent(inout)::dpdpints
end subroutine twoe_dpdp


end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishelldc,ishelldd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(3,3,3,3,3,3,9),intent(inout)::ddddints
double precision,dimension(3)::qkkck,wiipi,psssints
double precision,dimension(3,3)::ppssints,ppssintsmp1
double precision,dimension(3,3,3)::dpssints
double precision,dimension(3,3,3,3)::ddssints,dppsints
double precision,dimension(3,3,3,3,3)::ddpsints,dsdpints
double precision,dimension(3,3,3,3,3,3)::ddppintsm,ddppintsmp1,dpdpints
double precision,dimension(3,3,3,3,3,3,3)::dddpints


do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishelldd)
wiipi(i)=wvector(i)-qvector(i)
end do

call twoe_dddp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishelldd,pvector,qvector,wvector,dddpints)

do ip=1,3
do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddddints(i,j,k,l,m,n,isquish(io,ip))=qkkck(ip)*dddpints(i,j,k,l,m,n,io)
end do
end do
end do
end do
end do
end do
end do
end do



factor=1.0d0/(2.0d0*zc_plus_zd)
factor2=factor*(za_plus_zb)/zabcd
do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddddints(i,j,k,l,m,n,isquish(io,io))=ddddints(i,j,k,l,m,n,isquish(io,io))+factor*dddsintsm_global(i,j,k,l,m,n)-&
factor2*dddsintsmp1_global(i,j,k,l,m,n)
end do
end do
end do
end do
end do
end do
end do

call twoe_dddp(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishelldd,pvector,qvector,wvector,dddpints)

do ip=1,3
do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddddints(i,j,k,l,m,n,isquish(io,ip))=ddddints(i,j,k,l,m,n,isquish(io,ip))+wiipi(ip)*dddpints(i,j,k,l,m,n,io)
end do
end do
end do
end do
end do
end do
end do
end do


call twoe_ddpp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishelldd,pvector,qvector,wvector,ddppintsm)
call twoe_ddpp(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishelldd,pvector,qvector,wvector,ddppintsmp1)

do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddddints(i,j,k,l,m,n,isquish(io,m))=ddddints(i,j,k,l,m,n,isquish(io,m))+factor*ddppintsm(i,j,k,l,n,io)-factor2*ddppintsmp1(i,j,k,l,n,io)
end do
end do
end do
end do
end do
end do
end do

do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddddints(i,j,k,l,m,n,isquish(io,n))=ddddints(i,j,k,l,m,n,isquish(io,n))+factor*ddppintsm(i,j,k,l,m,io)-factor2*ddppintsmp1(i,j,k,l,m,io)
end do
end do
end do
end do
end do
end do
end do


factor=1.0d0/(2.0d0*zabcd)

call twoe_dpdp(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishelldd,pvector,qvector,wvector,dpdpints)


do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddddints(i,j,k,l,m,n,isquish(io,k))=ddddints(i,j,k,l,m,n,isquish(io,k))+factor*dpdpints(i,j,l,m,n,io)
end do
end do
end do
end do
end do
end do
end do

do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddddints(i,j,k,l,m,n,isquish(io,l))=ddddints(i,j,k,l,m,n,isquish(io,l))+factor*dpdpints(i,j,k,m,n,io)
end do
end do
end do
end do
end do
end do
end do







call twoe_dpdp(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishelldb,ishellda,ishelldc,ishelldd,pvector,qvector,wvector,dpdpints)


do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddddints(i,j,k,l,m,n,isquish(io,i))=ddddints(i,j,k,l,m,n,isquish(io,i))+factor*dpdpints(k,l,j,m,n,io)
end do
end do
end do
end do
end do
end do
end do

do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do j=1,3
do i=1,3
ddddints(i,j,k,l,m,n,isquish(io,j))=ddddints(i,j,k,l,m,n,isquish(io,j))+factor*dpdpints(k,l,i,m,n,io)
end do
end do
end do
end do
end do
end do
end do










return
end


