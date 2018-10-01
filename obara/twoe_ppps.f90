subroutine twoe_ppps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,ishellpc,pvector,qvector,wvector,pppsints)
use gaussian_basis
implicit double precision (a-h,o-z)
interface
subroutine twoe_psss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellp,pvector,qvector,wvector,psssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:),intent(inout)::psssints
end subroutine twoe_psss

subroutine twoe_ppss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,pvector,qvector,wvector,ppssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::ppssints
end subroutine twoe_ppss


end interface



integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,intent(inout),dimension(3,3,3)::pppsints
double precision,dimension(3)::qkkck,wkkqk,psssints
double precision,dimension(3,3)::ppssints

pi=dacos(-1.0d0)



do i=1,3
qkkck(i)=qvector(i)-shell_coor(i,ishellpc)
wkkqk(i)=wvector(i)-qvector(i)
end do

call twoe_ppss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,pvector,qvector,wvector,ppssints)

ppssints_global=ppssints

do k=1,3
do j=1,3
do i=1,3
pppsints(i,j,k)=qkkck(k)*ppssints(i,j)
end do
end do
end do


! add in delta(jk) term
factor=1.0d0/(2.0d0*zabcd)
do j=1,3
do i=1,3
pppsints(i,j,j)=pppsints(i,j,j)+ factor * psssints_mp1_global(i)
end do
end do


call twoe_ppss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,pvector,qvector,wvector,ppssints)
ppssintsmp1_global=ppssints
do k=1,3
do j=1,3
do i=1,3
pppsints(i,j,k)=pppsints(i,j,k)+wkkqk(k)*ppssints(i,j)
end do
end do
end do


call twoe_psss(mindex+1,zabcd,za_plus_zb,zc_plus_zd,ishellpb,pvector,qvector,wvector,psssints)

! add in delta(ik) term
do j=1,3
do i=1,3
pppsints(i,j,i)=pppsints(i,j,i)+factor*psssints(j)
end do
end do






return
end


