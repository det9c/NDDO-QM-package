subroutine projector(hessian,i3n)
use tables
use constants
implicit double precision (a-h,o-z)


interface
subroutine outerp(n,alpha,v1,v2,mat)
double precision,dimension(n,1),intent(in)::v1,v2
double precision,dimension(n,n),intent(inout)::mat
double precision,intent(in)::alpha
integer,intent(in)::n
end subroutine outerp
end interface


integer,intent(in)::i3n
double precision,dimension(i3n,i3n),intent(inout)::hessian
double precision,dimension(i3n)::xvec,yvec,zvec
double precision,dimension(i3n,i3n)::scratch,p

! compute translational eigenvectors
xvec=zero
yvec=zero
zvec=zero

do i=1,numat
istart=3*i-2
xvec(istart)=dsqrt( atmass( zeff(i) ) ) 
yvec(istart+1)=dsqrt( atmass( zeff(i) ) )
zvec(istart+2)=dsqrt( atmass( zeff(i) ) )
end do

xnorm=zero
ynorm=zero
znorm=zero

! normalize vectors
do i=1,i3n
xnorm=xnorm+xvec(i)**2
ynorm=ynorm+yvec(i)**2
znorm=znorm+zvec(i)**2
end do

xvec=xvec/dsqrt(xnorm)
yvec=yvec/dsqrt(ynorm)
zvec=zvec/dsqrt(znorm)
!compute outer product
p=zero
call outerp(i3n,1.0d0,xvec,xvec,p)
call outerp(i3n,1.0d0,yvec,yvec,p)
call outerp(i3n,1.0d0,zvec,zvec,p)



do i=1,numat
istart=3*i-2
zz=dsqrt(atmass( zeff(i) ))/.529177d0 ! convert to bohr since hessian in au
xvec(istart)=0.0d0
xvec(istart+1)=-z(i)*zz
xvec(istart+2)=y(i)*zz
yvec(istart)=z(i)*zz
yvec(istart+1)=0.0d0
yvec(istart+2)=-x(i)*zz
zvec(istart)=-y(i)*zz
zvec(istart+1)=x(i)*zz
zvec(istart+2)=0.0d0
end do

xnorm=zero
ynorm=zero
znorm=zero

! normalize vectors
do i=1,i3n
xnorm=xnorm+xvec(i)**2
ynorm=ynorm+yvec(i)**2
znorm=znorm+zvec(i)**2
end do

xvec=xvec/dsqrt(xnorm)
yvec=yvec/dsqrt(ynorm)
zvec=zvec/dsqrt(znorm)

call outerp(i3n,1.0d0,xvec,xvec,p)
call outerp(i3n,1.0d0,yvec,yvec,p)
call outerp(i3n,1.0d0,zvec,zvec,p)

p=-1.0d0*p
do i=1,i3n
p(i,i)=1.0d0+p(i,i)
end do


call dgemm( 'T', 'N', i3n,i3n,i3n, one,p &
 , i3n,hessian , i3n,zero,scratch,i3n )
call dgemm( 'N', 'N', i3n,i3n,i3n,one,scratch &
 , i3n,p , i3n,zero,hessian,i3n )

return
end subroutine projector





