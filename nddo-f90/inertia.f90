subroutine inertia(reori)
! reori=true means to put in principle axis orientation, otherwise
! just compute moments of inertia and return
use tables
use constants
implicit double precision(a-h,o-z)
logical,intent(in)::reori
double precision,dimension(3,3)::tensor,vectors,tvec
double precision,dimension(numat)::xold,yold,zold
convert=16.85803902D0 ! conversion from amu*ang**2 to cm(-1)

totmass=zero
sum1=zero
sum2=zero
sum3=zero
do i=1,numat
atmssi=atmass(zeff(i))
totmass=totmass+atmssi
sum1=sum1+x(i)*atmssi
sum2=sum2+y(i)*atmssi
sum3=sum3+z(i)*atmssi
end do
sum1=sum1/totmass
sum2=sum2/totmass
sum3=sum3/totmass




do i=1,numat
x(i)=x(i)-sum1
y(i)=y(i)-sum2
z(i)=z(i)-sum3
end do

tensor=zero
do i=1,numat
zcharge=atmass(zeff(i))
tensor(1,1)=tensor(1,1)+zcharge* ( y(i)**2 + z(i)**2)
tensor(2,2)=tensor(2,2)+zcharge* ( x(i)**2 + z(i)**2)
tensor(3,3)=tensor(3,3)+zcharge* ( x(i)**2 + y(i)**2)
tensor(2,1)=tensor(2,1)-zcharge*x(i)*y(i)
tensor(3,1)=tensor(3,1)-zcharge*x(i)*z(i)
tensor(3,2)=tensor(3,2)-zcharge*z(i)*y(i)
end do
tensor(1,2)=tensor(2,1)
tensor(1,3)=tensor(3,1)
tensor(2,3)=tensor(3,2)
call eig(tensor,vectors,3,3,0)
!call matprt(vectors,3,3,3,3)

t1=tensor(1,1)
t2=tensor(2,2)
t3=tensor(3,3)

rbohr=1.0d0/.52917**2
if(.not. reori)then
write(*,*)''
write(*,10)'Moments of Inertia in AMU*BOHR**2::',tensor(1,1)*rbohr,tensor(2,2)*rbohr,tensor(3,3)*rbohr
write(*,10)'Rotational Constants (cm-1)::',convert /tensor(1,1),convert/tensor(2,2),convert/tensor(3,3)
return
end if

xold=x
yold=y
zold=z


tvec=transpose(vectors)
do i=1,numat
x(i)=tvec(3,1)*xold(i)+tvec(3,2)*yold(i)+tvec(3,3)*zold(i)
y(i)=tvec(2,1)*xold(i)+tvec(2,2)*yold(i)+tvec(2,3)*zold(i)
z(i)=tvec(1,1)*xold(i)+tvec(1,2)*yold(i)+tvec(1,3)*zold(i)
end do





tensor=zero
do i=1,numat
zcharge=atmass(zeff(i))
tensor(1,1)=tensor(1,1)+zcharge* ( y(i)**2 + z(i)**2)
tensor(2,2)=tensor(2,2)+zcharge* ( x(i)**2 + z(i)**2)
tensor(3,3)=tensor(3,3)+zcharge* ( x(i)**2 + y(i)**2)
tensor(2,1)=tensor(2,1)-zcharge*x(i)*y(i)
tensor(3,1)=tensor(3,1)-zcharge*x(i)*z(i)
tensor(3,2)=tensor(3,2)-zcharge*z(i)*y(i)
end do
tensor(1,2)=tensor(2,1)
tensor(1,3)=tensor(3,1)
tensor(2,3)=tensor(3,2)
call eig(tensor,vectors,3,3,0)
!call matprt(tensor,3,3,3,3)


!print*,'Principle axis orientation used in vibrational frequency calculation'
!do i=1,numat
!write(*,20)zeff(i),x(i),y(i),z(i)
!end do

t1=abs(tensor(1,1)-t1)
t2=abs(tensor(2,2)-t2)
t3=abs(tensor(3,3)-t3)
err=1D-8
!if( (t1.gt.err) .or. (t2.gt.err) .or. (t3.gt.err) )then
!write(*,*)'Moments of inertia in principle axis coordinate system &
!differ by more than ',err,' when using input orientation.  Something is &
!wrong with principle axis transformation!! '
!stop
!end if

write(*,10)'Moments of Inertia in AMU*BOHR**2',tensor(1,1)*rbohr,tensor(2,2)*rbohr,tensor(3,3)*rbohr
write(*,10)'Rotational Constants (cm-1)',convert /tensor(1,1),convert/tensor(2,2),convert/tensor(3,3)




10 format(A,f10.5,f10.5,f10.5)
20 format(i3,f10.5,f10.5,f10.5)

end subroutine inertia
