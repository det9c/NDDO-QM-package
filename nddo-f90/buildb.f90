subroutine buildb
use constants
use tables
use scratch_array
! this will build the B matrix from Wilson, Decius,Cross
!
implicit double precision (a-h,o-z)
double precision,dimension(3)::eji,ejk,eij,ekl
integer,allocatable,dimension(:)::indx
double precision,allocatable,dimension(:,:)::binv

if(allocated(bmatrix))deallocate(bmatrix)
if(allocated(bmat4int))deallocate(bmat4int)
numcart=3*ninput
numint=numcart-6
if(numcart==6)numint=1
allocate(bmatrix(numint,numcart))
allocate(bmat4int(numint,numcart))
bmatrix=zero

! put in all the derivatives wrt bond lengths

irowb=0
do iatom=2,ninput
irowb=irowb+1
jatom=ref(iatom,1) ! atom to which atom iatom is bonded 
rij=(x(jatom)-x(iatom))**2+(y(jatom)-y(iatom))**2+(z(jatom)-z(iatom))**2
rij=dsqrt(rij)
icol=3*iatom-2
jcol=3*jatom-2
bmatrix(irowb,icol)=(x(iatom)-x(jatom))/rij
bmatrix(irowb,icol+1)=(y(iatom)-y(jatom))/rij
bmatrix(irowb,icol+2)=(z(iatom)-z(jatom))/rij
bmatrix(irowb,jcol)=-bmatrix(irowb,icol)
bmatrix(irowb,jcol+1)=-bmatrix(irowb,icol+1)
bmatrix(irowb,jcol+2)=-bmatrix(irowb,icol+2)
end do

! put in derivatives for angles
do iatom=3,ninput
irowb=irowb+1
jatom=ref(iatom,1)
katom=ref(iatom,2)
rji=(x(jatom)-x(iatom))**2+(y(jatom)-y(iatom))**2+(z(jatom)-z(iatom))**2
rjk=(x(jatom)-x(katom))**2+(y(jatom)-y(katom))**2+(z(jatom)-z(katom))**2
rji=dsqrt(rji)
rjk=dsqrt(rjk)
eji(1)=(x(iatom)-x(jatom))/rji
eji(2)=(y(iatom)-y(jatom))/rji
eji(3)=(z(iatom)-z(jatom))/rji
ejk(1)=(x(katom)-x(jatom))/rjk
ejk(2)=(y(katom)-y(jatom))/rjk
ejk(3)=(z(katom)-z(jatom))/rjk
icol=3*iatom-2
jcol=3*jatom-2
kcol=3*katom-2
phi=q(iatom,2)
cosphi=dcos(phi)
sinphi=dsin(phi)
bmatrix(irowb,icol)=( cosphi*eji(1)-ejk(1) ) / ( rji *sinphi)
bmatrix(irowb,icol+1)=( cosphi*eji(2)-ejk(2) ) / ( rji *sinphi)
bmatrix(irowb,icol+2)=( cosphi*eji(3)-ejk(3) ) / ( rji *sinphi)
bmatrix(irowb,jcol)=(  (rji-rjk*cosphi)*eji(1) + (rjk-rji*cosphi)*ejk(1) ) / (rji*rjk*sinphi)
bmatrix(irowb,jcol+1)=(  (rji-rjk*cosphi)*eji(2) + (rjk-rji*cosphi)*ejk(2) ) / (rji*rjk*sinphi)
bmatrix(irowb,jcol+2)=(  (rji-rjk*cosphi)*eji(3) + (rjk-rji*cosphi)*ejk(3) ) / (rji*rjk*sinphi)
bmatrix(irowb,kcol)=( cosphi*ejk(1)-eji(1) ) / ( rjk *sinphi)
bmatrix(irowb,kcol+1)=( cosphi*ejk(2)-eji(2) ) / ( rjk *sinphi)
bmatrix(irowb,kcol+2)=( cosphi*ejk(3)-eji(3) ) / ( rjk *sinphi)
end do


if(ninput.gt.3)then
do iatom=4,ninput
irowb=irowb+1
jatom=ref(iatom,1)
katom=ref(iatom,2)
latom=ref(iatom,3)
rij=(x(jatom)-x(iatom))**2+(y(jatom)-y(iatom))**2+(z(jatom)-z(iatom))**2
rjk=(x(jatom)-x(katom))**2+(y(jatom)-y(katom))**2+(z(jatom)-z(katom))**2
rkl=(x(latom)-x(katom))**2+(y(latom)-y(katom))**2+(z(latom)-z(katom))**2
rij=dsqrt(rij)
rjk=dsqrt(rjk)
rkl=dsqrt(rkl)
eij(1)=(x(jatom)-x(iatom))/rij
eij(2)=(y(jatom)-y(iatom))/rij
eij(3)=(z(jatom)-z(iatom))/rij
ejk(1)=(x(katom)-x(jatom))/rjk
ejk(2)=(y(katom)-y(jatom))/rjk
ejk(3)=(z(katom)-z(jatom))/rjk
ekl(1)=(x(latom)-x(katom))/rkl
ekl(2)=(y(latom)-y(katom))/rkl
ekl(3)=(z(latom)-z(katom))/rkl
icol=3*iatom-2
jcol=3*jatom-2
kcol=3*katom-2
lcol=3*latom-2
phi=q(iatom,2)
sinphi=dsin(phi)
cosphi=dcos(phi)
cosphi3=  -ekl(1)*ejk(1) - ekl(2)*ejk(2)- ekl(3)*ejk(3)
angle=acos(cosphi3)
sinphi3=dsin(angle)

bmatrix(irowb,icol)=-( eij(2)*ejk(3)-eij(3)*ejk(2) )/( rij * sinphi*sinphi) 
bmatrix(irowb,icol+1)=-( eij(3)*ejk(1)-eij(1)*ejk(3) )/( rij * sinphi*sinphi)
bmatrix(irowb,icol+2)=-( eij(1)*ejk(2)-eij(2)*ejk(1) )/( rij * sinphi*sinphi)

bmatrix(irowb,jcol)= (  rjk - rij * cosphi) * ( eij(2)*ejk(3)-eij(3)*ejk(2) ) / (rjk*rij*sinphi*sinphi) &
+ cosphi3 * (  ekl(2)*ejk(3)-ekl(3)*ejk(2) ) / ( rjk*sinphi3*sinphi3)
bmatrix(irowb,jcol+1)= (  rjk - rij * cosphi) * ( eij(3)*ejk(1)-eij(1)*ejk(3) ) / (rjk*rij*sinphi*sinphi) &
+ cosphi3 * (  ekl(3)*ejk(1)-ekl(1)*ejk(3) ) / ( rjk*sinphi3*sinphi3)
bmatrix(irowb,jcol+2)= (  rjk - rij * cosphi) * ( eij(1)*ejk(2)-eij(2)*ejk(1) ) / (rjk*rij*sinphi*sinphi) &
+ cosphi3 * (  ekl(1)*ejk(2)-ekl(2)*ejk(1) ) / ( rjk*sinphi3*sinphi3)



bmatrix(irowb,kcol)= (  rjk - rkl * cosphi3) * ( ekl(2)*ejk(3)-ekl(3)*ejk(2) ) / (rjk*rkl*sinphi3*sinphi3) &
+ cosphi * (  eij(2)*ejk(3)-eij(3)*ejk(2) ) / ( rjk*sinphi*sinphi)
bmatrix(irowb,kcol+1)= (  rjk - rkl * cosphi3) * ( ekl(3)*ejk(1)-ekl(1)*ejk(3) ) / (rjk*rkl*sinphi3*sinphi3) &
+ cosphi * (  eij(3)*ejk(1)-eij(1)*ejk(3) ) / ( rjk*sinphi*sinphi)
bmatrix(irowb,kcol+2)= (  rjk - rkl * cosphi3) * ( ekl(1)*ejk(2)-ekl(2)*ejk(1) ) / (rjk*rkl*sinphi3*sinphi3) &
+ cosphi * (  eij(1)*ejk(2)-eij(2)*ejk(1) ) / ( rjk*sinphi*sinphi)



bmatrix(irowb,lcol)=-( ekl(2)*ejk(3)-ekl(3)*ejk(2) )/( rkl * sinphi3*sinphi3)
bmatrix(irowb,lcol+1)=-( ekl(3)*ejk(1)-ekl(1)*ejk(3) )/( rkl * sinphi3*sinphi3)
bmatrix(irowb,lcol+2)=-( ekl(1)*ejk(2)-ekl(2)*ejk(1) )/( rkl * sinphi3*sinphi3)

end do
end if

bmat4int=bmatrix


allocate(indx(numint))
if(allocated(scrmat))deallocate(scrmat)
allocate(scrmat(numint,numint))
call dgemm( 'N', 'T', numint, numint, numcart, one,bmatrix &
 , numint,bmatrix, numint,zero,scrmat,numint )
if(allocated(binv))deallocate(binv)
allocate(binv(numint,numint))
indx=0
call migs(scrmat,numint,binv,indx)
deallocate(scrmat)
deallocate(indx)
allocate(scrmat(numint,numcart))
call dgemm( 'N', 'N', numint, numcart, numint, one,binv &
 , numint,bmatrix, numint,zero,scrmat,numint )
bmatrix=scrmat
deallocate(scrmat)
!call matprt(bmatrix,numint,numcart,numint,numcart)



end subroutine buildb


