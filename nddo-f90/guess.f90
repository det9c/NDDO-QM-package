subroutine guess(density)
use constants
use tables
use indices
use control
implicit double precision(a-h,o-z)
include 'mpif.h'
!pdouble precision,dimension(:,:),intent(inout)::density
double precision,dimension(:,:),intent(inout)::density

ishift=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1
density=zero
do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)

ni=nbas(species(i))

   if(ni>4.and.nqp(species(i))<4)then
          term=eff_core(zeff(i))/four
         imax=4
    else
         term=eff_core(zeff(i))/float(ni)
         imax=nbas(species(i))
   end if
! correct for hydrogen with an sp basis set b/c we don't want to populate p orbitals
if(eff_core(zeff(i))==1.and.ni>1)then

term=one
imax=1
end if


m=ifirstbf(i)
!n=ilast(i)
n=m+imax-1
    do j=m,n
    density(j,j-ishift)=term
!       density(j+offset1(j))=term
    end do
 
 end do
end subroutine guess
