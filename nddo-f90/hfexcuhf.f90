subroutine hfexcuhf(acoulomb,bcoulomb,abcoulomb,ahfx,bhfx,ekinetic)
use scratch_array
use constants
use tables
use indices

implicit double precision (a-h,o-z)


interface
subroutine trace(S,N,out)
double precision,dimension(:,:),intent(in)::S
integer::N
double precision,intent(out)::out
end subroutine trace
end interface

double precision,intent(inout)::acoulomb,bcoulomb,abcoulomb,ahfx,bhfx,ekinetic
double precision,dimension(:,:),allocatable::apa,bpa,apb,bpb,apab,bpab
integer::si
acoulomb=zero
ahfx=zero
bcoulomb=zero
bhfx=zero
abcoulomb=zero
do i=1,numat
ni=nbas(species(i))
istart=ifirst(i)
iend=ilast(i)
istart2=ifirst2(i)
iend2=ilast2(i)
allocate(apa(ni,ni))
allocate(bpa(ni,ni))

do j=1,numat
nj=nbas(species(j))
jstart=ifirst(j)
jend=ilast(j)
jstart2=ifirst2(j)
jend2=ilast2(j)
allocate(apb(nj,nj))
allocate(bpb(nj,nj))
allocate(apab(ni,nj))
allocate(bpab(ni,nj))

apa(1:ni,1:ni)=adens(istart:iend,istart:iend)
apb(1:nj,1:nj)=adens(jstart:jend,jstart:jend)
apab(1:ni,1:nj)=adens(istart:iend,jstart:jend)

bpa(1:ni,1:ni)=bdens(istart:iend,istart:iend)
bpb(1:nj,1:nj)=bdens(jstart:jend,jstart:jend)
bpab(1:ni,1:nj)=bdens(istart:iend,jstart:jend)


do mu=1,ni
do nu=1,ni
 do la=1,nj
 do si=1,nj
value=twoe(jstart2+map(si,la),istart2+map(nu,mu))

acoulomb=acoulomb+apa(nu,mu)*apb(si,la)*value
bcoulomb=bcoulomb+bpa(nu,mu)*bpb(si,la)*value
abcoulomb=abcoulomb+apa(nu,mu)*bpb(si,la)*value
ahfx=ahfx+apab(mu,la)*apab(nu,si)*value
bhfx=bhfx+bpab(mu,la)*bpab(nu,si)*value
end do
end do
end do
end do

deallocate(apb)
deallocate(bpb)
deallocate(apab)
deallocate(bpab)

end do
deallocate(apa)
deallocate(bpa)
end do

allocate(scrmat(num_basis,num_basis))
call dgemm( 'N', 'N', num_basis, num_basis,num_basis, one,s3 &
 , num_basis,H , num_basis,zero,scrmat,num_basis )
call trace(scrmat,num_basis,ekinetic)

acoulomb=acoulomb/two
bcoulomb=bcoulomb/two
ahfx=-ahfx/two
bhfx=-bhfx/two
deallocate(scrmat)
end subroutine hfexcuhf
