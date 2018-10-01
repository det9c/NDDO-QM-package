! once center part
subroutine build_g_se(density,gmat,gtmp,gmat2,gtmp2) !(density,g)
use tables
use indices
use control
use constants
implicit double precision (a-h,o-z)
interface

      subroutine pack1(indi,indj,ij)
integer,intent(in)::indi,indj
integer,intent(inout)::ij
end subroutine pack1
end interface
double precision,dimension(:,:),intent(in)::density
double precision,dimension(:,:),intent(inout)::gmat,gtmp,gmat2,gtmp2
double precision,dimension(9,9)::dmat1,dmat2
integer::dens_row1,dens_row2
include 'mpif.h'


!half=.5d0



gmat=0.
gmat2=0.
icounter=0
ishift2=ifirstbf(ifirst_atom_on_cpu(myrank+1))-1

do i=ifirst_atom_on_cpu(myrank+1),ilast_atom_on_cpu(myrank+1)
gtmp=0.
gtmp2=0.
do kk=1,num_neighbors(i)
jspot=ifirst_neighbor(i)+kk-1
j=neighbors(2,jspot)
icounter=icounter+1


dens_row1=denslist_map1(icounter)-1
nj=nbas(species(j))
ll=0
do i1=1,nj
do i2=1,nj
ll=ll+1
dmat1(i2,i1)=pair_dens1(dens_row1+ll)
!if(myrank.eq.0)print*,dens_row1,i,j
end do
end do
!if(myrank.eq.0)print*,'dma1',i,j,dmat1

dens_row2=denslist_map2(icounter)-1
ni=nbas(species(i))
nj=nbas(species(j))
ll=0
do i1=1,nj
do i2=1,ni
ll=ll+1
dmat2(i2,i1)=pair_dens2(dens_row2+ll)
end do
end do
!if(myrank.eq.0)print*,'dma2',i,j,dmat2
ishift=ifirstbf(i)-1
jshift=ifirstbf(j)-1


do irow=ifirst(icounter),ilast(icounter)
!print*,irow,ifirst(i),ilast(i),i
imu=map_pairs(1,tpair(1,irow))
inu=map_pairs(2,tpair(1,irow))
ila=map_pairs(1,tpair(2,irow))
isig=map_pairs(2,tpair(2,irow))
value=twoe(irow)



!!!if(io<0)exit
!coulomb contribution
double=1.0d0
scale1=1.
scale2=1.
call pack1(imu,inu,imu_nu)
call pack1(ila,isig,ila_sig)
if(imu.ne.inu)scale1=2.0d0
if(ila.ne.isig)scale2=2.0d0
if(imu_nu .eq. ila_sig)double=.5
gtmp(imu,inu-ishift2)=gtmp(imu,inu-ishift2)+dmat1(ila-jshift,isig-jshift)*scale2*value*double
gtmp(inu,imu-ishift2)=gtmp(imu,inu-ishift2)

if(i.eq.j)then
gtmp(ila,isig-ishift2)=gtmp(ila,isig-ishift2)+dmat1(imu-ishift,inu-ishift)*scale1*value*double
gtmp(isig,ila-ishift2)=gtmp(ila,isig-ishift2)
end if

!exchange contribution
call pack1(imu,ila,imu_la)
call pack1(imu,isig,imu_sig)
call pack1(inu,ila,inu_la)
call pack1(inu,isig,inu_sig)
factor=1.0
if(ila .eq. isig)factor=half
if(imu .eq. inu )factor=half*factor
scale3=1.
if(imu .eq. ila .and. inu.ne.isig)scale3=2.
gtmp(ila,imu-ishift2)=gtmp(ila,imu-ishift2)- half*dmat2(inu-ishift,isig-jshift)*value*scale3*factor
if(i.eq.j)gtmp(imu,ila-ishift2)=gtmp(ila,imu-ishift2)
scale3=1.
factor2=1.0d0
if(imu_nu .eq. ila_sig .and. imu .ne. inu)factor2=.5d0 !(uv|uv) ->uu uv vu vv so need half since uv vu is twice
if(imu .eq. isig .and.  inu.ne.ila)scale3=2.
gtmp(isig,imu-ishift2)=gtmp(isig,imu-ishift2)- half*dmat2(inu-ishift,ila-jshift)*value*scale3*factor*factor2
if(i.eq.j)gtmp(imu,isig-ishift2)=gtmp(isig,imu-ishift2)
scale3=1.
if(inu .eq. ila .and.  imu.ne.isig)scale3=2.
gtmp(ila,inu-ishift2)=gtmp(ila,inu-ishift2)- half*dmat2(imu-ishift,isig-jshift)*value*scale3*factor*factor2
if(i.eq.j)gtmp(inu,ila-ishift2)=gtmp(ila,inu-ishift2)
scale3=1.
if(inu.eq.isig  .and.  imu.ne.ila)scale3=2.
gtmp(isig,inu-ishift2)=gtmp(isig,inu-ishift2)- half*dmat2(imu-ishift,ila-jshift)*value*scale3*factor
if(i.eq.j)gtmp(inu,isig-ishift2)=gtmp(isig,inu-ishift2)


end do
!print*,'i',i,ifirstbf(i),ilastbf(i)


end do ! end loop over j (or kk as indexed)




gtmp2=gtmp



!if(myrank.eq.0)then
!call matprt(gtmp,num_basis,icols_cpu,num_basis,icols_cpu)
!stop
!end if


!h=0.

do ii=ifirsthc2c(i),ilasthc2c(i)

gtmp(hpair(2,ii),hpair(1,ii)-ishift2)=gtmp(hpair(2,ii),hpair(1,ii)-ishift2)+h(ii)
!gtmp(hpair(1,ii),hpair(2,ii))=gtmp(hpair(2,ii),hpair(1,ii))
gtmp2(hpair(2,ii),hpair(1,ii)-ishift2)=gtmp(hpair(2,ii),hpair(1,ii)-ishift2)+h(ii)
!gtmp2(hpair(1,ii),hpair(2,ii))=gtmp2(hpair(2,ii),hpair(1,ii))
end do

do ii=ifirsthc1c(i),ilasthc1c(i)
gtmp(hpair(1,ii),hpair(2,ii)-ishift2)=gtmp(hpair(1,ii),hpair(2,ii)-ishift2)+h(ii)
gtmp(hpair(2,ii),hpair(1,ii)-ishift2)=gtmp(hpair(1,ii),hpair(2,ii)-ishift2)
gtmp2(hpair(1,ii),hpair(2,ii)-ishift2)=gtmp(hpair(1,ii),hpair(2,ii)-ishift2)+h(ii)
gtmp2(hpair(2,ii),hpair(1,ii)-ishift2)=gtmp2(hpair(1,ii),hpair(2,ii)-ishift2)
end do





gmat(1:num_basis,ifirstbf(i)-ishift2:ilastbf(i)-ishift2)=gtmp(1:num_basis,ifirstbf(i)-ishift2:ilastbf(i)-ishift2)
gmat2(1:num_basis,ifirstbf(i)-ishift2:ilastbf(i)-ishift2)=gtmp2(1:num_basis,ifirstbf(i)-ishift2:ilastbf(i)-ishift2)



!do jj=ifirstbf(i),ilastbf(i)
!do ii=ifirstbf(i),num_basis
!print*,ii,jj
!gmat(ii,jj)=gtmp(ii,jj)
!end do
!end do


!if(i.eq.1)return
end do ! end loop over i






!call mpi_reduce(g,glocal,iveclength,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
!g=glocal

!call mpi_barrier(mpi_comm_world,ierr)

return
end
