subroutine get_basis(x,y,z,numat,basis_set,izcore,ibufflen,debug_ab)
use gaussian_basis
implicit double precision (a-h,o-z)
include 'mpif.h'
character(80) word
character(80),dimension(numat):: basname
double precision,dimension(20,20)::scr_coefs,oscr_coefs
double precision,dimension(20)::scr_exp
double precision,dimension(:)::x,y,z
integer,intent(in)::numat,ibufflen
character(20),intent(in)::basis_set
integer,intent(in),dimension(numat)::izcore
logical,intent(in)::debug_ab
debug_ab2=debug_ab

      call mpi_comm_rank(mpi_comm_world,myrank,ierr)
      call mpi_comm_size(mpi_comm_world,nprocs,ierr)


if(myrank.eq.0)print*,'Two electron integral buffer allocation= ',ibufflen,' MB'
jmaxints=ibufflen*1000/8
if(myrank.eq.0)print*,'Buffer can hold=',1000*jmaxints,'integrals'  
allocate(twoeints(jmaxints*1000))
allocate(bfuncs(4,jmaxints*1000))
allocate(bfuncs_diff(4,1296))
if(myrank.eq.0)print*,'Memory allocation successful. Let''s proceed...'


pi=dacos(-1.0d0)
allocate(coor(3,numat))
allocate(zcore(numat))
natoms=numat
ncmax=20 ! maximum contraction length. number of primitives per function
nmax=10000 ! maximum number of shells in calculation. change if system size is huge
if(myrank.eq.0)print*,''
if(myrank.eq.0)print*,'Basis sets assigned to atomic centers:'
do i=1,numat
k=len_trim(periodic(izcore(i)))
basname(i)(1:k)=periodic(izcore(i))
basname(i)(k+1:k+1)=':'
basname(i)(k+2:)=basis_set
k=len_trim(basname(i))
if(myrank.eq.0)print*,basname(i)(1:k)


coor(1,i)=x(i)
coor(2,i)=y(i)
coor(3,i)=z(i)
zcore(i)=float(izcore(i))

end do
coor=coor/0.52917721067d0


open(unit=130,file='GENBAS')
ishells_total=0

allocate(shell_exp(ncmax,nmax))
allocate(shell_coefs(ncmax,nmax))
allocate(num_contr(nmax))
allocate(iangm(nmax))
allocate(shell_coor(3,nmax))
allocate(index_basis(6,nmax))
allocate(ibfcenter(nmax))

shell_exp=0.
shell_coefs=0.
num_contr=0
iangm=0
shell_coor=0.

do iname=1,numat 
rewind 130

do
 read(130,*,iostat=io)word
 if(io<0)exit
 
 if(word.eq.basname(iname))then
 read(130,*)word
 read(130,*)nsections
 allocate(iang(nsections))
 allocate(icol(nsections))
 allocate(iprim(nsections))
 read(130,*)(iang(i),i=1,nsections)
 read(130,*)(icol(i),i=1,nsections)
 read(130,*)(iprim(i),i=1,nsections)
 scr_exp=0.
 scr_coefs=0.
 ishells_on_atom=0
 do k=1,nsections
  scr_exp=0.
  scr_coefs=0.
  read(130,*)(scr_exp(i),i=1,iprim(k))
   do m=1,iprim(k)
   read(130,*)(scr_coefs(m,n),n=1,icol(k))
   end do
  ishells_on_atom=icol(k)+ishells_on_atom
  ishells_total=ishells_total+icol(k)
  istart=ishells_total-icol(k)+1
  jrow=0
  do j=1,icol(k)
  mrow=0
  igauss=0
  icolumn=istart+j-1
  do m=1,iprim(k)
  if(abs(scr_coefs(m,j)).gt.1d-10)then
    igauss=igauss+1
    mrow=mrow+1
!    jrow=jrow+1
    shell_coefs(mrow,icolumn)=scr_coefs(m,j)   
    shell_exp(mrow,icolumn)=scr_exp(m)
  end if
  end do
   num_contr(icolumn)=igauss
   iangm(icolumn)=iang(k)
   shell_coor(1,icolumn)=x(iname)
   shell_coor(2,icolumn)=y(iname)
   shell_coor(3,icolumn)=z(iname)
   ibfcenter(icolumn)=iname
end do





! accumulate ishells_on_atom array here

end do



exit ! exit the loop over genbas file after finding basis name and  processing
end if


if(allocated(iang))deallocate(iang)
if(allocated(icol))deallocate(icol)
if(allocated(iprim))deallocate(iprim)

end do !end loop over genbas file


end do ! end loop over atoms and basis names


!do i=1,ishells_total
!print*,i,iangm(i),num_contr(i)
!do j=1,20
!print*,shell_coefs(j,i),shell_exp(j,i)
!end do
!end do
!stop


! get total number of basis functions,s-shells,p-shells
allocate(isorb(ishells_total))
allocate(iporb(ishells_total))
allocate(idorb(ishells_total))

num_gauss=0
num_s_shell=0
num_p_shell=0
num_d_shell=0


index_basis=0
icounter=0



do i=1,ishells_total 
if(iangm(i).eq.0)then
num_gauss=num_gauss+1
num_s_shell=num_s_shell+1
isorb(num_s_shell)=i
icounter=icounter+1
index_basis(1,i)=icounter
end if 

if(iangm(i).eq.1)then
num_gauss=num_gauss+2*iangm(i)+1
num_p_shell=num_p_shell+1
iporb(num_p_shell)=i

do kk=1,3
icounter=icounter+1
index_basis(kk,i)=icounter
end do

end if

if(iangm(i).eq.2)then
num_gauss=num_gauss+2*iangm(i)+2
num_d_shell=num_d_shell+1
idorb(num_d_shell)=i
do kk=1,6
icounter=icounter+1
index_basis(kk,i)=icounter
end do

end if

end do






iss_pairs_total=num_s_shell*(num_s_shell+1)/2
allocate(iss_pairs(iss_pairs_total,2))
kcount=0
do i=1,num_s_shell
do j=1,i
kcount=kcount+1
iss_pairs(kcount,1)=isorb(i)
iss_pairs(kcount,2)=isorb(j)
end do
end do




ips_pairs_total=num_p_shell*num_s_shell
allocate(ips_pairs(ips_pairs_total,2))
kcount=0
do i=1,num_p_shell
do j=1,num_s_shell
kcount=kcount+1
ips_pairs(kcount,1)=iporb(i)
ips_pairs(kcount,2)=isorb(j)
end do
end do

ipp_pairs_total=num_p_shell*(num_p_shell+1)/2
allocate(ipp_pairs(ipp_pairs_total,2))
kcount=0
do i=1,num_p_shell
do j=1,i
kcount=kcount+1
ipp_pairs(kcount,1)=iporb(i)
ipp_pairs(kcount,2)=iporb(j)
end do
end do




ids_pairs_total=num_d_shell*num_s_shell
allocate(ids_pairs(ids_pairs_total,2))
kcount=0
do i=1,num_d_shell
do j=1,num_s_shell
kcount=kcount+1
ids_pairs(kcount,1)=idorb(i)
ids_pairs(kcount,2)=isorb(j)
end do
end do





idp_pairs_total=num_d_shell*num_p_shell
allocate(idp_pairs(idp_pairs_total,2))
kcount=0
do i=1,num_d_shell
do j=1,num_p_shell
kcount=kcount+1
idp_pairs(kcount,1)=idorb(i)
idp_pairs(kcount,2)=iporb(j)
end do
end do




idd_pairs_total=num_d_shell*(num_d_shell+1)/2
allocate(idd_pairs(idd_pairs_total,2))
kcount=0
do i=1,num_d_shell
do j=1,i
kcount=kcount+1
idd_pairs(kcount,1)=idorb(i)
idd_pairs(kcount,2)=idorb(j)
end do
end do





shell_coor=shell_coor/.52917721067d0






if(myrank.eq.0)then
print*,''
print*,'Total number of basis functions: ',num_gauss
print*,'Number of s shells:',num_s_shell
print*,'Number of p shells:',num_p_shell
print*,'Number of d shells:',num_d_shell
end if

! normalize primitives
do i=1,num_s_shell
do j=1,num_contr(isorb(i))
shell_coefs(j,isorb(i))=shell_coefs(j,isorb(i))*(2.0D0*shell_exp(j,isorb(i))/pi)**(.75)
end do
end do

do i=1,num_p_shell
do j=1,num_contr(iporb(i))
shell_coefs(j,iporb(i))=shell_coefs(j,iporb(i))*((128.0D0*shell_exp(j,iporb(i))**5)/pi**3)**(.25)
end do
end do

!!!NORMALIZE D ORBS HERE. Still need factor of 9 for xx,yy,zz orbitals
do i=1,num_d_shell
do j=1,num_contr(idorb(i))
shell_coefs(j,idorb(i))=shell_coefs(j,idorb(i))*((2048.0D0*shell_exp(j,idorb(i))**7)/pi**3)**(.25)
end do
end do






iveclength=num_gauss*(num_gauss+1)/2
allocate(ioffset(num_gauss))
do i=1,num_gauss
ioffset(i)=i*(i-1)/2
end do
allocate(ioffset2(iveclength))
do i=1,iveclength
ioffset2(i)=i*(i-1)/2
end do



! compute repulsion energy and electron count

total=0.
nelectrons=0
do i=1,numat-1
do j=i+1,numat
rsq=0
do k=1,3
rsq=rsq+(coor(k,j)-coor(k,i))**2
end do
r=dsqrt(rsq)
total=total+zcore(j)*zcore(i)/r
end do
end do
repulsion_nuc=total
if(myrank.eq.0)print*,'Nuclear repulsion energy: ',repulsion_nuc

do i=1,numat
nelectrons=nelectrons+izcore(i)
end do

if(myrank.eq.0)then
print*,'Total number of electrons: ',nelectrons
open(unit=45,file='charge')
read(45,*)dcharge
read(45,*)lindep_tol
close(45)
end if
nelectrons=nelectrons-dcharge
if(myrank.eq.0)print*,'System charge ',dcharge
if(myrank.eq.0)print*,'Total number of electrons: ',nelectrons



!divvy up the pairs across cores

if(num_s_shell .ne. 0)then
allocate(iss_pairs_cpu(nprocs))
allocate(ifirst_on_cpu_ss(nprocs))
allocate(ilast_on_cpu_ss(nprocs))
allocate(iquads_cpu(nprocs))
allocate(ifirst_quad_cpu(nprocs))
allocate(ilast_quad_cpu(nprocs))
allocate(iall_on_cpu(nprocs))
end if

if(myrank.eq.0)then
a=float(iss_pairs_total)/float(nprocs)
ieach=floor(a)
iss_pairs_cpu=ieach
j=ieach
k=mod(iss_pairs_total,nprocs)
if(k.ne.0)then
j=iss_pairs_total-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
iss_pairs_cpu(icount)=iss_pairs_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
!         print*,'Number of CPUs  =',nprocs
!         print*,'Each cpu will handle',ieach,' ss pairs'
!         if(k.ne.0)print*,'Runaway slave will handle',j, 'ss pairs'
ifirst_on_cpu_ss(1)=1
ilast_on_cpu_ss(1)=ifirst_on_cpu_ss(1)+iss_pairs_cpu(1)-1
do i=2,nprocs
ifirst_on_cpu_ss(i)=ilast_on_cpu_ss(i-1)+1
ilast_on_cpu_ss(i)=ifirst_on_cpu_ss(i)+iss_pairs_cpu(i)-1
end do
if(debug_ab)then
print*,'ss Distribution Among Processors'
print*,'Processor          First ss            Last ss         Total'
do i=1,nprocs
print*,i,'       ',ifirst_on_cpu_ss(i),'       ',ilast_on_cpu_ss(i),ilast_on_cpu_ss(i)-ifirst_on_cpu_ss(i)+1
end do
end if
end if
call mpi_bcast(iss_pairs_cpu,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ifirst_on_cpu_ss,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ilast_on_cpu_ss,nprocs,mpi_integer,0,mpi_comm_world,ierr)




if(num_p_shell .ne. 0)then
allocate(ips_pairs_cpu(nprocs))
allocate(ifirst_on_cpu_ps(nprocs))
allocate(ilast_on_cpu_ps(nprocs))
allocate(ipp_pairs_cpu(nprocs))
allocate(ifirst_on_cpu_pp(nprocs))
allocate(ilast_on_cpu_pp(nprocs))
else
return
end if


if(myrank.eq.0)then
a=float(ips_pairs_total)/float(nprocs)
ieach=floor(a)
ips_pairs_cpu=ieach
j=ieach
k=mod(ips_pairs_total,nprocs)
if(k.ne.0)then
j=ips_pairs_total-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
ips_pairs_cpu(icount)=ips_pairs_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
!         print*,'Number of CPUs  =',nprocs
!         print*,'Each cpu will handle',ieach,' ps pairs'
!         if(k.ne.0)print*,'Runaway slave will handle',j, 'ps pairs'
ifirst_on_cpu_ps(1)=1
ilast_on_cpu_ps(1)=ifirst_on_cpu_ps(1)+ips_pairs_cpu(1)-1
do i=2,nprocs
ifirst_on_cpu_ps(i)=ilast_on_cpu_ps(i-1)+1
ilast_on_cpu_ps(i)=ifirst_on_cpu_ps(i)+ips_pairs_cpu(i)-1
end do
if(debug_ab)then
print*,'ps Distribution Among Processors'
print*,'Processor          First ps            Last ps         Total'
do i=1,nprocs
print*,i,'       ',ifirst_on_cpu_ps(i),'       ',ilast_on_cpu_ps(i),ilast_on_cpu_ps(i)-ifirst_on_cpu_ps(i)+1
end do
end if
end if
call mpi_bcast(ips_pairs_cpu,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ifirst_on_cpu_ps,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ilast_on_cpu_ps,nprocs,mpi_integer,0,mpi_comm_world,ierr)





if(myrank.eq.0)then
a=float(ipp_pairs_total)/float(nprocs)
ieach=floor(a)
ipp_pairs_cpu=ieach
j=ieach
k=mod(ipp_pairs_total,nprocs)
if(k.ne.0)then
j=ipp_pairs_total-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
ipp_pairs_cpu(icount)=ipp_pairs_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
!         print*,'Number of CPUs  =',nprocs
!         print*,'Each cpu will handle',ieach,' pp pairs'
!         if(k.ne.0)print*,'Runaway slave will handle',j, 'pp pairs'
ifirst_on_cpu_pp(1)=1
ilast_on_cpu_pp(1)=ifirst_on_cpu_pp(1)+ipp_pairs_cpu(1)-1
do i=2,nprocs
ifirst_on_cpu_pp(i)=ilast_on_cpu_pp(i-1)+1
ilast_on_cpu_pp(i)=ifirst_on_cpu_pp(i)+ipp_pairs_cpu(i)-1
end do
if(debug_ab)then
print*,'pp Distribution Among Processors'
print*,'Processor          First pp            Last pp         Total'
do i=1,nprocs
print*,i,'       ',ifirst_on_cpu_pp(i),'       ',ilast_on_cpu_pp(i),ilast_on_cpu_pp(i)-ifirst_on_cpu_pp(i)+1
end do
end if
end if
call mpi_bcast(ipp_pairs_cpu,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ifirst_on_cpu_pp,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ilast_on_cpu_pp,nprocs,mpi_integer,0,mpi_comm_world,ierr)




call mpi_barrier(mpi_comm_world,ierr)

if(num_d_shell .ne. 0)then
allocate(ids_pairs_cpu(nprocs))
allocate(ifirst_on_cpu_ds(nprocs))
allocate(ilast_on_cpu_ds(nprocs))
allocate(idp_pairs_cpu(nprocs))
allocate(ifirst_on_cpu_dp(nprocs))
allocate(ilast_on_cpu_dp(nprocs))
allocate(idd_pairs_cpu(nprocs))
allocate(ifirst_on_cpu_dd(nprocs))
allocate(ilast_on_cpu_dd(nprocs))
else
return
end if



if(myrank.eq.0)then
a=float(ids_pairs_total)/float(nprocs)
ieach=floor(a)
ids_pairs_cpu=ieach
j=ieach
k=mod(ids_pairs_total,nprocs)
if(k.ne.0)then
j=ids_pairs_total-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
ids_pairs_cpu(icount)=ids_pairs_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
!         print*,'Number of CPUs  =',nprocs
!         print*,'Each cpu will handle',ieach,' ds pairs'
!         if(k.ne.0)print*,'Runaway slave will handle',j, 'ds pairs'
ifirst_on_cpu_ds(1)=1
ilast_on_cpu_ds(1)=ifirst_on_cpu_ds(1)+ids_pairs_cpu(1)-1
do i=2,nprocs
ifirst_on_cpu_ds(i)=ilast_on_cpu_ds(i-1)+1
ilast_on_cpu_ds(i)=ifirst_on_cpu_ds(i)+ids_pairs_cpu(i)-1
end do
if(debug_ab)then
print*,'ds Distribution Among Processors'
print*,'Processor          First ds            Last ds         Total'
do i=1,nprocs
print*,i,'       ',ifirst_on_cpu_ds(i),'       ',ilast_on_cpu_ds(i),ilast_on_cpu_ds(i)-ifirst_on_cpu_ds(i)+1
end do
end if
end if
call mpi_bcast(ids_pairs_cpu,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ifirst_on_cpu_ds,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ilast_on_cpu_ds,nprocs,mpi_integer,0,mpi_comm_world,ierr)


if(myrank.eq.0)then
a=float(idp_pairs_total)/float(nprocs)
ieach=floor(a)
idp_pairs_cpu=ieach
j=ieach
k=mod(idp_pairs_total,nprocs)
if(k.ne.0)then
j=idp_pairs_total-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
idp_pairs_cpu(icount)=idp_pairs_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
!         print*,'Number of CPUs  =',nprocs
!         print*,'Each cpu will handle',ieach,' dp pairs'
!         if(k.ne.0)print*,'Runaway slave will handle',j, 'dp pairs'
ifirst_on_cpu_dp(1)=1
ilast_on_cpu_dp(1)=ifirst_on_cpu_dp(1)+idp_pairs_cpu(1)-1
do i=2,nprocs
ifirst_on_cpu_dp(i)=ilast_on_cpu_dp(i-1)+1
ilast_on_cpu_dp(i)=ifirst_on_cpu_dp(i)+idp_pairs_cpu(i)-1
end do
if(debug_ab)then
print*,'dp Distribution Among Processors'
print*,'Processor          First dp            Last dp         Total'
do i=1,nprocs
print*,i,'       ',ifirst_on_cpu_dp(i),'       ',ilast_on_cpu_dp(i),ilast_on_cpu_dp(i)-ifirst_on_cpu_dp(i)+1
end do
end if
end if
call mpi_bcast(idp_pairs_cpu,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ifirst_on_cpu_dp,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ilast_on_cpu_dp,nprocs,mpi_integer,0,mpi_comm_world,ierr)


if(myrank.eq.0)then
a=float(idd_pairs_total)/float(nprocs)
ieach=floor(a)
idd_pairs_cpu=ieach
j=ieach
k=mod(idd_pairs_total,nprocs)
if(k.ne.0)then
j=idd_pairs_total-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
idd_pairs_cpu(icount)=idd_pairs_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
!         print*,'Number of CPUs  =',nprocs
!         print*,'Each cpu will handle',ieach,' dd pairs'
!         if(k.ne.0)print*,'Runaway slave will handle',j, 'dd pairs'
ifirst_on_cpu_dd(1)=1
ilast_on_cpu_dd(1)=ifirst_on_cpu_dd(1)+idd_pairs_cpu(1)-1
do i=2,nprocs
ifirst_on_cpu_dd(i)=ilast_on_cpu_dd(i-1)+1
ilast_on_cpu_dd(i)=ifirst_on_cpu_dd(i)+idd_pairs_cpu(i)-1
end do
if(debug_ab)then
print*,'dd Distribution Among Processors'
print*,'Processor          First dd            Last dd         Total'
do i=1,nprocs
print*,i,'       ',ifirst_on_cpu_dd(i),'       ',ilast_on_cpu_dd(i),ilast_on_cpu_dd(i)-ifirst_on_cpu_dd(i)+1
end do
end if
end if
call mpi_bcast(idd_pairs_cpu,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ifirst_on_cpu_dd,nprocs,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ilast_on_cpu_dd,nprocs,mpi_integer,0,mpi_comm_world,ierr)





call mpi_barrier(mpi_comm_world,ierr)
















return
end
