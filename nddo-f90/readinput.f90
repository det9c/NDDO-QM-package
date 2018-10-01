!!
! parallelize neighborlist and spread out  over cores
! spread overlap and hcore out over cores /done
!  speed up hcore subroutine with neighbor list in Zb<us|ss> loop
!  remove the allocate statment from  Zb<us|ss> loop in hcore
! add code to check 2ebuffer overrun in twoe_electron.f90 (itotal>max length) 
! don't hold all one center integrals in core on first process. so parallelize bottom part of twoe_electron.f90 
!  fix build_g_se routine to not use gtmp2. can't make it work. don't know why
  subroutine readin
 use constants
 use tables
 use control
 use scratch_array
 use indices
 implicit double precision (a-h,o-z)
  include 'mpif.h'
!****************************
! include headers for all subroutines called by this program. this is painful

interface

subroutine get_basis(x,y,z,numat,basis,zcore,ERIBUFF,debug)
double precision,dimension(:)::x,y,z
integer,intent(in)::numat,ERIBUFF
integer,intent(in),dimension(numat)::zcore
character(20),intent(in)::basis
logical,intent(in)::debug
end subroutine get_basis



subroutine square(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(inout)::matrix
double precision,dimension(:),intent(in)::vector
integer,intent(in)::idim1
end subroutine square

subroutine hfdft
end subroutine hfdft

subroutine geometry
end subroutine geometry


     subroutine two_electron
     end subroutine two_electron

subroutine hcore
end subroutine hcore

subroutine loadparam
end subroutine loadparam

subroutine rhf
end subroutine rhf


subroutine default(ikey)
integer,intent(in)::ikey
end subroutine default


subroutine pm3_param(k,ussa,uppa,betasa,betapa,zsa,zpa,gssa,gspa,gppa,    &
 gp2a,hspa,d1a,d2a,nbasis,alp,a11,a21,a31,a12,a22,a32,a13,a23,a33)
 double precision,intent(inout):: ussa,uppa,betasa,betapa,zsa,zpa  &
 ,gssa,gspa,gppa,gp2a,hspa,d1a,  &
 d2a,alp,a11,a21,a31,a12,a22,a32,a13,a23,a33
 integer,intent(inout)::nbasis
end subroutine pm3_param

 subroutine am1_param(k,ussa,uppa,betasa,betapa,zsa,zpa,gssa,gspa,gppa,    &
 gp2a,hspa,d1a,d2a,nbasis,alp,a11,a21,a31,a12,a22,a32,a13,a23,a33,a14,a24,a34)
 double precision,intent(inout):: ussa,uppa,betasa,betapa,zsa,zpa  &
 ,gssa,gspa,gppa,gp2a,hspa,d1a,  &
 d2a,alp,a11,a21,a31,a12,a22,a32,a13,a23,a33,a14,a24,a34
 integer,intent(inout)::nbasis
end subroutine am1_param

 subroutine mndod_param(k,ussa,uppa,udda,betasa,betapa,betada,zsa,zpa,zda,gssa,gspa,gppa,    &
 gp2a,hspa,nbasis,alp,gdd,zsn,zpn,zdn,a11,a21,a31,a12,a22,a32,a13,a23,a33,a14,a24,a34)
 implicit double precision (a-h,o-z)
 double precision,intent(inout):: ussa,uppa,udda,betasa,betapa,betada,zsa,zpa,zda  &
 ,gssa,gspa,gppa,gp2a,hspa,alp,gdd,zsn,zpn,zdn,a11,a21,a31,a12,a22,a32,a13,a23,a33,a14,a24,a34
 integer,intent(inout)::nbasis
end subroutine mndod_param

subroutine find
end subroutine find

subroutine ptchg
end subroutine ptchg

subroutine parout
end subroutine parout

SUBROUTINE INIGHD (NI,ns,np,nd,es,ep,ed)
integer,intent(in)::NI,ns,np,nd
double precision,intent(in)::es,ep,ed
end subroutine inighd

SUBROUTINE WSTORE (NI,MODE,ip)
integer,intent(in)::mode,ip,NI
end subroutine wstore

subroutine om1(k,as,ap,ad,bs,bp,bd)
implicit double precision (a-h,o-z)
double precision,intent(inout):: as,ap,ad,bs,bp,bd
integer,intent(in):: k
end subroutine om1

subroutine inertia(reori)
logical,intent(in)::reori
end subroutine inertia

subroutine atoms_cpu(natoms)
integer,intent(in)::natoms
end subroutine atoms_cpu


end interface
!******************************************

character(80) title ,word
character(10) date,time,pretty_time
character(20)ch
double precision ,dimension(:), allocatable ::oldx,oldy,oldz
double precision,dimension(:,:),allocatable::oldq
integer,dimension(:,:),allocatable::oldref,oldopt
integer ,dimension(:), allocatable ::oldzeff,oldntype,num_type
integer ,dimension(:), allocatable ::oldoptx,oldopty,oldoptz,optx,opty,optz
double precision::a2,a4,a6
integer::a1,a3,a5,a7,a8,a9,a10
logical::unique





allocate(neighbors(2,5000000))
allocate(iatom_pairs_cpu(nprocs))
allocate(ifirst_pair_on_cpu(nprocs))
allocate(ilast_pair_on_cpu(nprocs))
allocate(ifirst_atom_on_cpu(nprocs))
allocate(ilast_atom_on_cpu(nprocs))
allocate(atoms_on_cpu(nprocs))
allocate(list_length(nprocs))
allocate(ifirst_occ_on_cpu(nprocs))
allocate(ilast_occ_on_cpu(nprocs))
allocate(occupied_on_cpu(nprocs))
list_length=0



 call date_and_time(date,time)

if(master) write(*,*)'Calculation started at: ',time(1:2),':',time(3:4),':',time(5:6) & 
 ,' on ',date(5:6) ,'/', date(7:8),'/',date(1:4)


open(unit=20,file='INPUT')
read(20,*)title

if(master)write(*,*)'Job title :',title
if(master)write(*,*)''
! read the keyword section

read(20,*)word
k=len_trim(word)
ikey=1
do i=1,k
ch=word(i:i)
if(ch==' ')then
ikey=ikey+1
end if
end do
if(master)print*,'Found ',ikey,' keywords defined by user.  Defaults chosen for remaining variables.'
ibegin=1
allocate(keyword(ikey))
do i=1,ikey
  icount=0
  ch=' '
  do j=ibegin,k
     if(word(j:j)==' ')then
     ibegin=j+1
     exit
     else
     icount=icount+1
     ch(icount:icount)=word(j:j)
     end if
  end do
 keyword(i)=ch

end do

! now, set the defaults

call default(ikey)


! read the atomic numbers and coordinates. the use of the allocatable statement
! here is a fancy schmancy way of using fortran 90 to read a file in which
! the length is not yet known (i.e. a geometry input section)

numat=0


if(XYZ)then
allocate( x(1))
allocate( optx(1))
allocate( y(1))
allocate( opty(1))
allocate( z(1))
allocate( optz(1))
allocate( zeff(1))
if(master)WRITE(*,*)'Cartesian coordinate system will be used.'
 do
! read(20,*,iostat=io)a,bb,i1,c,i2,d,i3
read(20,*,iostat=io)a,bb,c,d
 if(io<0)exit
!bb=bb*0.52917721067d0
!c=c*0.52917721067d0
!d=d*0.52917721067d0




 numat=numat+1
allocate(oldx(numat))
allocate(oldy(numat))
allocate(oldz(numat))
allocate(oldzeff(numat))
allocate(oldoptx(numat))
allocate(oldopty(numat))
allocate(oldoptz(numat))
 oldx=x
 oldy=y
 oldz=z
 oldzeff=zeff
 oldoptx=optx
 oldopty=opty
 oldoptz=optz
 deallocate(x)
 deallocate(y)
 deallocate(z)
 deallocate(zeff)
 deallocate(optx)
 deallocate(opty)
 deallocate(optz)
 allocate( x(numat))
 allocate( y(numat))
 allocate( z(numat))
 allocate( zeff(numat))
 allocate( optx(numat))
 allocate( opty(numat))
 allocate( optz(numat))
 x = oldx
 y = oldy
 z = oldz
 zeff = oldzeff
 optx=oldoptx
 opty=oldopty
 optz=oldoptz
 zeff(numat)=a
 x(numat)=bb
 y(numat)=c
 z(numat)=d
 optx(numat)=i1
 opty(numat)=i2
 optz(numat)=i3
 deallocate(oldx)
 deallocate(oldy)
 deallocate(oldz)
 deallocate(oldzeff)
 deallocate(oldoptx)
 deallocate(oldopty)
 deallocate(oldoptz)
 end do
ninput=numat
! if not xyz input, read the zmat
else
if(master)write(*,*)'Geometry input with internal coordinates'
numat=1
allocate( q(1,3))
allocate( opt(1,3))
allocate( ref(1,3))
allocate( zeff(1))
read(20,*,iostat=io)a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
 zeff(numat)=a1
 q(numat,1)=a2
q(numat,2)=a4
q(numat,3)=a6
opt(numat,1)=a3
opt(numat,2)=a5
opt(numat,3)=a7
ref(numat,1)=a8
ref(numat,2)=a9
ref(numat,3)=a10
 do

 read(20,*,iostat=io)a1,a2,a3,a4,a5,a6,a7,a8,a9,a10

 if(io<0)exit
 numat=numat+1
allocate(oldq(numat,3))
allocate(oldopt(numat,3))
allocate(oldref(numat,3))
allocate(oldzeff(numat))
 oldq(1:numat-1,1:3)=q(1:numat-1,1:3)
 oldopt(1:numat-1,1:3)=opt(1:numat-1,1:3)
 oldref(1:numat-1,1:3)=ref(1:numat-1,1:3)
 oldzeff=zeff
 deallocate(q)
 deallocate(opt)
 deallocate(ref)
 deallocate(zeff)
 allocate( q(numat,3))
 allocate( opt(numat,3))
 allocate( ref(numat,3))
 allocate( zeff(numat))
 q(1:numat-1,1:3) = oldq(1:numat-1,1:3)
 opt(1:numat-1,1:3) = oldopt(1:numat-1,1:3)
 ref(1:numat-1,1:3) = oldref(1:numat-1,1:3)
 zeff = oldzeff

 zeff(numat)=a1
 q(numat,1)=a2
q(numat,2)=a4
q(numat,3)=a6
opt(numat,1)=a3
opt(numat,2)=a5
opt(numat,3)=a7
ref(numat,1)=a8
ref(numat,2)=a9
ref(numat,3)=a10
 deallocate(oldq)
 deallocate(oldopt)
 deallocate(oldref)
 deallocate(oldzeff)
 end do



 end if
close(1)
close(20)
if(.not.XYZ)then


allocate(x(numat))
allocate(y(numat))
allocate(z(numat))
! make sure user has input proper connectivity of first three atoms
if(numat>2.and.ref(3,1)/=2)then
if(master)print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
if(master)print*,'NOTICE: Ordering of first two atoms in Z matrix has been changed by program'
if(master)print*,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
a1=zeff(2)
zeff(2)=zeff(1)
zeff(1)=a1
do j=1,3
do i=3,numat
if(ref(i,j)==1)then
ref(i,j)=-1
end if
end do
end do








do j=1,3
do i=3,numat
if(ref(i,j)==2)then
ref(i,j)=1
end if
end do
end do

do j=1,3
do i=3,numat

if(ref(i,j)==-1)then
ref(i,j)=2
end if
end do
end do
end if
factor=pi/180.0d0
do i=1,numat
q(i,2)=factor*q(i,2)
q(i,3)=factor*q(i,3)
end do
ninput=numat





!make array of internal coordinate addresses
numint=3*numat-6
allocate(intadd(numint,2))
allocate(inttype(numint))
icount=0
do i=2,ninput
icount=icount+1
intadd(icount,1)=i
intadd(icount,2)=1
inttype(icount)='Bond '
end do
do i=3,ninput
icount=icount+1
intadd(icount,1)=i
intadd(icount,2)=2
inttype(icount)='Angle'
end do
do i=4,ninput
icount=icount+1
intadd(icount,1)=i
intadd(icount,2)=3
inttype(icount)='Dihed'
end do
! end of internal coordinate addresses in q array








call geometry
end if

if(master)write(*,*)'Found ',numat,'entries in the geometry section'



! print input coordinates
if(master)write(*,*)'           CARTESIAN ATOMIC COORDINATES ON INPUT: (ANGSTROMS):'
if(master)write(*,*)'           ---------------------------------------------------'
if(master)write(*,*)'ATOM','       ATOMIC NO.        ','    X','                   Y','                  Z'

            imax=numat  
            if(numat.gt.100)imax=100

                do i=1,imax
                if(master) WRITE(*,20)i,ZEFF(I),X(I),Y(I),Z(I) 
                end do
                20 format(i5,i15,f20.5,f20.5,f20.5)
                

if(master)then
                do j=1,numat
!                 WRITE(*,20)j,ZEFF(j),X(j),Y(j),Z(j)
                 WRITE(7,71)periodic(zeff(j)),X(j),Y(j),Z(j)
                end do
                71 format(A,f20.5,f20.5,f20.5)
end if


                
if(master)print*,'End Geometry-------------------------------------------------'

if(pbc)then
if(master)print*,'---------------------------------------------'
if(master)PRINT*,' periodic boundary conditions will be used'
if(master)print*,'---------------------------------------------'
open(unit=1,file='cell')
read(1,*)cell
if(master)call matprt(cell,3,3,3,3)
allocate(fracs(3,numat))
close(1)
end if





! save original atom labels and remove dummy/effective atoms from zmatrix
allocate(zstore(ninput)) 
zstore=zeff
if(.not.XYZ)call buildb
call nodummy
!call inertia(.false.)

if(.not.XYZ .and. ninput.eq.2)then
! for some reason the three lines below this block-if will not execute when the
! geometry is in internals, even though it is the same execution
! statement for xyz input format.  i have no idea what the hell 
! is going on so make the user use XYZ for diatomics until i fix it.
if(master)write(*,*)'*******************************************************'
if(master)write(*,*)'*******************************************************'
if(master)write(*,*)'*******************************************************'
if(master)write(*,*)
if(master)write(*,*)'Due to an unknown program bug which for the life of me I can not'
if(master)write(*,*)'fix, please input diatomics in XYZ coordinates and not internals!'
if(master)write(*,*)
if(master)write(*,*)'*******************************************************'
if(master)write(*,*)'*******************************************************'
if(master)write(*,*)'*******************************************************'
stop
end if
! determine number of unique atoms and assign types to the atoms in zeff




itype=1
allocate(ntype(1)) ! does not execute in internals ???????????????????????
ntype(1)=zeff(1)


do i=2,numat
   j=zeff(i)
   unique=.true.

   do k=1,itype
   if(j.eq.ntype(k))then
   unique=.false.
   exit
   end if
   end do
if(unique)then
allocate(oldntype(itype))
oldntype=ntype
itype=itype+1

deallocate(ntype)
allocate(ntype(itype))
ntype=oldntype
ntype(itype)=j
deallocate(oldntype)
end if
end do

if(.not. abinitio)then !start nddo block
! allocate arrays for parameters 
allocate(uss(itype))
allocate(upp(itype))
allocate(betas(itype))
allocate(betap(itype))
allocate(zs(itype))
allocate(zp(itype))
allocate(gss(itype))
allocate(gsp(itype))
allocate(gpp(itype))
allocate(gp2(itype))
allocate(hsp(itype))
allocate(d1(itype))
allocate(d2(itype))
allocate(nbas(itype))
allocate(p0(itype))
allocate(p1(itype))
allocate(p2(itype))
allocate(alpha(itype))
allocate(g1(4,itype))
allocate(g2(4,itype))
allocate(g3(4,itype))
g1=zero
g2=zero
g3=zero
allocate(hyfsp(itype))
allocate(alpha3(itype))
allocate(ac(2,itype))
allocate(bc(2,itype))


if(METHOD=='MNDOD')then
allocate(udd(itype))
allocate(betad(itype))
allocate(gdd(itype))
allocate(p3(itype))
allocate(p4(itype))
allocate(p5(itype))
allocate(p6(itype))
allocate(p7(itype))
allocate(d3(itype))
allocate(d4(itype))
allocate(d5(itype))
allocate(d0(itype))
d0=zero
allocate(zetad(itype))
allocate(F0DD(itype))
allocate(F2DD(itype))
allocate(F4DD(itype))
allocate(F0SD(itype))
allocate(G2SD(itype))
allocate(F0PD(itype))
allocate(F2PD(itype))
allocate(G1PD(itype))
allocate(G3PD(itype))
allocate(IF0SD(itype))
allocate(IG2SD(itype))
call fbinom
allocate(repd(52,itype))
allocate(hpp(itype))
allocate(zsone(itype))
allocate(zpone(itype))
allocate(zdone(itype))
allocate(hyfpd(itype))
end if
end if ! end nddo block 



allocate(species(numat))
do i=1,numat
k=zeff(i)
do j=1,itype
if(k==ntype(j))exit
end do
species(i)=j
end do

if(.not. abinitio)then
! get the parameters
do i=1,itype
if(method=='PM3')then
call pm3_param(ntype(i),uss(i),upp(i),betas(i),betap(i),zs(i),zp(i),gss(i),gsp(i),gpp(i),    &
 gp2(i),hsp(i),d1(i),d2(i),nbas(i),alpha(i),g1(1,i),g2(1,i),g3(1,i),g1(2,i),g2(2,i),g3(2,i) &
,g1(3,i),g2(3,i),g3(3,i))

elseif(method=='AM1')then
call am1_param(ntype(i),uss(i),upp(i),betas(i),betap(i),zs(i),zp(i),gss(i),gsp(i),gpp(i),    &
 gp2(i),hsp(i),d1(i),d2(i),nbas(i),alpha(i),g1(1,i),g2(1,i),g3(1,i),g1(2,i),g2(2,i),g3(2,i) &
,g1(3,i),g2(3,i),g3(3,i),g1(4,i),g2(4,i),g3(4,i))

elseif(method=='MNDOD')then
call mndod_param(ntype(i),uss(i),upp(i),udd(i),betas(i),betap(i),betad(i),zs(i),zp(i),zetad(i),gss(i),gsp(i),gpp(i),    &
 gp2(i),hsp(i),nbas(i),alpha(i),gdd(i),zsone(i),zpone(i),zdone(i),g1(1,i),g2(1,i),g3(1,i),g1(2,i) &
,g2(2,i),g3(2,i) ,g1(3,i),g2(3,i),g3(3,i),g1(4,i),g2(4,i),g3(4,i))
end if
end do
end if ! end ab initio block 


allocate(num_type(itype))

do i=1,itype
k=ntype(i)
icount=0
do j=1,numat
if(zeff(j)==k)icount=icount+1
end do
num_type(i)=icount
end do

if(master)print*,''
if(master)print*,'There are',itype,'unique atomic species in this molecule:'
do i=1,itype
if(master)write(*,*)'Number of ',periodic(ntype(i)),' atoms = ',num_type(i)
end do
if(master)write(*,*)'There are ',nsparkle,' effective atoms in this system.'

if(.not. abinitio)then
do i=1,nsparkle
if(master)print*,zeffsp(i),xsparkle(i),ysparkle(i),zsparkle(i)
end do
end if

call atoms_cpu(numat)


! determine number of basis fns, one center pairs, and size of required matrices

if(.not. abinitio)then
allocate(pairs(itype))
num_pairs=0
num_basis=0
do i=1,itype
pairs(i)=nbas(i)*(nbas(i)+1)/2 
num_basis=num_type(i) * nbas(i)+num_basis
num_pairs=num_type(i) * pairs(i)+num_pairs
end do
tt=float(num_pairs)
tt=tt*(tt+1)/2
end if

if(master)print*,''
if(master)print*,'     SUMMARY OF INPUT KEYWORDS AND SIZE        '
if(master)print*,'     OF ONE AND TWO ELECTRON MATRICES          '
if(master)print*,'*********************************************************************'
if(master)write(*,*)'Keywords set by user: ',keyword
if(master)write(*,*)'Total number of functions in computational basis =',num_basis
if(master)write(*,*)'Total number of one center AO pairs =',num_pairs
if(master)write(*,*)'Required memory for one-electron integral vectors = ' &
,8 * num_basis*(num_basis+1)/2,' bytes'
if(master)write(*,*)'Required memory for two-electron integral vector  =', &
nint(tt)*8,' bytes'
if(master)print*,'End Summary**********************************************'
if(master)print*,''
if(.not. abinitio)then
if(master)then
write(*,*)'*****************************************************************'
write(*,*)'            ',method,' Hamiltonian will be used                  '
write(*,*)'*****************************************************************'
end if
else
if(master)write(*,*)'*****************************************************************'
if(master)write(*,*)'            First Principles Hamiltonian will be used                  '
if(master)write(*,*)'*****************************************************************'
end if

if(.not. abinitio)then ! 1009 block
! get dimension of vectors to store unique integrals
ndim1=num_basis*(num_basis+1)/2
ndim2=num_pairs*(num_pairs+1)/2
!allocate(twoe(ndim2))
!!allocate(S(ndim1))
jmaxints=ERIBUFF*1000/8
if(myrank.eq.0)print*,'Buffer can hold=',1000*jmaxints,'integrals'
allocate(twoe(jmaxints*1000))
allocate(tpair(2,jmaxints*1000))
! get indices of bf's assigned to each atom
allocate(ifirstbf(numat))
allocate(ilastbf(numat))
allocate(ifirst_neighbor(numat))
ifirstbf(1)=1
ilastbf(1)=ifirstbf(1)+nbas(species(1))-1
do i=2,numat
ifirstbf(i)=ilastbf(i-1)+1
ilastbf(i)=ifirstbf(i)+nbas(species(i))-1
end do
!do i=1,numat
!print*,'bf',i,ifirstbf(i),ilastbf(i)
!end do


allocate(num_neighbors(numat))
allocate(ifirstbf2(numat))
allocate(ilastbf2(numat))
ifirstbf2(1)=1
ilastbf2(1)=ifirstbf2(1)+pairs(species(1))-1
do i=2,numat
ifirstbf2(i)=ilastbf2(i-1)+1
ilastbf2(i)=ifirstbf2(i)+pairs(species(i))-1
end do




!do i=1,numat
!print*,'bf2',i,ifirstbf2(i),ilastbf2(i)
!end do







!  compute storage indices for one and two electron integrals
! element (i,j) is stored as i+j(j-1)/2 where i is strictly
! less than or equal to j. in other words, i will always work
! with the upper triangle of the matrix
allocate(offset1(num_basis))
allocate(offset2(num_pairs))
do i=1,num_basis
offset1(i)=i*(i-1)/2
end do
do i=1,num_pairs
offset2(i)=i*(i-1)/2
end do






! read in external parameters if necessary
if(parameters)call loadparam

!  compute the one center integrals for an SPD basis analytically
if(method=='MNDOD')then
 do i=1,itype
 if(nbas(i)>4)then
 call inighd(i,nqs(ntype(i)),nqp(ntype(i)),nqd(ntype(i)),zsone(i),zpone(i),zdone(i))
 end if
 end do
end if
! determine charge separations and additive terms

call ptchg
if(myrank.eq.0)call parout
! thiels d orbital package uses a scaled value of d3 and d5 for computational reasons.
d3=d3*sqrt(two)
d5=d5*sqrt(two)



if(master)write(*,*)''
if(master)write(*,*)'Computing one and two electron integrals...'
call cpusec(time1)
!S=zero
twoe=zero
! put 1's on the diagonal of the overlap matrix
!do i=1,num_basis
!S(i+offset1(i))=one
!end do
! get memory for hcore matrix

!allocate(H(ndim1))
!pick 100MB arbitrarily
hdim_max=12500000
if(master)print*,'Hcore buffer can hold max of',hdim_max,' integrals. Requires',hdim_max*8/1000000,'MB on each core.' 
allocate(H(hdim_max))
allocate(hpair(2,hdim_max))
H=zero



!map is a global matrix of storage indices (see tablesmodule.f90)
! for the integrals - don't screw with it
!map=map-1

! build list that maps 2eri pairs to basis function indices. will need to build fock matrix 
allocate(map_pairs(2,num_pairs))

ii=0
do i=1,numat
   icount=0
   do j=1,nbas(species(i))
   do k=1,j
   icount=icount+1
      do jcol=1,9
      do icol=1,9
      if(map(icol,jcol).eq.icount)goto 40
      end do
      end do
 40   ii=ii+1
      map_pairs(1,ii)=ifirstbf(i)+icol-1
      map_pairs(2,ii)=ifirstbf(i)+jcol-1
   end do
   end do
end do


map=map-1



call two_electron
call hcore
call cpusec(time2)
if(master)write(*,200)time2-time1
if(master)write(*,*)
200 format(' Elapsed time computing integrals (sec) =  ',F10.3)
!
if(master)write(*,*)'Number of electrons for this system = ',nelectrons
if(master)write(*,*)'System charge =',sys_charge
if(master)print*,'Spin Multiplicity = ',int(mult)
iunpair=int(mult)-1
nclose=nelectrons-iunpair
numalpha=nclose/2+iunpair
numbeta=nclose/2
if(master)print*,'# Alpha electrons =',numalpha
if(master)print*,'# Beta  electrons =',numbeta


!call occupied_cpu

end if ! this is end of block labeled 1009 above


!ab initio
if(abinitio)then

call cpusec(time1)
if(master)then
open(unit=45,file='charge')
write(45,*)sys_charge
write(45,*)lindep_tol
close(45)
end if
call get_basis(x,y,z,numat,basis_set,zeff,ERIBUFF,debug)
call pick2
call cpusec(time2)
if(master)print*,'@readinput: Elapsed integral time',time2-time1
!print*,'Entering SCF...'
call scf_abinitio
call diff_pick2
goto 123
end if



if(myrank.eq.0)write(*,*)''
call cpusec(time1)
if(.not.OPTIMIZE)THEN
    TRIAL=.false.
    keep=.true.
    if(restricted)then

       if(CGSCF)then
          call rhf_conjugate
          else
        call rhf
        end if



    else
        call uhf
     end if
    call cpusec(time2)
if(master)    write(*,201)time2-time1
    write(*,*)
   201 format(' Elapsed time for SCF calculation (sec) =  ',F10.3)
!call dipole
   TRIAL=.false.
   keep=.false.


END IF

! the dft integrations have to be before numerical gradient part !!!!!!!!





if(GRADIENT)then
   if(restricted)then
   call find
   else
   call finduhf
   end if


end if

if(OPTIMIZE)then
save_tree=.true.
   if(XYZ)then
! print out initial values of optimized coordinates
print*,'The following atoms have at least one optimizable coordinate:'
do i=1,numat
if(optx(i).eq.1 .or. opty(i).eq.1  .or. optz(i).eq.1)then
WRITE(*,21)i,ZEFF(I),X(I),optx(i),Y(I),opty(i),Z(I),optz(i) 
end if
end do
21 format(i5,i15,f20.5,i2,f20.5,i2,f20.5,i2)

! get total number of optimizable coords.
jopt=0
do i=1,numat
if(optx(i).eq.1)jopt=jopt+1
if(opty(i).eq.1)jopt=jopt+1
if(optz(i).eq.1)jopt=jopt+1
end do
print*,'Total number of optimizable cartesian coordinates',jopt





   allocate(g(3*numat))
   call bfgs(x,y,z,numat,jopt,OPTTOL,optx,opty,optz)
deallocate(ifirst)
deallocate(ifirst2)
deallocate(ilast)
deallocate(ilast2)
call two_electron
call hcore

print*,'Computing energy, dipole moment, and gradient for optimized geometry...'
save_tree=.false.
if(restricted)then
call rhf
call find
else
call uhf
call finduhf
end if

   write(*,*)'OPTIMAL GEOMETRY'
   do i=1,numat
      write(*,300)zeff(i),x(i),y(i),z(i)
   end do
! call dipole
else
     nopt=0
    do j=1,3
    do i=1,ninput
    if(opt(i,j)==1)nopt=nopt+1
    end do
    end do

    allocate(internal(2,nopt))
    icount=0
    icount2=0
    do j=1,3
    do i=1,ninput
    if(opt(i,j)==1)then
    icount=icount+1
    internal(1,icount)=i
    internal(2,icount)=j
    if(j.le.2)icount2=icount2+1
    end if
    end do
    end do
    print*,'there are ',icount2,'optimizable bonds and angles'

if(.not. TS)then
   !call conjugate(ninput,nopt,OPTTOL,q,internal,icount2)
   !else
   call bfgsq(ninput,nopt,OPTTOL,q,internal,icount2,numint,intadd)   
   !end if

else
call make_hess2
end if


call geometry
call nodummy
deallocate(ifirst)
deallocate(ifirst2)
deallocate(ilast)
deallocate(ilast2)
call two_electron
call hcore
save_tree=.false.
print*,'Computing energy, dipole moment and forces for optimized geometry...'
if (restricted)then
call rhf
call find
else
print*,trial,keep
trial=.false.
keep=.true.
call uhf
trial=.false.
keep=.false.
call finduhf
end if

! for srp
open(unit=12,file='opt')
write(12,*)'1'
write(12,*)rsrp
close(12)
open(unit=10,file='geometry')
rewind 10
do i=2,ninput
write(10,*)q(i,1)
end do
do i=3,ninput
write(10,*)q(i,2)*180.0d0/pi
end do
do i=4,ninput
write(10,*)q(i,3)*180.0d0/pi
end do
close(10)
! end for srp


     write(*,*)'OPTIMAL GEOMETRY'
     do i=1,ninput
     write(*,301)zstore(i),q(i,1),opt(i,1),q(i,2)*180.0d0/pi,opt(i,2),q(i,3)*180.0d0/pi,opt(i,3),ref(i,1),ref(i,2),ref(i,3)
     end do


call geometry

write(*,*)'           CARTESIAN ATOMIC COORDINATES ON EXIT: (ANGSTROMS):'
write(*,*)'           ---------------------------------------------------'
write(*,*)'ATOM','       ATOMIC NO.        ','    X','                   Y','                  Z'

                do i=1,numat
                 WRITE(*,20)i,ZEFF(I),X(I),Y(I),Z(I)
                end do


!call dipole
  

end if  


end if

if(DFTHF)then

! redetermine integrals and density if finite difference gradient has been done
if(GRADIENT)then
deallocate(ifirst)
deallocate(ifirst2)
deallocate(ilast)
deallocate(ilast2)
call two_electron
call hcore
call rhf
end if
   if(restricted)then
   call hfdft
   else
!   call hfdftuhf
   end if
end if



if(FREQUENCY)THEN
if(allocated(g))deallocate(g)
allocate(g(3*numat))
call make_hess
END IF


if(MODES)then
call make_hess2
end if


if(IRC)then
if(allocated(g))deallocate(g)
allocate(g(3*numat))
call mep
!call mep2
end if


 call date_and_time(date,time)


123 continue
if(master) write(*,*)'Calculation ended at: ',time(1:2),':',time(3:4),':',time(5:6) &
 ,' on ',date(5:6) ,'/', date(7:8),'/',date(1:4)

return

300 format(i3,f10.5,f10.5,f10.5)
301 format(i3,f10.5,i3,f10.5,i3,f10.5,i3,i3,i3,i3)
end subroutine readin




