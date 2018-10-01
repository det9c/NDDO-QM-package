subroutine mep
use constants
use tables
use scratch_array
use control
use indices
implicit double precision (a-h,o-z)

interface
subroutine scfopt(x0,y0,z0,gout,natoms,ff)
double precision,intent(inout) ::ff
double precision, dimension(natoms),intent(inout) ::x0,y0,z0
double precision, dimension(3*natoms),intent(inout) ::gout
end subroutine scfopt
end interface


double precision,allocatable,dimension(:)::gplus,stepx,stepy,stepz,amass,k1,k2,k3,k4
double precision,dimension(numat)::xmw,ymw,zmw,xold,yold,zold,xinit,yinit,zinit
double precision,allocatable,dimension(:,:)::energies
i3n=3*numat
numint=3*numat-6
allocate(gplus(i3n))
allocate(k1(i3n))
allocate(k2(i3n))
allocate(k3(i3n))
allocate(k4(i3n))
allocate(stepx(numat))
allocate(stepy(numat))
allocate(stepz(numat))
allocate(amass(i3n))

! control variables
deltainit=5.0D-4
deltas=1.0D-2
nsteps=50
open(unit=12,file='mep')
read(12,*)deltainit,deltas,nsteps
close(12)


ilen=2*nsteps+1
allocate(eneRgies(ilen,3))
energies=zero


! if want to make sam file allocate some arrays
if(SAM)then
allocate(rpxyz(i3n,ilen))
allocate(rpdxyz(i3n,ilen))
iupper=i3n*(i3n+1)/2
allocate(rpd2xyz(iupper,ilen))
end if



! evaluate initial hessian and get transition vector
call make_hess


print*,'-----------------------------------------------'
print*,'|          Intrinsic Reaction Coordinate       |'
print*,'|                                              |'
print*,'|                                              |'
print*,'-----------------------------------------------'
print*,'initial variables:'
print*,'deltas=',deltas
print*,'nsteps=',nsteps
print*,'initial delta',deltainit
if(euler)print*,'Euler integration'
if(rk4)print*,'RK4 integration'

! pick up atomic masses
do i=1,numat
istart=3*i-2
amass(istart)=dsqrt( atmass( zeff(i) ) )
amass(istart+1)=dsqrt( atmass( zeff(i) ) )
amass(istart+2)=dsqrt( atmass( zeff(i) ) )
end do
!amass=1.0d0/dsqrt(amass)

amass=1.0d0/amass


!print*,'Transition eigenvector:'
rnorm=zero
   do j=1,i3n
      print*,j,transvec2(j)
    rnorm=rnorm+transvec2(j)*transvec2(j)
   end do
  print*,'transvec norm is',dsqrt(rnorm)
jmax=1
! get initial energy and gradient
!call geometry
!call nodummy
xinit=x
yinit=y
zinit=z
open(unit=1,file='atoms')
read(1,*)iii,jjj
close(1)
rsave=dsqrt((x(iii)-x(jjj))**2+(y(iii)-y(jjj))**2+(z(iii)-z(jjj))**2)
do k=1,numat
xinit(k)=dsqrt( atmass( zeff(k) ) )*xinit(k) !save initial mass weighted coords to compute path length
yinit(k)=dsqrt( atmass( zeff(k) ) )*yinit(k)
zinit(k)=dsqrt( atmass( zeff(k) ) )*zinit(k)
end do
call scfopt(x,y,z,gplus,numat,ff)
rnorm=ddot(i3n,gplus,1,gplus,1)
rnorm=dsqrt(rnorm)
energies(nsteps+1,1)=zero
energies(nsteps+1,2)=ff
energies(nsteps+1,3)=zero

do k=1,numat
print*,x(k),y(k),z(k)
end do

if(SAM)then 
 ircpt=nsteps+1
 call make_hess_sam
 do k=1,numat
  istart=3*k-2
  rpxyz(istart,ircpt)=x(k)
  rpxyz(istart+1,ircpt)=y(k)
  rpxyz(istart+2,ircpt)=z(k)
 end do
 call dcopy(i3n,gplus,1,rpdxyz(1,ircpt),1)
end if


print*,'grad'
do k=1,numat
istart=3*k-2
print*,gplus(istart),gplus(istart+1),gplus(istart+2)
end do



transvec2=transvec2*deltainit ! scale initial transition vector

do iphase=1,2            
if(iphase.eq.1)then
idir=1.0d0
!mass weight initial coordinates
do k=1,numat
istart=3*k-2
xmw(k)=xinit(k)+transvec2(istart)
ymw(k)=yinit(k)+transvec2(istart+1)
zmw(k)=zinit(k)+transvec2(istart+2)
end do

else
idir=-1.0d0
!mass weight initial coordinates
do k=1,numat
istart=3*k-2
xmw(k)=xinit(k)-transvec2(istart)
ymw(k)=yinit(k)-transvec2(istart+1)
zmw(k)=zinit(k)-transvec2(istart+2)
end do

end if









SAVE_TREE=.true.
rlast=rsave
do i=1,nsteps
if(mod(i,10).eq.0)print*,i

  if(i.eq.1)then
 densityout=.true.
 densityin=.false.
 else
 densityout=.true.
 densityin=.true.
 end if

tot1=zero
tot2=zero
tot3=zero
do j=1,numat
tot1=tot1+(xinit(j)-xmw(j))**2+(yinit(j)-ymw(j))**2+(zinit(j)-zmw(j))**2
end do
tot1=dsqrt(tot1)


   do j=1,numat
      x(j)=xmw(j)/dsqrt( atmass( zeff(j) ) )
      y(j)=ymw(j)/dsqrt( atmass( zeff(j) ) )
      z(j)=zmw(j)/dsqrt( atmass( zeff(j) ) )
   end do
!diis=.false.
call scfopt(x,y,z,gplus,numat,ff)
!get gradient norm/
 
rnorma=ddot(i3n,gplus,1,gplus,1)
rnorma=dsqrt(rnorma)

!print*,'rnorm',rnorma
rhere=dsqrt((x(iii)-x(jjj))**2+(y(iii)-y(jjj))**2+(z(iii)-z(jjj))**2)

write(*,49)i,rhere,rhere-rlast,rnorma
rlast=rhere
49 format(i2,f10.5,f10.5,f10.5)




if(SAM)then
ircpt=(nsteps+1)+idir*i
!densityin=.true.
!densityout=.false.
!if(nsteps.gt.5)scftol=1d-3
call make_hess_sam
!densityout=.true.
do k=1,numat
istart=3*k-2
rpxyz(istart,ircpt)=x(k)
rpxyz(istart+1,ircpt)=y(k)
rpxyz(istart+2,ircpt)=z(k)
end do
call dcopy(i3n,gplus,1,rpdxyz(1,ircpt),1)
end if






if(idir.eq.1)then !if forward
energies(nsteps+1+i,1)=deltas*float(i)
energies(nsteps+1+i,2)=ff
energies(nsteps+1+i,3)=tot1
else
energies(nsteps+1-i,1)=-deltas*float(i)
energies(nsteps+1-i,2)=ff
energies(nsteps+1-i,3)=-tot1
end if

!mass weight gragient
gplus=gplus*amass
rnorm=ddot(i3n,gplus,1,gplus,1)
rnorm=dsqrt(rnorm)
rnorm=one/rnorm
gplus=gplus*rnorm
if(EULER)then
gplus=-gplus
gplus=gplus*deltas

elseif(RK4)then
call dcopy(numat,xmw,1,xold,1)
call dcopy(numat,ymw,1,yold,1)
call dcopy(numat,zmw,1,zold,1)
gplus=-gplus
k1=gplus*deltas ! k1

do k=1,numat
istart=3*k-2
xmw(k)=xold(k)+half*k1(istart)
ymw(k)=yold(k)+half*k1(istart+1)
zmw(k)=zold(k)+half*k1(istart+2)
end do
   do k=1,numat
      x(k)=xmw(k)/dsqrt( atmass( zeff(k) ) )
      y(k)=ymw(k)/dsqrt( atmass( zeff(k) ) )
      z(k)=zmw(k)/dsqrt( atmass( zeff(k) ) )
   end do
call scfopt(x,y,z,gplus,numat,ff)
gplus=gplus*amass
rnorm=ddot(i3n,gplus,1,gplus,1)
rnorm=dsqrt(rnorm)
rnorm=one/rnorm
k2=-gplus*rnorm*deltas
goto 123
do k=1,numat
istart=3*k-2
xmw(k)=xold(k)+half*k2(istart)
ymw(k)=yold(k)+half*k2(istart+1)
zmw(k)=zold(k)+half*k2(istart+2)
end do
   do k=1,numat
      x(k)=xmw(k)/dsqrt( atmass( zeff(k) ) )
      y(k)=ymw(k)/dsqrt( atmass( zeff(k) ) )
      z(k)=zmw(k)/dsqrt( atmass( zeff(k) ) )
   end do

call scfopt(x,y,z,gplus,numat,ff)
gplus=gplus*amass
rnorm=ddot(i3n,gplus,1,gplus,1)
rnorm=dsqrt(rnorm)
rnorm=one/rnorm
k3=-gplus*rnorm*deltas

do k=1,numat
istart=3*k-2
xmw(k)=xold(k)+k3(istart)
ymw(k)=yold(k)+k3(istart+1)
zmw(k)=zold(k)+k3(istart+2)
end do
   do k=1,numat
      x(k)=xmw(k)/dsqrt( atmass( zeff(k) ) )
      y(k)=ymw(k)/dsqrt( atmass( zeff(k) ) )
      z(k)=zmw(k)/dsqrt( atmass( zeff(k) ) )
   end do

call scfopt(x,y,z,gplus,numat,ff)
gplus=gplus*amass
rnorm=ddot(i3n,gplus,1,gplus,1)
rnorm=dsqrt(rnorm)
rnorm=one/rnorm
k4=-gplus*rnorm*deltas

gplus=(one/6.0d0)*k1  + (one/three) * k2 + (one/three) * k3 + (one/6.0d0) * k4



123 gplus=k2



call dcopy(numat,xold,1,xmw,1)
call dcopy(numat,yold,1,ymw,1)
call dcopy(numat,zold,1,zmw,1)
end if







!do j=1,numat
!istart=3*j-2
!x(j)=x(j)+gplus(istart)
!y(j)=y(j)+gplus(istart+1)
!z(j)=z(j)+gplus(istart+2)
!end do


do j=1,numat
istart=3*j-2
xmw(j)=xmw(j)+gplus(istart)
ymw(j)=ymw(j)+gplus(istart+1)
zmw(j)=zmw(j)+gplus(istart+2)
end do


end do

write(*,*)'ATOM','       ATOMIC NO.        ','    X','                   Y','                  Z'
                do j=1,numat
!                 WRITE(*,20)j,ZEFF(j),X(j),Y(j),Z(j)
                 WRITE(*,21)periodic(zeff(j)),X(j),Y(j),Z(j)
                end do
                20 format(i5,i15,f20.5,f20.5,f20.5)
                21 format(A,f20.5,f20.5,f20.5)
print*,'norm for this point',rnorma




end do




open(unit=13,file='plot')
do i=1,ilen
write(13,*)energies(i,1),energies(i,3),energies(i,2)*627.51d0/27.21d0
end do
close(13)


if(SAM)then
open(unit=13,file='sam.irc')
iloop=2*nsteps+1
do i=1,iloop
write(13,*)i
do j=1,numat
istart=3*j-2
write(13,30)rpxyz(istart,i),rpxyz(istart+1,i),rpxyz(istart+2,i)
end do
write(13,*)i
do j=1,numat
istart=3*j-2
write(13,30)rpdxyz(istart,i),rpdxyz(istart+1,i),rpdxyz(istart+2,i)
end do
write(13,*)i
do j=1,numat*(i3n+1)/2
istart=3*j-2
write(13,30)rpd2xyz(istart,i),rpd2xyz(istart+1,i),rpd2xyz(istart+2,i)
end do

end do

30 format(f20.15,f20.15,f20.15)
end if






return
end subroutine mep






