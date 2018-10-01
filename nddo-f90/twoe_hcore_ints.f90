 subroutine twoe_hcore_ints(i,j,wout)
  use constants
  use indices
  use control
  use scratch_array
  use tables
 implicit double precision (a-h,o-z)
 include 'mpif.h'
 integer,intent(in)::i,j
 double precision,dimension(:,:),intent(inout)::wout

 double precision,dimension(22) :: twoelcl(22)
 double precision,dimension(10,10)::twoemol
 double precision,dimension(5)::overlcl
 double precision,dimension(4,4)::overmol
 double precision,dimension(3,2)::coor
 double precision,dimension(15,45)::YY
 double precision,allocatable,dimension(:,:)::WW,W2,dlocal,smatrix
 logical,allocatable,dimension(:,:)::r2cent
 double precision,dimension(8,2)::po
 double precision,dimension(6,2)::dd
  double precision,dimension(10)::betai,betaj 
 
!********************************************
!         this section lists the interfaces for all the subroutines

interface 
   subroutine three_core(angular)
double precision,intent(inout)::angular
end subroutine three_core

  subroutine local(twoelcl,ni,nj,d1a,d2a,d1b,d2b,p0a,p0b,p1a,p1b,p2a,p2b,rsq,r) 
  double precision,intent(in) ::d1a,d2a,d1b,d2b,p0a,p0b,p1a,p1b,p2a,p2b,rsq,r
  integer,intent(in) ::ni,nj
  double precision , dimension(22),intent(inout) ::twoelcl
  end subroutine local

  subroutine quadrant(x1,y1,z1,phi,sinth)
  double precision, intent(inout)::phi,x1,y1,z1,sinth
  end subroutine quadrant

  subroutine rotate(twoemol,theta,phi,twoelcl,overlcl,overmol)
  double precision,dimension(:,:),intent(inout) ::twoemol,overmol
  double precision,intent(in)::theta,phi
  double precision,dimension(:),intent(in)::twoelcl,overlcl
  end subroutine rotate

  subroutine overlap(iatom,jatom,iorbs,jorbs,rau,overlcl)
  integer,intent(in)::iatom,jatom,iorbs,jorbs
  double precision,dimension(:),intent(inout)::overlcl
  double precision,intent(in)::rau
  end subroutine overlap

 subroutine repel(nuclear,r,iatom,jatom,ssss)
double precision,intent(inout)::nuclear,r
double precision,intent(in)::ssss
integer,intent(in)::iatom,jatom
 end subroutine repel

subroutine rearrange(WW,ipair,jpair)
double precision,intent(inout),dimension(:,:)::WW
integer,intent(in)::ipair,jpair
end subroutine rearrange

      SUBROUTINE WSTORE (ww,NI,MODE,ip)
            integer,intent(in)::mode,ip,NI
          double precision,dimension(:,:),intent(inout)::WW
      end subroutine wstore

subroutine doverlap(iatom,jatom,rau,dlocal)
double precision,intent(in)::rau
integer,intent(in)::iatom,jatom
double precision,intent(inout),dimension(:,:)::dlocal
end subroutine doverlap

end interface
!  end  section 
!********************************************8   




!***********************************************



wout=zero

 ni=nbas(species(i))

        if(pbc )then
             sxij=fracs(1,i)-fracs(1,j)
             syij=fracs(2,i)-fracs(2,j)
             szij=fracs(3,i)-fracs(3,j)
             sxij=sxij-anint(sxij)
             syij=syij-anint(syij)
             szij=szij-anint(szij)
             rxij=cell(1,1)*sxij+cell(1,2)*syij+cell(1,3)*szij
             ryij=cell(2,1)*sxij+cell(2,2)*syij+cell(2,3)*szij
             rzij=cell(3,1)*sxij+cell(3,2)*syij+cell(3,3)*szij
       end if
 if(.not. pbc)then
   rsq=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
 else
 rsq=rxij*rxij + ryij*ryij + rzij*rzij
 end if
   r=dsqrt(rsq)
   if(r>cutoff)then
   go to 90
   end if
   r=r/autoang
   rau=r
   rsq=r*r
   nj=nbas(species(j))

if(METHOD=='MNDOD')then
allocate(r2cent(45,45))
call bdata4(r2cent)
end if



! compute the 22 nonzero integrals for an sp basis set in local coordinates where 
! atom i is always at the origin and atom j is up the positive z axiz which fixes
! the phase in the local coordinates
twoelcl=zero
! do the sp portion
if(ni>4)then
isp=4
else
isp=ni
end if
if(nj>4)then
jsp=4
else
jsp=nj
end if
 call local(twoelcl,isp,jsp,d1(species(i)),d2(species(i)),d1(species(j)),d2(species(j)) &
,p0(species(i)),p0(species(j)),p1(species(i)),p1(species(j)),p2(species(i)),p2(species(j)),rsq,r)
!  the 22 local integrals are stored in twoelcl in the order given in local.f90
!  now we need to compute the rotation matrix from local to molecular coordinates
! compute   r,theta,phi, and orientation (quadrant) of two atoms
 r=r*autoang

if(.not. pbc)then
 x1=x(j)-x(i)
 y1=y(j)-y(i)
 z1=z(j)-z(i)
else
x1=-rxij
y1=-ryij
z1=-rzij
end if
 costh=z1/r
if(abs(costh).gt.one)costh=dsign(one,costh)
 theta=dacos(costh)
 sinth=dsin(theta)

! special case if two atoms both lie on the x,y,or z axes. Because of my not so cleverly
! chosen convention for rotation, i have to handle these cases separately.
!print*,'atom ji',j,i
!if on + or - z axis,then
 if(abs(sinth)<1D-10)sinth=zero
 if(abs(sinth)==zero)then  ! there is a problem here with machine precision. 1d-20 effectively turns
! this check off
     phi=zero
     if(z1.gt.zero)then
     theta=zero
     else
     theta=pi
     end if
! if not on z axis then...
else
         cosphi=x1/(r*sinth)
         if(abs(cosphi).gt.one)cosphi=dsign(one,cosphi)
         phi=dacos(cosphi)
         cosphi=dcos(phi)
!         print*,'so phi is',phi*180.0d0/pi

! if on y axis or in plane of y axis then
         if(abs(cosphi)<1D-10)cosphi=zero
         if(cosphi==zero)then
              if(y1.gt.zero)then
              phi=pi/two
              else
              phi=three*pi/two
              end if
         elseif(y1==zero.and.z1==zero)then
               if(x1<zero)then
                  phi=pi
                  else
                  phi=zero
                  end if
         else
         call quadrant(x1,y1,z1,phi,sinth)
          end if


end if


! now we have r,theta, and phi in local coordinates.  use
! these angles to rotate the integrals counter clockwise about z axis thru phi
! followed by a counter clockwise rotation thru theta about y axis
 
! twoemol will return a 10 x 10 matrix of two electron repulsion integrals between
! atoms i and j in molecular coordinates.  these will then be placed into the two electron
! integral vector
! do the overlaps for this pair as well
!parallel call overlap(i,j,isp,jsp,rau,overlcl)

call rotate(twoemol,theta,phi,twoelcl,overlcl,overmol)

! do the d orbital integrals if necessary)
if( ni>4  .or. nj>4)then
! pass the coordinates in a 2 column array
coor(1,1)=x(i)
coor(2,1)=y(i)
coor(3,1)=z(i)
coor(1,2)=x(j)
coor(2,2)=y(j)
coor(3,2)=z(j)
dd=zero
po=zero

po(1,1)=p0(species(i))
po(2,1)=p1(species(i))
po(3,1)=p2(species(i))
po(4,1)=p3(species(i))
po(5,1)=p4(species(i))
po(6,1)=p5(species(i))
po(7,1)=p6(species(i))
po(8,1)=p7(species(i))
dd(1,1)=zero
dd(2,1)=d1(species(i))
dd(3,1)=d2(species(i))
dd(4,1)=d3(species(i))
dd(5,1)=d4(species(i))
dd(6,1)=d5(species(i))

po(1,2)=p0(species(j))
po(2,2)=p1(species(j))
po(3,2)=p2(species(j))
po(4,2)=p3(species(j))
po(5,2)=p4(species(j))
po(6,2)=p5(species(j))
po(7,2)=p6(species(j))
po(8,2)=p7(species(j))
dd(1,2)=zero
dd(2,2)=d1(species(j))
dd(3,2)=d2(species(j))
dd(4,2)=d3(species(j))
dd(5,2)=d4(species(j))
dd(6,2)=d5(species(j))
 call rotmat(1,2,ni,nj,2,coor,rau,yy)
 allocate(ww(pairs(species(i)),pairs(species(j))))
ww=zero
 CALL REPPD (2,1,rau,WW,pairs(species(j)),pairs(species(i)),1,dd,po)
 CALL ROTD  (WW,YY,pairs(species(j)),pairs(species(i)),r2cent)
call rearrange(WW,pairs(species(i)),pairs(species(j)))


!parallel
!if(allocated(dlocal))deallocate(dlocal)
!allocate(dlocal(ni,nj))
!call doverlap(i,j,rau,dlocal)
end if









! compute the core-core repulsion for this pair which requires two electron integral (ss|ss)
!parallel call repel(t1,r,i,j,twoemol(1,1))
!enuc=enuc+t1
! put the overlap chunk into the overlap vector
!allocate(smatrix(ni,nj))
!betai(1)=betas(species(i))
!if(ni>1)then
!betai(2)=betap(species(i))
!betai(3)=betap(species(i))
!betai(4)=betap(species(i))
!end if
!if(ni>4)then
!betai(5)=betad(species(i))
!betai(6)=betad(species(i))
!betai(7)=betad(species(i))
!betai(8)=betad(species(i))
!betai(9)=betad(species(i))
!end if

!betaj(1)=betas(species(j))
!if(nj>1)then
!betaj(2)=betap(species(j))
!betaj(3)=betap(species(j))
!betaj(4)=betap(species(j))
!end if
!if(nj>4)then
!betaj(5)=betad(species(j))
!betaj(6)=betad(species(j))
!betaj(7)=betad(species(j))
!betaj(8)=betad(species(j))
!betaj(9)=betad(species(j))
!end if



!smatrix(1:isp,1:jsp)=overmol(1:isp,1:jsp)




itop=isp*(isp+1)/2
jtop=jsp*(jsp+1)/2
allocate(w2(pairs(species(i)),pairs(species(j))))
W2(1:itop,1:jtop)=twoemol(1:itop,1:jtop)
if(method=='MNDOD')then
  if(nj>4)then
   w2(1:itop,jtop+1:pairs(species(j)))=ww(1:itop,jtop+1:pairs(species(j)))
!   smatrix(1:isp,jsp+1:nj)=dlocal(1:isp,jsp+1:nj)
   end if
  if(ni>4)then
  w2(itop+1:pairs(species(i)),1:jtop)=ww(itop+1:pairs(species(i)),1:jtop)
!  smatrix(isp+1:ni,1:jsp)=dlocal(isp+1:ni,1:jsp)
  end if
  if(ni>4 .and. nj >4)then
  w2(itop+1:pairs(species(i)),jtop+1:pairs(species(j)))=ww(itop+1:pairs(species(i)),jtop+1:pairs(species(j)))
!smatrix(isp+1:ni,jsp+1:nj)=dlocal(isp+1:ni,jsp+1:nj)
  end if
 if(allocated(ww))deallocate(ww)
end if


wout=w2
deallocate(w2)



if(method=='MNDOD')deallocate(r2cent)


90 continue


return




  end subroutine twoe_hcore_ints

