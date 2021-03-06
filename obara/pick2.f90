subroutine pick2
use gaussian_basis
implicit double precision(a-h,o-z)
include 'mpif.h'
double precision,dimension(3)::psints,pvector,psbuffer,dkpsints,dkpsbuffer,dnpsints,scr3,dnpsbuffer,qvector,psssints
double precision,dimension(3)::wvector,psss_buffer
double precision,dimension(3,3)::ppints,ppbuffer,dkppints,dkppbuffer,scr4,dnppints,dnppbuffer,psps_buffer,pspsints,&
ppss_buffer,ppssints,dsss_buffer,dsssints
double precision,dimension(3,3,3,3)::pppp_buffer,ppppints,ddbuffer,ddints,dppsints,dpps_buffer,dsppints,dspp_buffer
double precision,dimension(3,3,3)::ppps_buffer,pppsints,dpssints,dpss_buffer,dspsints,dsps_buffer
integer,dimension(1296)::ipppp_write
integer,dimension(4,600)::inds
double precision,dimension(600)::rbuffer
double precision,dimension(3,3)::dsints,dsbuffer,dkdsbuffer,dkdsints,dndsbuffer,dndsints,scrdsn
double precision,dimension(3,3,3)::dpints,dpbuffer,dkdpbuffer,dkdpints,dndpbuffer,dndpints,scrdpn
double precision,dimension(3,3,3,3)::dkddbuffer,dkddints,dnddbuffer,dnddints,scrddn,ddss_buffer,ddssints, &
dsds_buffer,dsdsints 
double precision,dimension(3,3,3,3,3)::dpppints,dppp_buffer,ddpsints,ddps_buffer,dsdpints,dsdp_buffer
double precision,dimension(3,3,3,3,3,3)::ddppints,ddpp_buffer,dpdpints,dpdp_buffer,dddsints,ddds_buffer
double precision,dimension(3,3,3,3,3,3,3)::dddpints,dddp_buffer
double precision,dimension(3,3,3,3,3,3,9)::ddddints,dddd_buffer

double precision,dimension(6)::imapdl,imapdr
logical::iwrite

! routine will compute all matrices (overlap, hcore integrals) depending on 2 basis functions
! name is for historical reasons




interface
subroutine overlap_ss(dist,r,t,ssout)
double precision,intent(in)::dist,r,t
double precision,intent(out)::ssout
end subroutine overlap_ss

subroutine buildp(expa,expb,icola,icolb,pvector)
double precision,intent(in)::expa,expb
integer,intent(in)::icola,icolb
double precision,dimension(3),intent(inout)::pvector
end subroutine buildp

subroutine overlap_ps(distsq,sexp,pexp,ishellp,ishells,psints,pvector)
double precision,dimension(:),intent(inout)::psints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,pexp
integer,intent(in)::ishellp,ishells
end subroutine overlap_ps

subroutine overlap_pp(distsq,pexpa,pexpb,ishellpa,ishellpb,ppints,pvector)
double precision,dimension(:,:),intent(inout)::ppints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexpa,pexpb
integer,intent(in)::ishellpa,ishellpb
end subroutine overlap_pp

subroutine overlap_ds(distsq,sexp,dexp,idshell,isshell,dsints,pvector)
double precision,dimension(:,:),intent(inout)::dsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexp
integer,intent(in)::isshell,idshell
end subroutine overlap_ds


subroutine overlap_dp(distsq,pexp,dexpb,ishellp,ishelld,dpints,pvector)
double precision,dimension(:,:,:),intent(inout)::dpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpb
integer,intent(in)::ishellp,ishelld
end subroutine overlap_dp

subroutine overlap_dd(distsq,dexpa,dexpb,ishella,ishellb,ddints,pvector)
double precision,dimension(:,:,:,:),intent(inout)::ddints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,dexpa,dexpb
integer,intent(in)::ishella,ishellb
end subroutine overlap_dd


subroutine kinetic_ss(dist,r,t,ss,dkssout)
double precision,intent(in)::dist,r,t,ss
double precision,intent(out)::dkssout
end subroutine kinetic_ss

subroutine kinetic_ps(distsq,sexp,pexp,ishellp,ishells,psints,pvector,dkspints)
double precision,dimension(:),intent(in)::psints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,pexp
integer,intent(in)::ishellp,ishells
double precision,intent(inout),dimension(:)::dkspints
end subroutine kinetic_ps

subroutine kinetic_pp(distsq,pexpa,pexpb,ishellpa,ishellpb,ppints,pvector,dkppints)
double precision,dimension(:,:),intent(inout)::dkppints
double precision,dimension(:,:),intent(in)::ppints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexpa,pexpb
integer,intent(in)::ishellpa,ishellpb
end subroutine kinetic_pp

subroutine kinetic_ds(distsq,sexp,dexp,ishelld,ishells,dsints,pvector,dkdsints)
double precision,dimension(:,:),intent(inout)::dkdsints
double precision,dimension(:,:),intent(in)::dsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexp
integer,intent(in)::ishelld,ishells
end subroutine kinetic_ds

subroutine kinetic_dp(distsq,pexp,dexpa,ishellp,ishelld,dpints,pvector,dkdpints)
double precision,dimension(:,:,:),intent(inout)::dkdpints
double precision,dimension(:,:,:),intent(in)::dpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpa
integer,intent(in)::ishelld,ishellp
end subroutine kinetic_dp


subroutine kinetic_dd(distsq,dexpa,dexpb,ishella,ishellb,ddints,pvector,dkddints)
double precision,dimension(:,:,:,:),intent(inout)::dkddints
double precision,dimension(:,:,:,:),intent(in)::ddints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,dexpb,dexpa
integer,intent(in)::ishella,ishellb
end subroutine kinetic_dd






subroutine nuclear_ss(distsq,r,t,ss,pvector,mindex,snsout)
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,r,t,ss
double precision,intent(inout)::snsout
integer,intent(in)::mindex
end subroutine nuclear_ss 

subroutine nuclear_ps(distsq,sexp,pexp,ishellp,ishells,pvector,mindex,dnpsints)
double precision,dimension(:),intent(inout)::dnpsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,pexp
integer,intent(in)::ishellp,ishells,mindex
end subroutine nuclear_ps


subroutine nuclear_pp(distsq,pexpa,pexpb,ishellpa,ishellpb,pvector,mindex,dnppints)
double precision,dimension(:,:),intent(inout)::dnppints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexpa,pexpb
integer,intent(in)::ishellpa,ishellpb,mindex
end subroutine nuclear_pp

subroutine nuclear_ds(distsq,sexp,dexpb,ishelld,ishells,pvector,mindex,dndsints)
double precision,dimension(:,:),intent(inout)::dndsints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,sexp,dexpb
integer,intent(in)::ishells,ishelld,mindex
end subroutine nuclear_ds

subroutine nuclear_dp(distsq,pexp,dexpb,ishellp,ishelld,pvector,mindex,dndpints)
double precision,dimension(:,:,:),intent(inout)::dndpints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,pexp,dexpb
integer,intent(in)::ishellp,ishelld,mindex
end subroutine nuclear_dp 


subroutine nuclear_dd(distsq,dexpa,dexpb,ishella,ishellb,pvector,mindex,dnddints)
double precision,dimension(:,:,:,:),intent(inout)::dnddints
double precision,dimension(:),intent(in)::pvector
double precision,intent(in)::distsq,dexpa,dexpb
integer,intent(in)::ishella,ishellb,mindex
end subroutine nuclear_dd


subroutine buildk(distsq,expa,expb,icola,icolb,dkab,za_plus_zb)
double precision,intent(in)::distsq,expa,expb
integer,intent(in)::icola,icolb
double precision,intent(inout)::dkab,za_plus_zb
end subroutine buildk

subroutine twoe_ssss(mindex,zabcd,za_plus_zb,zc_plus_zd,pvector,qvector,ssss)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
double precision,dimension(:),intent(in)::pvector,qvector
double precision,intent(inout)::ssss
end subroutine twoe_ssss

subroutine twoe_psss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellp,pvector,qvector,wvector,psssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:),intent(inout)::psssints
end subroutine twoe_psss

subroutine twoe_psps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpc,pvector,qvector,wvector,pspsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::pspsints
end subroutine twoe_psps

subroutine twoe_ppss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,pvector,qvector,wvector,ppssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::ppssints
end subroutine twoe_ppss


subroutine twoe_ppps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,ishellpc,pvector,qvector,wvector,pppsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,intent(inout),dimension(3,3,3)::pppsints
end subroutine twoe_ppps


subroutine twoe_pppp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellpa,ishellpb,ishellpc,ishellpd,pvector,qvector,wvector,ppppints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellpa,ishellpb,ishellpc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::ppppints
end subroutine twoe_pppp

subroutine twoe_dsss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,pvector,qvector,wvector,dsssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:),intent(inout)::dsssints
end subroutine twoe_dsss


subroutine twoe_dpss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellp,pvector,qvector,wvector,dpssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:),intent(inout)::dpssints
end subroutine twoe_dpss



subroutine twoe_dsps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellp,pvector,qvector,wvector,dspsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellp
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:),intent(inout)::dspsints
end subroutine twoe_dsps

subroutine twoe_dpps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpc,pvector,qvector,wvector,dppsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellpb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dppsints
end subroutine twoe_dpps

subroutine twoe_dspp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpc,ishellpd,pvector,qvector,wvector,dsppints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellpd,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dsppints
end subroutine twoe_dspp

subroutine twoe_dppp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishelld,ishellpb,ishellpc,ishellpd,pvector,qvector,wvector,dpppints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishelld,ishellpd,ishellpc,ishellpb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:),intent(inout)::dpppints
double precision,dimension(3)::qkkck,wiipi,psssints
end subroutine twoe_dppp


subroutine twoe_ddss(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,pvector,qvector,wvector,ddssints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::ddssints
end subroutine twoe_ddss

subroutine twoe_dsds(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,pvector,qvector,wvector,dsdsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:),intent(inout)::dsdsints
end subroutine twoe_dsds

subroutine twoe_ddps(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishellpc,pvector,qvector,wvector,ddpsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishellpc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:),intent(inout)::ddpsints
end subroutine twoe_ddps


subroutine twoe_dsdp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldc,ishellpd,pvector,qvector,wvector,dsdpints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:),intent(inout)::dsdpints
end subroutine twoe_dsdp

subroutine twoe_ddpp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishellpc,ishellpd,pvector,qvector,wvector,ddppints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishellpc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:),intent(inout)::ddppints
end subroutine twoe_ddpp

subroutine twoe_dpdp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishellpb,ishelldc,ishellpd,pvector,qvector,wvector,dpdpints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishellpb,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:),intent(inout)::dpdpints
end subroutine twoe_dpdp

subroutine twoe_ddds(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,pvector,qvector,wvector,dddsints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishelldc
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:),intent(inout)::dddsints
end subroutine twoe_ddds

subroutine twoe_dddp(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishellpd,pvector,qvector,wvector,dddpints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishelldc,ishellpd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(:,:,:,:,:,:,:),intent(inout)::dddpints
end subroutine twoe_dddp

subroutine twoe_dddd(mindex,zabcd,za_plus_zb,zc_plus_zd,ishellda,ishelldb,ishelldc,ishelldd,pvector,qvector,wvector,ddddints)
integer,intent(in)::mindex
double precision,intent(in)::zabcd,za_plus_zb,zc_plus_zd
integer,intent(in)::ishellda,ishelldb,ishelldc,ishelldd
double precision,dimension(:),intent(in)::pvector,qvector,wvector
double precision,dimension(3,3,3,3,3,3,9),intent(inout)::ddddints
end subroutine twoe_dddd



subroutine pack(indi,indj,ij)
integer,intent(in)::indi,indj
integer,intent(inout)::ij
end subroutine pack




end interface


double precision,dimension(3)::ps_batch
integer::mu,nu


34 format(i6,i6,i6,i6,f15.10)

imapdl(1)=1
imapdl(2)=1
imapdl(3)=2
imapdl(4)=1
imapdl(5)=2
imapdl(6)=3

imapdr(1)=1
imapdr(2)=2
imapdr(3)=2
imapdr(4)=3
imapdr(5)=3
imapdr(6)=3



!do i=1,ishells_total
!do j=1,20
!print*,i,j,shell_coefs(j,i),shell_exp(j,i)
!end do
!end do




call cpusec(timein)
if(myrank.eq.0)print*,''
if(myrank.eq.0)print*,'Computing overlap,kinetic, nuclear attraction and 2-electron integrals....' 

pi=dacos(-1.0d0)
cutoff_ov=1.0d-15
cutoff_hc=1.0d-15
cutoff_2e=1.0d-15
dorbscale=1.0d0/dsqrt(3.0d0)  !!!9.0d0**(-0.25d0)
thresh=1.0d-15
itotal=0
rdist=189.0d0
! set up interplation grid for auxiliary integrals
call gen_grid







if(myrank.eq.0)open(unit=1,file='overlap')
if(myrank.eq.0)open(unit=2,file='hcore')
!!open(unit=3,file='two_eri')


if(myrank.eq.0)then
! loop over (s|s) batches
do ishell=1,num_s_shell
do jshell=1,ishell
distsq=(shell_coor(1,isorb(jshell))-shell_coor(1,isorb(ishell)))**2+(shell_coor(2,isorb(jshell))-shell_coor(2,isorb(ishell)))**2+(shell_coor(3,isorb(jshell))-shell_coor(3,isorb(ishell)))**2
sum=0.
sumkin=0.
sumnuc=0.
do mu=1,num_contr(isorb(ishell))
do nu=1,num_contr(isorb(jshell))
call buildp(shell_exp(nu,isorb(jshell)),shell_exp(mu,isorb(ishell)),isorb(jshell),isorb(ishell),pvector)
call overlap_ss(distsq,shell_exp(nu,isorb(jshell)),shell_exp(mu,isorb(ishell)),ssintegral)
call kinetic_ss(distsq,shell_exp(nu,isorb(jshell)),shell_exp(mu,isorb(ishell)),ssintegral,dkssintegral)

total=0.
do inuc=1,natoms
call nuclear_ss(distsq,shell_exp(nu,isorb(jshell)),shell_exp(mu,isorb(ishell)),ssintegral,pvector,0,dnssintegral)
total=total+dnssintegral
end do

sum=sum+ssintegral*shell_coefs(nu,isorb(jshell))*shell_coefs(mu,isorb(ishell))
sumkin=sumkin+dkssintegral*shell_coefs(nu,isorb(jshell))*shell_coefs(mu,isorb(ishell))
sumnuc=sumnuc+total*shell_coefs(nu,isorb(jshell))*shell_coefs(mu,isorb(ishell))

end do
end do
if(dabs(sum).gt.cutoff_ov)write(1,237)index_basis(1,isorb(ishell)),index_basis(1,isorb(jshell)),sum
if(dabs(sumkin+sumnuc).gt.cutoff_hc)write(2,237)index_basis(1,isorb(ishell)),index_basis(1,isorb(jshell)),sumkin,sumnuc
!print*,'ss',index_basis(1,isorb(ishell)),index_basis(1,isorb(jshell)),sum,sumkin,sumnuc
end do 
end do
!!!!!!!!!!!
!
!speed it up box o stuff
!
! piiai in these integrals may be useful in more than one place.  global or something.  for example
! in the loop over the nuclei, piiai never changes and only need piici 
!
! precompute distance over shells so won't keep computing distsq
!
!!!!!!!!!!!!!!


! loop over (p|s) batches
 do ishell=1,num_p_shell
 do jshell=1,num_s_shell
  distsq=(shell_coor(1,isorb(jshell))-shell_coor(1,iporb(ishell)))**2+&
  (shell_coor(2,isorb(jshell))-shell_coor(2,iporb(ishell)))**2&
  +(shell_coor(3,isorb(jshell))-shell_coor(3,iporb(ishell)))**2
 sum=0.
 psbuffer=0.
 dkpsbuffer=0.
 dnpsbuffer=0.
 do mu=1,num_contr(iporb(ishell))
 do nu=1,num_contr(isorb(jshell))
 call buildp(shell_exp(nu,isorb(jshell)),shell_exp(mu,iporb(ishell)),isorb(jshell),iporb(ishell),pvector)
 call overlap_ps(distsq,shell_exp(nu,isorb(jshell)),shell_exp(mu,iporb(ishell)),iporb(ishell),isorb(jshell),psints,pvector)
 psbuffer=psbuffer+psints*shell_coefs(nu,isorb(jshell))*shell_coefs(mu,iporb(ishell))
 call kinetic_ps(distsq,shell_exp(nu,isorb(jshell)),shell_exp(mu,iporb(ishell)),iporb(ishell),isorb(jshell),psints,pvector,dkpsints)
 dkpsbuffer=dkpsbuffer+dkpsints*shell_coefs(nu,isorb(jshell))*shell_coefs(mu,iporb(ishell))

scr3=0.
do inuc=1,natoms
 call nuclear_ps(distsq,shell_exp(nu,isorb(jshell)),shell_exp(mu,iporb(ishell)),iporb(ishell),isorb(jshell),pvector,0,dnpsints)
scr3=scr3+dnpsints
end do 
dnpsbuffer=dnpsbuffer+scr3*shell_coefs(nu,isorb(jshell))*shell_coefs(mu,iporb(ishell))

 end do
 end do

 do nn=1,3
! print*,'ps',index_basis(nn,iporb(ishell)),index_basis(1,isorb(jshell)),psbuffer(nn),dkpsbuffer(nn),dnpsbuffer(nn)
if(dabs(psbuffer(nn)).gt.cutoff_ov)write(1,237)index_basis(nn,iporb(ishell)),index_basis(1,isorb(jshell)),psbuffer(nn)
if(dabs(dkpsbuffer(nn)+dnpsbuffer(nn)).gt.cutoff_hc)write(2,237)index_basis(nn,iporb(ishell)),index_basis(1,isorb(jshell)),dkpsbuffer(nn),dnpsbuffer(nn)
 end do

 end do
 end do


! loop over (p|p) batches
do ishell=1,num_p_shell
do jshell=1,ishell
 distsq=(shell_coor(1,iporb(ishell))-shell_coor(1,iporb(jshell)))**2+(shell_coor(2,iporb(ishell))-shell_coor(2,iporb(jshell)  ))**2+(shell_coor(3,iporb(ishell))-shell_coor(3,iporb(jshell)))**2
  ppbuffer=0.
  dkppbuffer=0.
  dnppbuffer=0.
 do mu=1,num_contr(iporb(ishell))
 do nu=1,num_contr(iporb(jshell))
 call buildp(shell_exp(mu,iporb(ishell)),shell_exp(nu,iporb(jshell)),iporb(ishell),iporb(jshell),pvector)

 call overlap_pp(distsq,shell_exp(mu,iporb(ishell)),shell_exp(nu,iporb(jshell)),iporb(ishell),iporb(jshell),ppints,pvector)

 call kinetic_pp(distsq,shell_exp(mu,iporb(ishell)),shell_exp(nu,iporb(jshell)),iporb(ishell),iporb(jshell),ppints,pvector,dkppints)

 ppbuffer=ppbuffer+ppints*shell_coefs(mu,iporb(ishell))*shell_coefs(nu,iporb(jshell))
 dkppbuffer=dkppbuffer+dkppints*shell_coefs(mu,iporb(ishell))*shell_coefs(nu,iporb(jshell))


scr4=0.
do inuc=1,natoms
 call nuclear_pp(distsq,shell_exp(mu,iporb(ishell)),shell_exp(nu,iporb(jshell)),iporb(ishell),iporb(jshell),pvector,0,dnppints)
scr4=scr4+dnppints
end do
dnppbuffer=dnppbuffer+scr4*shell_coefs(mu,iporb(ishell))*shell_coefs(nu,iporb(jshell))




 end do
 end do




if(ishell.ne.jshell)then
do nn=1,3
do mm=1,3
!print*,'pp',index_basis(nn,iporb(ishell)),index_basis(mm,iporb(jshell)),dkppbuffer(mm,nn),dnppbuffer(mm,nn)
write(1,237)index_basis(mm,iporb(ishell)),index_basis(nn,iporb(jshell)),ppbuffer(mm,nn)
write(2,237)index_basis(mm,iporb(ishell)),index_basis(nn,iporb(jshell)),dkppbuffer(mm,nn),dnppbuffer(mm,nn)
end do
end do

else
do nn=1,3
do mm=1,nn
!print*,'pp',index_basis(nn,iporb(ishell)),index_basis(mm,iporb(jshell)),dkppbuffer(mm,nn),dnppbuffer(mm,nn)
write(1,237)index_basis(mm,iporb(ishell)),index_basis(nn,iporb(jshell)),ppbuffer(mm,nn)
write(2,237)index_basis(mm,iporb(ishell)),index_basis(nn,iporb(jshell)),dkppbuffer(mm,nn),dnppbuffer(mm,nn)
end do
end do
end if


end do
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute (d|s) one electron integrals
!
!

! loop over (d|s) batches
 do ishell=1,num_d_shell
 do jshell=1,num_s_shell
  distsq=(shell_coor(1,isorb(jshell))-shell_coor(1,idorb(ishell)))**2+&
  (shell_coor(2,isorb(jshell))-shell_coor(2,idorb(ishell)))**2&
  +(shell_coor(3,isorb(jshell))-shell_coor(3,idorb(ishell)))**2
 sum=0.
 dsbuffer=0.
 dkdsbuffer=0.
 dndsbuffer=0.
 do mu=1,num_contr(idorb(ishell))
 do nu=1,num_contr(isorb(jshell))
 call buildp(shell_exp(nu,isorb(jshell)),shell_exp(mu,idorb(ishell)),isorb(jshell),idorb(ishell),pvector)
 call overlap_ds(distsq,shell_exp(nu,isorb(jshell)),shell_exp(mu,idorb(ishell)),idorb(ishell),isorb(jshell),dsints,pvector)
 dsbuffer=dsbuffer+dsints*shell_coefs(nu,isorb(jshell))*shell_coefs(mu,idorb(ishell))
 call kinetic_ds(distsq,shell_exp(nu,isorb(jshell)),shell_exp(mu,idorb(ishell)),idorb(ishell),isorb(jshell),dsints,pvector,dkdsints)
 dkdsbuffer=dkdsbuffer+dkdsints*shell_coefs(nu,isorb(jshell))*shell_coefs(mu,idorb(ishell))

scrdsn=0.0
do inuc=1,natoms
 call nuclear_ds(distsq,shell_exp(nu,isorb(jshell)),shell_exp(mu,idorb(ishell)),idorb(ishell),isorb(jshell),pvector,0,&
dndsints)
scrdsn=scrdsn+dndsints
end do
dndsbuffer=dndsbuffer+scrdsn*shell_coefs(nu,isorb(jshell))*shell_coefs(mu,idorb(ishell))

 end do
 end do

! put in rest of normalization for xx,yy,zz 
do i=1,3
dsbuffer(i,i)=dsbuffer(i,i)*dorbscale
dkdsbuffer(i,i)=dkdsbuffer(i,i)*dorbscale
dndsbuffer(i,i)=dndsbuffer(i,i)*dorbscale
end do

 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
! print*,'ds',index_basis(icount,idorb(ishell)),index_basis(1,isorb(jshell)),dsbuffer(i,j)!,dkpsbuffer(nn),dnpsbuffer(nn)
if(dabs(dsbuffer(i,j)).gt.cutoff_ov)write(1,237)index_basis(icount,idorb(ishell)),index_basis(1,isorb(jshell)),dsbuffer(i,j)
!print*,'kinetic ds',index_basis(icount,idorb(ishell)),index_basis(1,isorb(jshell)),dkdsbuffer(i,j),dndsbuffer(i,j)
if(dabs(dkdsbuffer(i,j)+dndsbuffer(i,j)).gt.cutoff_hc)write(2,237)index_basis(icount,idorb(ishell)),index_basis(1,isorb(jshell)),dkdsbuffer(i,j),dndsbuffer(i,j)
 end do
 end do
 end do
 end do






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! done with (d|s) one electron integrals
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!    (d|p) integrals
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop over (d|p) batches
do ishell=1,num_p_shell
do jshell=1,num_d_shell
distsq=(shell_coor(1,iporb(ishell))-shell_coor(1,idorb(jshell)))**2+(shell_coor(2,iporb(ishell))-shell_coor(2,idorb(jshell) ))**2+(shell_coor(3,iporb(ishell))-shell_coor(3,idorb(jshell)))**2
  dpbuffer=0.
  dkdpbuffer=0.
  dndpbuffer=0.
 do mu=1,num_contr(iporb(ishell))
 do nu=1,num_contr(idorb(jshell))
 call buildp(shell_exp(mu,iporb(ishell)),shell_exp(nu,idorb(jshell)),iporb(ishell),idorb(jshell),pvector)

 call overlap_dp(distsq,shell_exp(mu,iporb(ishell)),shell_exp(nu,idorb(jshell)),iporb(ishell),idorb(jshell),dpints,pvector)

 call kinetic_dp(distsq,shell_exp(mu,iporb(ishell)),shell_exp(nu,idorb(jshell)),iporb(ishell),idorb(jshell),dpints,pvector &
,dkdpints)

 dpbuffer=dpbuffer+dpints*shell_coefs(mu,iporb(ishell))*shell_coefs(nu,idorb(jshell))
 dkdpbuffer=dkdpbuffer+dkdpints*shell_coefs(mu,iporb(ishell))*shell_coefs(nu,idorb(jshell))

scrdpn=0.
do inuc=1,natoms
 call nuclear_dp(distsq,shell_exp(mu,iporb(ishell)),shell_exp(nu,idorb(jshell)),iporb(ishell),idorb(jshell),pvector,0,&
dndpints)
scrdpn=scrdpn+dndpints
end do
dndpbuffer=dndpbuffer+scrdpn*shell_coefs(mu,iporb(ishell))*shell_coefs(nu,idorb(jshell))

 end do
 end do

do kk=1,3
do mm=1,3
dpbuffer(mm,mm,kk)=dpbuffer(mm,mm,kk)*dorbscale
dkdpbuffer(mm,mm,kk)=dkdpbuffer(mm,mm,kk)*dorbscale
dndpbuffer(mm,mm,kk)=dndpbuffer(mm,mm,kk)*dorbscale
end do
end do




do kk=1,3
icount=0
do nn=1,3
do mm=1,nn
icount=icount+1
!print*,index_basis(kk,iporb(ishell)),index_basis(icount,idorb(jshell)),dpbuffer(mm,nn,kk)
!print*,'dp kin',index_basis(kk,iporb(ishell)),index_basis(icount,idorb(jshell)),dkdpbuffer(mm,nn,kk),dndpbuffer(mm,nn,kk)
write(1,237)index_basis(kk,iporb(ishell)),index_basis(icount,idorb(jshell)),dpbuffer(mm,nn,kk)
write(2,237)index_basis(kk,iporb(ishell)),index_basis(icount,idorb(jshell)),dkdpbuffer(mm,nn,kk),dndpbuffer(mm,nn,kk)
end do
end do
end do


end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   done with (d|p) integrals
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    (d|d) integrals
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! loop over (d|d) batches
do ishell=1,num_d_shell
do jshell=1,ishell
distsq=(shell_coor(1,idorb(ishell))-shell_coor(1,idorb(jshell)))**2+(shell_coor(2,idorb(ishell))-&
shell_coor(2,idorb(jshell) ))**2+(shell_coor(3,idorb(ishell))-shell_coor(3,idorb(jshell)))**2




  ddbuffer=0.
  dkddbuffer=0.
  dnddbuffer=0.
 do mu=1,num_contr(idorb(ishell))
 do nu=1,num_contr(idorb(jshell))
 call buildp(shell_exp(mu,idorb(ishell)),shell_exp(nu,idorb(jshell)),idorb(ishell),idorb(jshell),pvector)

 call overlap_dd(distsq,shell_exp(mu,idorb(ishell)),shell_exp(nu,idorb(jshell)),idorb(ishell),idorb(jshell),ddints,pvector)

 call kinetic_dd(distsq,shell_exp(mu,idorb(ishell)),shell_exp(nu,idorb(jshell)),idorb(ishell),idorb(jshell),ddints,&
pvector,dkddints)

 ddbuffer=ddbuffer+ddints*shell_coefs(mu,idorb(ishell))*shell_coefs(nu,idorb(jshell))
 dkddbuffer=dkddbuffer+dkddints*shell_coefs(mu,idorb(ishell))*shell_coefs(nu,idorb(jshell))


scrddn=0.
do inuc=1,natoms
 call nuclear_dd(distsq,shell_exp(mu,idorb(ishell)),shell_exp(nu,idorb(jshell)),idorb(ishell),idorb(jshell),pvector,0,&
dnddints)
scrddn=scrddn+dnddints
end do
dnddbuffer=dnddbuffer+scrddn*shell_coefs(mu,idorb(ishell))*shell_coefs(nu,idorb(jshell))


 end do
 end do

do kk=1,3
do mm=1,3
do ll=1,3
ddbuffer(ll,ll,mm,kk)=ddbuffer(ll,ll,mm,kk)*dorbscale
dkddbuffer(ll,ll,mm,kk)=dkddbuffer(ll,ll,mm,kk)*dorbscale
dnddbuffer(ll,ll,mm,kk)=dnddbuffer(ll,ll,mm,kk)*dorbscale
ddbuffer(mm,kk,ll,ll)=ddbuffer(mm,kk,ll,ll)*dorbscale
dkddbuffer(mm,kk,ll,ll)=dkddbuffer(mm,kk,ll,ll)*dorbscale
dnddbuffer(mm,kk,ll,ll)=dnddbuffer(mm,kk,ll,ll)*dorbscale
end do
end do
end do




if(ishell.ne.jshell)then
do kl=1,6
do ij=1,6
ic=imapdl(ij)
ic2=imapdr(ij)
ic3=imapdl(kl)
ic4=imapdr(kl)
!print*,'dd kin',index_basis(kk,idorb(ishell)),index_basis(ll,idorb(jshell)),dkddbuffer(ic,ic2,ic3,ic4),dnddbuffer(ic,ic2,ic3,ic4)

if(dabs(ddbuffer(ic,ic2,ic3,ic4)).gt.cutoff_ov)write(1,237)index_basis(ij,idorb(ishell)),index_basis(kl,idorb(jshell)),&
ddbuffer(ic,ic2,ic3,ic4)

if(dabs(dkddbuffer(ic,ic2,ic3,ic4)+dnddbuffer(ic,ic2,ic3,ic4)).gt.cutoff_hc)write(2,237)index_basis(ij,idorb(ishell)),&
index_basis(kl,idorb(jshell)),dkddbuffer(ic,ic2,ic3,ic4),dnddbuffer(ic,ic2,ic3,ic4)

!if(index_basis(ij,idorb(ishell)).eq.16 .and. index_basis(kl,idorb(jshell)).eq.5)then
!print*,'<16|Z|5> =',dnddbuffer(ic,ic2,ic3,ic4)
!stop
!end if

end do 
end do


else

do kl=1,6
do ij=1,kl
ic=imapdl(ij)
ic2=imapdr(ij)
ic3=imapdl(kl)
ic4=imapdr(kl)
!print*,'dd kin',index_basis(kk,idorb(ishell)),index_basis(ll,idorb(jshell)),dkddbuffer(ic,ic2,ic3,ic4),dnddbuffer(ic,ic2,ic3\
!,ic4)

if(dabs(ddbuffer(ic,ic2,ic3,ic4)).gt.cutoff_ov)write(1,237)index_basis(ij,idorb(ishell)),index_basis(kl,idorb(jshell)),&
ddbuffer(ic,ic2,ic3,ic4)

if(dabs(dkddbuffer(ic,ic2,ic3,ic4)+dnddbuffer(ic,ic2,ic3,ic4)).gt.cutoff_hc)write(2,237)index_basis(ij,idorb(ishell)),&
index_basis(kl,idorb(jshell)),dkddbuffer(ic,ic2,ic3,ic4),dnddbuffer(ic,ic2,ic3,ic4)

 end do
 end do

end if


237 format(i6,i6,f20.15,f20.15)


end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   done with (d|d) integrals
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




if(myrank.eq.0)print*,'One electron integrals are done.'

end if !end for if(myrank.eq.0)above 

call mpi_barrier(mpi_comm_world,ierr)

! compute 4 index integrals

!print*,'ssss',myrank,iss_pairs_cpu(myrank+1),ifirst_on_cpu_ss(myrank+1),ilast_on_cpu_ss(myrank+1)
! loop over (ss|ss) batches
!allocate(iss_pairs_cpu(nprocs))
!allocate(ifirst_on_cpu_ss(nprocs))
!allocate(ilast_on_cpu_ss(nprocs))

iquad=iss_pairs_total*(iss_pairs_total+1)/2
a=float(iquad)/float(nprocs)
ieach=floor(a)
iquads_cpu=ieach
j=ieach
k=mod(iquad,nprocs)
if(k.ne.0)then
j=iquad-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
iquads_cpu(icount)=iquads_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
ifirst_quad_cpu(1)=1
ilast_quad_cpu(1)=ifirst_quad_cpu(1)+iquads_cpu(1)-1
iall_on_cpu(1)=ilast_quad_cpu(1)-ifirst_quad_cpu(1)+1
do i=2,nprocs
ifirst_quad_cpu(i)=ilast_quad_cpu(i-1)+1
ilast_quad_cpu(i)=ifirst_quad_cpu(i)+iquads_cpu(i)-1
iall_on_cpu(i)=ilast_quad_cpu(i)-ifirst_quad_cpu(i)+1
end do
if(myrank.eq.0)then
if(debug_ab2)then
print*,'ss|ss Distribution Among Processors'
print*,'Processor          First ss            Last ss         Total'
do i=1,nprocs
print*,i,'       ',ifirst_quad_cpu(i),'       ',ilast_quad_cpu(i),iall_on_cpu(i)
end do
end if
end if
call mpi_barrier(mpi_comm_world,ierr)

!!do ii=ifirst_quad_cpu(myrank+1),ilast_quad_cpu(myrank+1)

 ii=ifirst_quad_cpu(myrank+1)
 ifind=0
 do j=1,iss_pairs_total
 do i=1,j
 ifind=ifind+1
 if(ifind.eq.ii)then
 irowstart=i
 jcolstart=j
 goto 10
 end if
 end do
 end do
 10 continue

icurrent=0
do mmm=jcolstart,iss_pairs_total
kkk=irowstart
do nnn=kkk,mmm
icurrent=icurrent+1
if(icurrent.gt.iall_on_cpu(myrank+1))goto 140


iloop=nnn
kloop=mmm







ishell=iss_pairs(iloop,1)
jshell=iss_pairs(iloop,2)

distsqab=(shell_coor(1,jshell)-shell_coor(1,ishell))**2+(shell_coor(2,jshell)&
-shell_coor(2,ishell))**2+(shell_coor(3,jshell)-shell_coor(3,ishell))**2

!if(distsqab.gt.rdist)goto 2000

kshell=iss_pairs(kloop,1)
lshell=iss_pairs(kloop,2)


distsqcd=(shell_coor(1,kshell)-shell_coor(1,lshell))**2+(shell_coor(2,kshell)-shell_coor(2,lshell))**2+(shell_coor(3,kshell)-shell_coor(3,lshell))**2


sint=0.

!if(distsqcd.gt.rdist)goto 2000
do mu=1,num_contr(ishell)
do nu=1,num_contr(jshell)

call buildk(distsqab,shell_exp(mu,ishell),shell_exp(nu,jshell),ishell,jshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,ishell),shell_exp(nu,jshell),ishell,jshell,pvector)

do la=1,num_contr(kshell)
do isig=1,num_contr(lshell)

call buildk(distsqcd,shell_exp(la,kshell),shell_exp(isig,lshell),kshell,lshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
!call buildp(shell_exp(la,kshell),shell_exp(isig,lshell),kshell,lshell,qvector)
!zabcd=za_plus_zb+zc_plus_zd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,ishell)*shell_coefs(nu,jshell)*shell_coefs(la,kshell)*shell_coefs(isig,lshell)
if(dabs(dterm*kabkcd).gt.thresh)then
zabcd=za_plus_zb+zc_plus_zd
call buildp(shell_exp(mu,ishell),shell_exp(nu,jshell),ishell,jshell,pvector)
call buildp(shell_exp(la,kshell),shell_exp(isig,lshell),kshell,lshell,qvector)
call twoe_ssss(0,zabcd,za_plus_zb,zc_plus_zd,pvector,qvector,ssss)
sint=sint+ssss*dterm
end if
end do
end do
end do
end do

!print*,'ssss',index_basis(1,ishell),index_basis(1,jshell),index_basis(1,kshell),index_basis(1,lshell),sint
if(dabs(sint).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if
bfuncs(1,itotal)=index_basis(1,ishell)
bfuncs(2,itotal)=index_basis(1,jshell)
bfuncs(3,itotal)=index_basis(1,kshell)
bfuncs(4,itotal)=index_basis(1,lshell)
twoeints(itotal)=sint
end if
!write(3,*)index_basis(1,ishell),index_basis(1,jshell),index_basis(1,kshell),index_basis(1,lshell),sint
!write(30,*)index_basis(1,ishell),index_basis(1,jshell),index_basis(1,kshell),index_basis(1,lshell),sint
2000 continue
end do
irowstart=1
end do

140 continue
call mpi_barrier(mpi_comm_world,ierr)




if(ips_pairs_total.eq.0)then
itotalints_all_cpu=0
call mpi_reduce(itotal,itotalints_all_cpu,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
if(myrank.eq.0)print*,'There are ',itotalints_all_cpu,'integrals'
close(1)
close(2)
close(3)
call mpi_barrier(mpi_comm_world,ierr)
return
end if


if(myrank.eq.0)print*,'psss'
!twoeints=0.
!**********************************************************
! loop over (ps|ss) batches



kstart=ifirst_on_cpu_ps(myrank+1)-1

do iloop=1,ips_pairs_cpu(myrank+1)

kstart=kstart+1
ipshell=ips_pairs(kstart,1)
jsshell=ips_pairs(kstart,2)



distsqab=(shell_coor(1,jsshell)-shell_coor(1,ipshell))**2+(shell_coor(2,jsshell)&
-shell_coor(2,ipshell))**2+(shell_coor(3,jsshell)-shell_coor(3,ipshell))**2

!if(distsqab.gt.rdist)goto 2001

do kloop=1,iss_pairs_total
ksshell=iss_pairs(kloop,1)
lsshell=iss_pairs(kloop,2)

distsqcd=(shell_coor(1,ksshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,ksshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,ksshell)-shell_coor(3,lsshell))**2
!if(distsqcd.gt.rdist)goto 2001

psss_buffer=0.
iwrite=.false.
do mu=1,num_contr(ipshell)
do nu=1,num_contr(jsshell)

call buildk(distsqab,shell_exp(mu,ipshell),shell_exp(nu,jsshell),ipshell,jsshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,ipshell),shell_exp(nu,jsshell),ipshell,jsshell,pvector)

do la=1,num_contr(ksshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
!call buildp(shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,qvector)
!zabcd=za_plus_zb+zc_plus_zd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,ipshell)*shell_coefs(nu,jsshell)*shell_coefs(la,ksshell)*shell_coefs(isig,lsshell)

if(dabs(dterm*kabkcd).gt.thresh)then
iwrite=.true.
zabcd=za_plus_zb+zc_plus_zd
call buildp(shell_exp(mu,ipshell),shell_exp(nu,jsshell),ipshell,jsshell,pvector)
call buildp(shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,qvector)
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd
call twoe_psss(0,zabcd,za_plus_zb,zc_plus_zd,ipshell,pvector,qvector,wvector,psssints)
psss_buffer=psss_buffer+psssints*dterm
end if
end do
end do
end do
end do

if(iwrite)then
 do nn=1,3
! print*,'psss',index_basis(nn,iporb(ishell)),index_basis(1,isorb(jshell)),index_basis(1,isorb(kshell)),index_basis(1,isorb(lshell)),psss_buffer(nn)
!if(dabs(psss_buffer(nn)).gt.cutoff_2e)write(3,*)index_basis(nn,ipshell),index_basis(1,jsshell),index_basis(1,ksshell),index_basis(1,lsshell),psss_buffer(nn)
!write(31,*)index_basis(nn,ipshell),index_basis(1,jsshell),index_basis(1,ksshell),index_basis(1,lsshell),psss_buffer(nn)
if(dabs(psss_buffer(nn)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if
bfuncs(1,itotal)=index_basis(nn,ipshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(1,ksshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=psss_buffer(nn)
end if



end do
end if
2001 continue
end do
end do










if(myrank.eq.0)print*,'psps'
!twoeints=0.
!!**********************************************************************

!loop over psps batches

!**********************************************************

iquad=ips_pairs_total*(ips_pairs_total+1)/2
ifirst_quad_cpu=0
ilast_quad_cpu=0
iall_on_cpu=0
a=float(iquad)/float(nprocs)
ieach=floor(a)
iquads_cpu=ieach
j=ieach
k=mod(iquad,nprocs)

if(k.ne.0)then
j=iquad-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
iquads_cpu(icount)=iquads_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
ifirst_quad_cpu(1)=1
ilast_quad_cpu(1)=ifirst_quad_cpu(1)+iquads_cpu(1)-1
iall_on_cpu(1)=ilast_quad_cpu(1)-ifirst_quad_cpu(1)+1
do i=2,nprocs
ifirst_quad_cpu(i)=ilast_quad_cpu(i-1)+1
ilast_quad_cpu(i)=ifirst_quad_cpu(i)+iquads_cpu(i)-1
iall_on_cpu(i)=ilast_quad_cpu(i)-ifirst_quad_cpu(i)+1
end do
if(myrank.eq.0)then
if(debug_ab2)then
print*,'ps|ps Distribution Among Processors'
print*,'Processor          First ps            Last ps         Total'
do i=1,nprocs
print*,i,'       ',ifirst_quad_cpu(i),'       ',ilast_quad_cpu(i),iall_on_cpu(i)
end do
end if
end if
call mpi_barrier(mpi_comm_world,ierr)

 ii=ifirst_quad_cpu(myrank+1)
 ifind=0
 do j=1,ips_pairs_total
 do i=1,j
 ifind=ifind+1
 if(ifind.eq.ii)then
 irowstart=i
 jcolstart=j
 goto 80
 end if
 end do
 end do
 80 continue

icurrent=0
do mmm=jcolstart,ips_pairs_total
kkk=irowstart
do nnn=kkk,mmm
icurrent=icurrent+1
if(icurrent.gt.iall_on_cpu(myrank+1))goto 141


iloop=nnn
kloop=mmm





ipshell=ips_pairs(iloop,1)
jsshell=ips_pairs(iloop,2)


distsqab=(shell_coor(1,jsshell)-shell_coor(1,ipshell))**2+(shell_coor(2,jsshell)&
-shell_coor(2,ipshell))**2+(shell_coor(3,jsshell)-shell_coor(3,ipshell))**2

!if(distsqab.gt.rdist)goto 2002

kpshell=ips_pairs(kloop,1)
lsshell=ips_pairs(kloop,2)


distsqcd=(shell_coor(1,kpshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,kpshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,kpshell)-shell_coor(3,lsshell))**2

!if(distsqcd.gt.rdist)goto 2002

psps_buffer=0.
iwrite=.false.
do mu=1,num_contr(ipshell)
do nu=1,num_contr(jsshell)

call buildk(distsqab,shell_exp(mu,ipshell),shell_exp(nu,jsshell),ipshell,jsshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,ipshell),shell_exp(nu,jsshell),ipshell,jsshell,pvector)

do la=1,num_contr(kpshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
!call buildp(shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,qvector)
!zabcd=za_plus_zb+zc_plus_zd


kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,ipshell)*shell_coefs(nu,jsshell)*shell_coefs(la,kpshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
iwrite=.true.
zabcd=za_plus_zb+zc_plus_zd
call buildp(shell_exp(mu,ipshell),shell_exp(nu,jsshell),ipshell,jsshell,pvector)
call buildp(shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,qvector)
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd
call twoe_psps(0,zabcd,za_plus_zb,zc_plus_zd,ipshell,kpshell,pvector,qvector,wvector,pspsints)
psps_buffer=psps_buffer+pspsints*dterm
end if
end do
end do
end do
end do

if(iwrite)then

if(iloop.ne.kloop)then
 do nn=1,3
 do mm=1,3
!if(dabs(psps_buffer(mm,nn)).gt.cutoff_2e)write(3,*)index_basis(mm,ipshell),index_basis(1,jsshell),index_basis(nn,kpshell)&
!,index_basis(1,lsshell),psps_buffer(mm,nn)
!write(32,*)index_basis(mm,ipshell),index_basis(1,jsshell),index_basis(nn,kpshell)&
!,index_basis(1,lsshell),psps_buffer(mm,nn)

if(dabs(psps_buffer(mm,nn)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if
bfuncs(1,itotal)=index_basis(mm,ipshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(nn,kpshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=psps_buffer(mm,nn)
end if



end do
end do
else
 do nn=1,3
 do mm=1,nn
!if(dabs(psps_buffer(mm,nn)).gt.cutoff_2e)write(3,*)index_basis(mm,ipshell),index_basis(1,jsshell),index_basis(nn,kpshell)&
!,index_basis(1,lsshell),psps_buffer(mm,nn)
!write(32,*)index_basis(mm,ipshell),index_basis(1,jsshell),index_basis(nn,kpshell)&
!,index_basis(1,lsshell),psps_buffer(mm,nn)
if(dabs(psps_buffer(mm,nn)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if
bfuncs(1,itotal)=index_basis(mm,ipshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(nn,kpshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=psps_buffer(mm,nn)
end if



end do
end do
end if

end if !iwrite block

2002 continue
end do
irowstart=1
end do
141 continue




if(myrank.eq.0)print*,'ppss' 
!twoeints=0.
!!**********************************************************************

!loop over ppss batches

kstart=ifirst_on_cpu_pp(myrank+1)-1

do iloop=1,ipp_pairs_cpu(myrank+1)
kstart=kstart+1
ipshell=ipp_pairs(kstart,1)
jpshell=ipp_pairs(kstart,2)



distsqab=(shell_coor(1,jpshell)-shell_coor(1,ipshell))**2+(shell_coor(2,jpshell)&
-shell_coor(2,ipshell))**2+(shell_coor(3,jpshell)-shell_coor(3,ipshell))**2

!if(distsqab.gt.rdist)goto 2003

do kloop=1,iss_pairs_total
ksshell=iss_pairs(kloop,1)
lsshell=iss_pairs(kloop,2)


distsqcd=(shell_coor(1,ksshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,ksshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,ksshell)-shell_coor(3,lsshell))**2

!if(distsqcd.gt.rdist)goto 2003

ppss_buffer=0.
iwrite=.false.
do mu=1,num_contr(ipshell)
do nu=1,num_contr(jpshell)

call buildk(distsqab,shell_exp(mu,ipshell),shell_exp(nu,jpshell),ipshell,jpshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,ipshell),shell_exp(nu,jpshell),ipshell,jpshell,pvector)

do la=1,num_contr(ksshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
!call buildp(shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,qvector)
!zabcd=za_plus_zb+zc_plus_zd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,ipshell)*shell_coefs(nu,jpshell)*shell_coefs(la,ksshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
iwrite=.true.
zabcd=za_plus_zb+zc_plus_zd
call buildp(shell_exp(mu,ipshell),shell_exp(nu,jpshell),ipshell,jpshell,pvector)
call buildp(shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,qvector)
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd
call twoe_ppss(0,zabcd,za_plus_zb,zc_plus_zd,ipshell,jpshell,pvector,qvector,wvector,ppssints)
ppss_buffer=ppss_buffer+ppssints*dterm
end if
end do
end do
end do
end do

if(iwrite)then
if(ipshell.ne.jpshell)then
 do nn=1,3
 do mm=1,3
!if(dabs(ppss_buffer(mm,nn)).gt.cutoff_2e)write(3,*)index_basis(mm,ipshell),index_basis(nn,jpshell),&
!index_basis(1,ksshell),index_basis(1,lsshell),ppss_buffer(mm,nn)
!write(33,*)index_basis(mm,ipshell),index_basis(nn,jpshell),&
!index_basis(1,ksshell),index_basis(1,lsshell),ppss_buffer(mm,nn)
if(dabs(ppss_buffer(mm,nn)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(mm,ipshell)
bfuncs(2,itotal)=index_basis(nn,jpshell)
bfuncs(3,itotal)=index_basis(1,ksshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ppss_buffer(mm,nn)
end if



end do
end do
else

 do nn=1,3
 do mm=1,nn
!if(dabs(ppss_buffer(mm,nn)).gt.cutoff_2e)write(3,*)index_basis(mm,ipshell),index_basis(nn,jpshell),&
!index_basis(1,ksshell),index_basis(1,lsshell),ppss_buffer(mm,nn)
!write(33,*)index_basis(mm,ipshell),index_basis(nn,jpshell),&
!index_basis(1,ksshell),index_basis(1,lsshell),ppss_buffer(mm,nn)

if(dabs(ppss_buffer(mm,nn)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(mm,ipshell)
bfuncs(2,itotal)=index_basis(nn,jpshell)
bfuncs(3,itotal)=index_basis(1,ksshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ppss_buffer(mm,nn)
end if




end do
end do
end if

end if !iwrite
2003 continue
end do
end do










if(myrank.eq.0)print*,'ppps'
!twoeints=0.
!!**********************************************************************
!twoeints=0.
!loop over ppps batches
kstart=ifirst_on_cpu_pp(myrank+1)-1

do iloop=1,ipp_pairs_cpu(myrank+1)
kstart=kstart+1
ipshell=ipp_pairs(kstart,1)
jpshell=ipp_pairs(kstart,2)
distsqab=(shell_coor(1,jpshell)-shell_coor(1,ipshell))**2+(shell_coor(2,jpshell)&
-shell_coor(2,ipshell))**2+(shell_coor(3,jpshell)-shell_coor(3,ipshell))**2

!if(distsqab.gt.rdist)goto 2004

do kloop=1,ips_pairs_total
kpshell=ips_pairs(kloop,1)
lsshell=ips_pairs(kloop,2)


distsqcd=(shell_coor(1,kpshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,kpshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,kpshell)-shell_coor(3,lsshell))**2

!if(distsqcd.gt.rdist)goto 2004

ppps_buffer=0.
iwrite=.false.
do mu=1,num_contr(ipshell)
do nu=1,num_contr(jpshell)

call buildk(distsqab,shell_exp(mu,ipshell),shell_exp(nu,jpshell),ipshell,jpshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,ipshell),shell_exp(nu,jpshell),ipshell,jpshell,pvector)

do la=1,num_contr(kpshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
!call buildp(shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,qvector)
!zabcd=za_plus_zb+zc_plus_zd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,ipshell)*shell_coefs(nu,jpshell)*shell_coefs(la,kpshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
iwrite=.true.
zabcd=za_plus_zb+zc_plus_zd
call buildp(shell_exp(mu,ipshell),shell_exp(nu,jpshell),ipshell,jpshell,pvector)
call buildp(shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,qvector)
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_ppps(0,zabcd,za_plus_zb,zc_plus_zd,ipshell,jpshell,kpshell,pvector,qvector,wvector,pppsints)
ppps_buffer=ppps_buffer+pppsints*dterm
end if
end do
end do
end do
end do


if(iwrite)then
if(ipshell.ne.jpshell)then
do nn=1,3
do mm=1,3
 do ll=1,3
! print*,'ppps',index_basis(nn,iporb(ishell)),index_basis(mm,iporb(jshell)),index_basis(ll,iporb(kshell))&
!,index_basis(1,isorb(lshell)),ppps_buffer(ll,mm,nn)
!if(dabs(ppps_buffer(ll,mm,nn)).gt.cutoff_2e)write(3,*)&
! index_basis(ll,ipshell),index_basis(mm,jpshell),index_basis(nn,kpshell)&
!,index_basis(1,lsshell),ppps_buffer(ll,mm,nn)
!write(34,*) index_basis(ll,ipshell),index_basis(mm,jpshell),index_basis(nn,kpshell)&
!,index_basis(1,lsshell),ppps_buffer(ll,mm,nn)

if(dabs(ppps_buffer(ll,mm,nn)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ll,ipshell)
bfuncs(2,itotal)=index_basis(mm,jpshell)
bfuncs(3,itotal)=index_basis(nn,kpshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ppps_buffer(ll,mm,nn)
end if



end do
end do
end do

else
do ll=1,3
do nn=1,3
 do mm=1,nn
! print*,'ppps',index_basis(nn,iporb(ishell)),index_basis(mm,iporb(jshell)),index_basis(ll,iporb(kshell))&
!,index_basis(1,isorb(lshell)),ppps_buffer(ll,mm,nn)
!if(dabs(ppps_buffer(mm,nn,ll)).gt.cutoff_2e)write(3,*)&
! index_basis(mm,ipshell),index_basis(nn,jpshell),index_basis(ll,kpshell)&
!,index_basis(1,lsshell),ppps_buffer(mm,nn,ll)
!write(34,*) index_basis(mm,ipshell),index_basis(nn,jpshell),index_basis(ll,kpshell)&
!,index_basis(1,lsshell),ppps_buffer(mm,nn,ll)

if(dabs(ppps_buffer(mm,nn,ll)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(mm,ipshell)
bfuncs(2,itotal)=index_basis(nn,jpshell)
bfuncs(3,itotal)=index_basis(ll,kpshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ppps_buffer(mm,nn,ll)
end if



end do
end do
end do
end if

end if !iwrite
2004 continue
end do
end do





if(myrank.eq.0)print*,'pppp'
!loop over pppp batches
!twoeints=0.

iquad=ipp_pairs_total*(ipp_pairs_total+1)/2
ifirst_quad_cpu=0
ilast_quad_cpu=0
iall_on_cpu=0
a=float(iquad)/float(nprocs)
ieach=floor(a)
iquads_cpu=ieach
j=ieach
k=mod(iquad,nprocs)
if(k.ne.0)then
j=iquad-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
iquads_cpu(icount)=iquads_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
ifirst_quad_cpu(1)=1
ilast_quad_cpu(1)=ifirst_quad_cpu(1)+iquads_cpu(1)-1
iall_on_cpu(1)=ilast_quad_cpu(1)-ifirst_quad_cpu(1)+1
do i=2,nprocs
ifirst_quad_cpu(i)=ilast_quad_cpu(i-1)+1
ilast_quad_cpu(i)=ifirst_quad_cpu(i)+iquads_cpu(i)-1
iall_on_cpu(i)=ilast_quad_cpu(i)-ifirst_quad_cpu(i)+1
end do
if(myrank.eq.0)then
if(debug_ab2)then
print*,'pp|pp Distribution Among Processors'
print*,'Processor          First pp            Last pp         Total'
do i=1,nprocs
print*,i,'       ',ifirst_quad_cpu(i),'       ',ilast_quad_cpu(i),iall_on_cpu(i)
end do
end if
end if
call mpi_barrier(mpi_comm_world,ierr)

 ii=ifirst_quad_cpu(myrank+1)
 ifind=0
 do j=1,ipp_pairs_total
 do i=1,j
 ifind=ifind+1
 if(ifind.eq.ii)then
 irowstart=i
 jcolstart=j
 goto 81
 end if
 end do
 end do
 81 continue

icurrent=0
do mmm=jcolstart,ipp_pairs_total
kkk=irowstart
do nnn=kkk,mmm
icurrent=icurrent+1
if(icurrent.gt.iall_on_cpu(myrank+1))goto 142


iloop=nnn
kloop=mmm





ipshell=ipp_pairs(iloop,1)
jpshell=ipp_pairs(iloop,2)
distsqab=(shell_coor(1,jpshell)-shell_coor(1,ipshell))**2+(shell_coor(2,jpshell)&
-shell_coor(2,ipshell))**2+(shell_coor(3,jpshell)-shell_coor(3,ipshell))**2

!if(distsqab.gt.rdist)goto 2005
   
kpshell=ipp_pairs(kloop,1)
lpshell=ipp_pairs(kloop,2)

distsqcd=(shell_coor(1,kpshell)-shell_coor(1,lpshell))**2+(shell_coor(2,kpshell)&
-shell_coor(2,lpshell))**2+(shell_coor(3,kpshell)-shell_coor(3,lpshell))**2

!if(distsqcd.gt.rdist)goto 2005

pppp_buffer=0.
iwrite=.false.
do mu=1,num_contr(ipshell)
do nu=1,num_contr(jpshell)

call buildk(distsqab,shell_exp(mu,ipshell),shell_exp(nu,jpshell),ipshell,jpshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,ipshell),shell_exp(nu,jpshell),ipshell,jpshell,pvector)

do la=1,num_contr(kpshell)
do isig=1,num_contr(lpshell)

call buildk(distsqcd,shell_exp(la,kpshell),shell_exp(isig,lpshell),kpshell,lpshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
!call buildp(shell_exp(la,kpshell),shell_exp(isig,lpshell),kpshell,lpshell,qvector)
!zabcd=za_plus_zb+zc_plus_zd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,ipshell)*shell_coefs(nu,jpshell)*shell_coefs(la,kpshell)*shell_coefs(isig,lpshell)

if(dabs(dterm*kabkcd).gt.thresh)then
iwrite=.true.
zabcd=za_plus_zb+zc_plus_zd
call buildp(shell_exp(mu,ipshell),shell_exp(nu,jpshell),ipshell,jpshell,pvector)
call buildp(shell_exp(la,kpshell),shell_exp(isig,lpshell),kpshell,lpshell,qvector)
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_pppp(0,zabcd,za_plus_zb,zc_plus_zd,ipshell,jpshell,kpshell,lpshell,pvector,qvector,wvector,ppppints)
pppp_buffer=pppp_buffer+ppppints*dterm
end if
end do
end do
end do
end do


if(iwrite)then
irow=1
ipppp_write=0
 do nn=1,3
 do mm=1,3
 do ll=1,3
 do kk=1,3
indi=index_basis(kk,ipshell)
indj=index_basis(ll,jpshell)
indk=index_basis(mm,kpshell)
indl=index_basis(nn,lpshell)
CALL PACK(INDI,INDJ,IJ)
CALL PACK(INDK,INDL,KL)
CALL PACK2G(IJ,KL,JOFFSET)
do ii=1,irow
if(joffset.eq.ipppp_write(ii))goto 25
end do
ipppp_write(irow)=joffset
irow=irow+1
! if(dabs(pppp_buffer(kk,ll,mm,nn)).gt.cutoff_2e)write(3,*)index_basis(kk,ipshell),index_basis(ll,jpshell),&
! index_basis(mm,kpshell),index_basis(nn,lpshell),pppp_buffer(kk,ll,mm,nn)

if(dabs(pppp_buffer(kk,ll,mm,nn)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(kk,ipshell)
bfuncs(2,itotal)=index_basis(ll,jpshell)
bfuncs(3,itotal)=index_basis(mm,kpshell)
bfuncs(4,itotal)=index_basis(nn,lpshell)
twoeints(itotal)=pppp_buffer(kk,ll,mm,nn)
end if



!write(35,*)index_basis(kk,ipshell),index_basis(ll,jpshell),&
! index_basis(mm,kpshell),index_basis(nn,lpshell),pppp_buffer(kk,ll,mm,nn)

25 continue
 end do
 end do
 end do
 end do

end if !iwrite

2005 continue
end do
irowstart=1
end do

142 continue



!itotalints_all_cpu=0
!call mpi_reduce(itotal,itotalints_all_cpu,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
!if(myrank.eq.0)print*,'There are ',itotalints_all_cpu,'integrals'
!close(1)
!close(2)
!close(3)
!open(unit=12,file='twoedump')
!do kkk=1,itotal
!write(12,34)bfuncs(1,kkk),bfuncs(2,kkk),bfuncs(3,kkk),bfuncs(4,kkk),twoeints(kkk)
!end do
!close(12)
!call mpi_barrier(mpi_comm_world,ierr)
!return







if(num_d_shell.eq.0)then
itotalints_all_cpu=0
call mpi_reduce(itotal,itotalints_all_cpu,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
if(myrank.eq.0)print*,'There are ',itotalints_all_cpu,'integrals'
close(1)
close(2)
close(3)
call mpi_barrier(mpi_comm_world,ierr)
return
end if

if(myrank.eq.0)print*,'dsss'
! start of d integrals!!!!!!!!!!!!!!!!!!
!**********************************************************
! loop over (ds|ss) batches
kstart=ifirst_on_cpu_ds(myrank+1)-1

do iloop=1,ids_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=ids_pairs(kstart,1)
jsshell=ids_pairs(kstart,2)
distsqab=(shell_coor(1,jsshell)-shell_coor(1,idshell))**2+(shell_coor(2,jsshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jsshell)-shell_coor(3,idshell))**2
!if(distsqab.gt.rdist)goto 2006

do kloop=1,iss_pairs_total
ksshell=iss_pairs(kloop,1)
lsshell=iss_pairs(kloop,2)
distsqcd=(shell_coor(1,ksshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,ksshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,ksshell)-shell_coor(3,lsshell))**2
!if(distsqcd.gt.rdist)goto 2006

dsss_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jsshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)

do la=1,num_contr(ksshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
!call buildp(shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,qvector)
!zabcd=za_plus_zb+zc_plus_zd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jsshell)*shell_coefs(la,ksshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
zabcd=za_plus_zb+zc_plus_zd
call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)
call buildp(shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,qvector)
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd
call twoe_dsss(0,zabcd,za_plus_zb,zc_plus_zd,idshell,pvector,qvector,wvector,dsssints)
dsss_buffer=dsss_buffer+dsssints*dterm
end if
end do
end do
end do
end do


do i=1,3
dsss_buffer(i,i)=dsss_buffer(i,i)*dorbscale
end do

 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dsss',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(1,ksshell),index_basis(1,lsshell),dsss_buffer(i,j)
!write(20,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(1,ksshell),index_basis(1,lsshell),dsss_buffer(i,j)

!if(dabs(dsss_buffer(i,j)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(1,ksshell),index_basis(1,lsshell),dsss_buffer(i,j)

if(dabs(dsss_buffer(i,j)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(1,ksshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=dsss_buffer(i,j)
end if



 end do
 end do




2006 continue
end do
end do
 
if(myrank.eq.0)print*,'dpss',itotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (dp|ss)

kstart=ifirst_on_cpu_dp(myrank+1)-1

do iloop=1,idp_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=idp_pairs(kstart,1)
jpshell=idp_pairs(kstart,2)

distsqab=(shell_coor(1,jpshell)-shell_coor(1,idshell))**2+(shell_coor(2,jpshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jpshell)-shell_coor(3,idshell))**2

!if(distsqab.gt.rdist)goto 2007

do kloop=1,iss_pairs_total
ksshell=iss_pairs(kloop,1)
lsshell=iss_pairs(kloop,2)
distsqcd=(shell_coor(1,ksshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,ksshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,ksshell)-shell_coor(3,lsshell))**2

!if(distsqcd.gt.rdist)goto 2007

dpss_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jpshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,pvector)

do la=1,num_contr(ksshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
!call buildp(shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,qvector)
!zabcd=za_plus_zb+zc_plus_zd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jpshell)*shell_coefs(la,ksshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
zabcd=za_plus_zb+zc_plus_zd
call buildp(shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,pvector)
call buildp(shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,qvector)
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd
call twoe_dpss(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jpshell,pvector,qvector,wvector,dpssints)
dpss_buffer=dpss_buffer+dpssints*dterm
end if
end do
end do
end do
end do

do j=1,3
do i=1,3
dpss_buffer(i,i,j)=dpss_buffer(i,i,j)*dorbscale
end do
end do

 do k=1,3
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dpss',index_basis(icount,idshell),index_basis(k,jpshell),index_basis(1,ksshell),index_basis(1,lsshell),&
!dpss_buffer(i,j,k)
!write(21,*)index_basis(icount,idshell),index_basis(k,jpshell),index_basis(1,ksshell),index_basis(1,lsshell),&
!dpss_buffer(i,j,k)

!if(dabs(dpss_buffer(i,j,k)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(k,jpshell),index_basis(1,ksshell),index_basis(1,lsshell),&
!dpss_buffer(i,j,k)

if(dabs(dpss_buffer(i,j,k)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(k,jpshell)
bfuncs(3,itotal)=index_basis(1,ksshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=dpss_buffer(i,j,k)
end if



 end do
 end do
 end do

2007 continue
end do
end do

if(myrank.eq.0)print*,'dsps',itotal

!!! (ds|ps)
kstart=ifirst_on_cpu_ds(myrank+1)-1

do iloop=1,ids_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=ids_pairs(kstart,1)
jsshell=ids_pairs(kstart,2)

distsqab=(shell_coor(1,jsshell)-shell_coor(1,idshell))**2+(shell_coor(2,jsshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jsshell)-shell_coor(3,idshell))**2
!if(distsqab.gt.rdist)goto 2008

do kloop=1,ips_pairs_total
kpshell=ips_pairs(kloop,1)
lsshell=ips_pairs(kloop,2)
distsqcd=(shell_coor(1,kpshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,kpshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,kpshell)-shell_coor(3,lsshell))**2
!if(distsqcd.gt.rdist)goto 2008

dsps_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jsshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)

do la=1,num_contr(kpshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
!call buildp(shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,qvector)


kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jsshell)*shell_coefs(la,kpshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
zabcd=za_plus_zb+zc_plus_zd
call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)
call buildp(shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,qvector)
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_dsps(0,zabcd,za_plus_zb,zc_plus_zd,idshell,kpshell,pvector,qvector,wvector,dspsints)
dsps_buffer=dsps_buffer+dspsints*dterm
end if
end do
end do
end do
end do

do j=1,3
do i=1,3
dsps_buffer(i,i,j)=dsps_buffer(i,i,j)*dorbscale
end do
end do


 do k=1,3
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dsps',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(k,kpshell),index_basis(1,lsshell),&
!dsps_buffer(i,j,k)
!write(22,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(k,kpshell),index_basis(1,lsshell),&
!dsps_buffer(i,j,k)

!if(dabs(dsps_buffer(i,j,k)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(k,kpshell),index_basis(1,lsshell),&
!dsps_buffer(i,j,k)

if(dabs(dsps_buffer(i,j,k)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(k,kpshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=dsps_buffer(i,j,k)
end if

 end do
 end do
 end do

2008 continue
end do
end do


if(myrank.eq.0)print*,'dpps',itotal


!!!!!!!!!!!(dp|ps)

kstart=ifirst_on_cpu_dp(myrank+1)-1

do iloop=1,idp_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=idp_pairs(kstart,1)
jpshell=idp_pairs(kstart,2)

distsqab=(shell_coor(1,jpshell)-shell_coor(1,idshell))**2+(shell_coor(2,jpshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jpshell)-shell_coor(3,idshell))**2

!if(distsqab.gt.rdist)goto 2009

do kloop=1,ips_pairs_total
kpshell=ips_pairs(kloop,1)
lsshell=ips_pairs(kloop,2)
distsqcd=(shell_coor(1,kpshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,kpshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,kpshell)-shell_coor(3,lsshell))**2

!if(distsqcd.gt.rdist)goto 2009
dpps_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jpshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,pvector)

do la=1,num_contr(kpshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jpshell)*shell_coefs(la,kpshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,pvector)
call buildp(shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_dpps(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jpshell,kpshell,pvector,qvector,wvector,dppsints)
dpps_buffer=dpps_buffer+dppsints*dterm
end if
end do
end do
end do
end do

do l=1,3
do k=1,3
do i=1,3
dpps_buffer(i,i,k,l)=dpps_buffer(i,i,k,l)*dorbscale
end do
end do
end do

 do l=1,3
 do k=1,3
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dpps',index_basis(icount,idshell),index_basis(k,jpshell),index_basis(l,kpshell),index_basis(1,lsshell),&
!dpps_buffer(i,j,k,l)
!write(23,*)index_basis(icount,idshell),index_basis(k,jpshell),index_basis(l,kpshell),index_basis(1,lsshell),&
!dpps_buffer(i,j,k,l)

!if(dabs(dpps_buffer(i,j,k,l)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(k,jpshell),index_basis(l,kpshell),index_basis(1,lsshell),&
!dpps_buffer(i,j,k,l)
if(dabs(dpps_buffer(i,j,k,l)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(k,jpshell)
bfuncs(3,itotal)=index_basis(l,kpshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=dpps_buffer(i,j,k,l)
end if



 end do
 end do
 end do
 end do

2009 continue
end do
end do


if(myrank.eq.0)print*,'dspp',itotal

!!!!!!!!!!!(ds|pp)
kstart=ifirst_on_cpu_ds(myrank+1)-1

do iloop=1,ids_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=ids_pairs(kstart,1)
jsshell=ids_pairs(kstart,2)

distsqab=(shell_coor(1,jsshell)-shell_coor(1,idshell))**2+(shell_coor(2,jsshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jsshell)-shell_coor(3,idshell))**2

!if(distsqab.gt.rdist)goto 2010

do kloop=1,ipp_pairs_total
kpshell=ipp_pairs(kloop,1)
lpshell=ipp_pairs(kloop,2)
distsqcd=(shell_coor(1,kpshell)-shell_coor(1,lpshell))**2+&
(shell_coor(2,kpshell)-shell_coor(2,lpshell))**2&
+(shell_coor(3,kpshell)-shell_coor(3,lpshell))**2
!if(distsqcd.gt.rdist)goto 2010

dspp_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jsshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)

do la=1,num_contr(kpshell)
do isig=1,num_contr(lpshell)

call buildk(distsqcd,shell_exp(la,kpshell),shell_exp(isig,lpshell),kpshell,lpshell,dkcd,zc_plus_zd)
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jsshell)*shell_coefs(la,kpshell)*shell_coefs(isig,lpshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)
call buildp(shell_exp(la,kpshell),shell_exp(isig,lpshell),kpshell,lpshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_dspp(0,zabcd,za_plus_zb,zc_plus_zd,idshell,kpshell,lpshell,pvector,qvector,wvector,dsppints)
dspp_buffer=dspp_buffer+dsppints*dterm
end if
end do
end do
end do
end do

do l=1,3
do k=1,3
do i=1,3
dspp_buffer(i,i,k,l)=dspp_buffer(i,i,k,l)*dorbscale
end do
end do
end do


if(kpshell.ne.lpshell)then
 do l=1,3
 do k=1,3
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dspp',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(k,kpshell),index_basis(l,lpshell),&
!dspp_buffer(i,j,k,l)
!write(24,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(k,kpshell),index_basis(l,lpshell),&
!dspp_buffer(i,j,k,l)
!if(dabs(dspp_buffer(i,j,k,l)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(k,kpshell),index_basis(l,lpshell),&
!dspp_buffer(i,j,k,l)

if(dabs(dspp_buffer(i,j,k,l)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(k,kpshell)
bfuncs(4,itotal)=index_basis(l,lpshell)
twoeints(itotal)=dspp_buffer(i,j,k,l)
end if



 end do
 end do
 end do
 end do

else

 do l=1,3
 do k=1,l
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dspp',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(k,kpshell),index_basis(l,lpshell),&
!dspp_buffer(i,j,k,l)
!write(24,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(k,kpshell),index_basis(l,lpshell),&
!dspp_buffer(i,j,k,l)
!if(dabs(dspp_buffer(i,j,k,l)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(k,kpshell),index_basis(l,lpshell),&
!dspp_buffer(i,j,k,l)

if(dabs(dspp_buffer(i,j,k,l)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(k,kpshell)
bfuncs(4,itotal)=index_basis(l,lpshell)
twoeints(itotal)=dspp_buffer(i,j,k,l)
end if



 end do
 end do
 end do
 end do

end if



2010 continue
end do
end do



if(myrank.eq.0)print*,'dppp',itotal
!!!!!!!!!!!(dp|pp)
kstart=ifirst_on_cpu_dp(myrank+1)-1

do iloop=1,idp_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=idp_pairs(kstart,1)
jpshell=idp_pairs(kstart,2)

distsqab=(shell_coor(1,jpshell)-shell_coor(1,idshell))**2+(shell_coor(2,jpshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jpshell)-shell_coor(3,idshell))**2

!if(distsqab.gt.rdist)goto 2011

do kloop=1,ipp_pairs_total
kpshell=ipp_pairs(kloop,1)
lpshell=ipp_pairs(kloop,2)
distsqcd=(shell_coor(1,kpshell)-shell_coor(1,lpshell))**2+&
(shell_coor(2,kpshell)-shell_coor(2,lpshell))**2&
+(shell_coor(3,kpshell)-shell_coor(3,lpshell))**2

!if(distsqcd.gt.rdist)goto 2011

dppp_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jpshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,pvector)

do la=1,num_contr(kpshell)
do isig=1,num_contr(lpshell)

call buildk(distsqcd,shell_exp(la,kpshell),shell_exp(isig,lpshell),kpshell,lpshell,dkcd,zc_plus_zd)
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jpshell)*shell_coefs(la,kpshell)*shell_coefs(isig,lpshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,pvector)
call buildp(shell_exp(la,kpshell),shell_exp(isig,lpshell),kpshell,lpshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_dppp(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jpshell,kpshell,lpshell,pvector,qvector,wvector,dpppints)
dppp_buffer=dppp_buffer+dpppints*dterm
end if
end do
end do
end do
end do

do m=1,3
do l=1,3
do k=1,3
do i=1,3
dppp_buffer(i,i,k,l,m)=dppp_buffer(i,i,k,l,m)*dorbscale
end do
end do
end do
end do


if(kpshell.ne.lpshell)then
 do m=1,3
 do l=1,3
 do k=1,3
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dppp',index_basis(icount,idshell),index_basis(k,jpshell),index_basis(l,kpshell),index_basis(m,lpshell),&
!dppp_buffer(i,j,k,l,m)
!write(25,*)index_basis(icount,idshell),index_basis(k,jpshell),index_basis(l,kpshell),index_basis(m,lpshell),&
!dppp_buffer(i,j,k,l,m)
!if(dabs(dppp_buffer(i,j,k,l,m)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(k,jpshell),index_basis(l,kpshell),index_basis(m,lpshell),&
!dppp_buffer(i,j,k,l,m)

if(dabs(dppp_buffer(i,j,k,l,m)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(k,jpshell)
bfuncs(3,itotal)=index_basis(l,kpshell)
bfuncs(4,itotal)=index_basis(m,lpshell)
twoeints(itotal)=dppp_buffer(i,j,k,l,m)
end if



 end do
 end do
 end do
 end do
 end do


else

 do m=1,3
 do l=1,m
 do k=1,3
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dppp',index_basis(icount,idshell),index_basis(k,jpshell),index_basis(l,kpshell),index_basis(m,lpshell),&
!dppp_buffer(i,j,k,l,m)
!write(25,*)index_basis(icount,idshell),index_basis(k,jpshell),index_basis(l,kpshell),index_basis(m,lpshell),&
!dppp_buffer(i,j,k,l,m)
!if(dabs(dppp_buffer(i,j,k,l,m)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(k,jpshell),index_basis(l,kpshell),index_basis(m,lpshell),&
!dppp_buffer(i,j,k,l,m)

if(dabs(dppp_buffer(i,j,k,l,m)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(k,jpshell)
bfuncs(3,itotal)=index_basis(l,kpshell)
bfuncs(4,itotal)=index_basis(m,lpshell)
twoeints(itotal)=dppp_buffer(i,j,k,l,m)
end if



 end do
 end do
 end do
 end do
 end do

end if



2011 continue

end do
end do



if(myrank.eq.0)print*,'ddss',itotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (dd|ss)
kstart=ifirst_on_cpu_dd(myrank+1)-1

do iloop=1,idd_pairs_cpu(myrank+1)
kstart=kstart+1
idshell=idd_pairs(kstart,1)
jdshell=idd_pairs(kstart,2)

distsqab=(shell_coor(1,jdshell)-shell_coor(1,idshell))**2+(shell_coor(2,jdshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jdshell)-shell_coor(3,idshell))**2

!if(distsqab.gt.rdist)goto 2012

do kloop=1,iss_pairs_total
ksshell=iss_pairs(kloop,1)
lsshell=iss_pairs(kloop,2)
distsqcd=(shell_coor(1,ksshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,ksshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,ksshell)-shell_coor(3,lsshell))**2

!if(distsqcd.gt.rdist)goto 2012

ddss_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jdshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

do la=1,num_contr(ksshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jdshell)*shell_coefs(la,ksshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)
call buildp(shell_exp(la,ksshell),shell_exp(isig,lsshell),ksshell,lsshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_ddss(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jdshell,pvector,qvector,wvector,ddssints)
ddss_buffer=ddss_buffer+ddssints*dterm
end if
end do
end do
end do
end do

do l=1,3
do k=1,3
do i=1,3
ddss_buffer(i,i,k,l)=ddss_buffer(i,i,k,l)*dorbscale
ddss_buffer(k,l,i,i)=ddss_buffer(k,l,i,i)*dorbscale
end do
end do
end do


if(idshell.ne.jdshell)then
 do kl=1,6
 do ij=1,6
i=imapdl(ij)
j=imapdr(ij)
k=imapdl(kl)
l=imapdr(kl)
!print*,'ddss',index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(1,ksshell),index_basis(1,lsshell),&
!ddss_buffer(i,j,k,l)
!write(26,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(1,ksshell),index_basis(1,lsshell),&
!ddss_buffer(i,j,k,l)
!if(dabs(ddss_buffer(i,j,k,l)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(1,ksshell),index_basis(1,lsshell),&
!ddss_buffer(i,j,k,l)

if(dabs(ddss_buffer(i,j,k,l)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(1,ksshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ddss_buffer(i,j,k,l)
end if




 end do
 end do

else
 do kl=1,6
 do ij=1,kl
i=imapdl(ij)
j=imapdr(ij)
k=imapdl(kl)
l=imapdr(kl)
!print*,'ddss',index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(1,ksshell),index_basis(1,lsshell),&
!ddss_buffer(i,j,k,l)
!write(26,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(1,ksshell),index_basis(1,lsshell),&
!ddss_buffer(i,j,k,l)
!if(dabs(ddss_buffer(i,j,k,l)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(1,ksshell),index_basis(1,lsshell),&
!ddss_buffer(i,j,k,l)

if(dabs(ddss_buffer(i,j,k,l)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(1,ksshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ddss_buffer(i,j,k,l)
end if




 end do
 end do
end if


2012 continue
end do
end do



if(myrank.eq.0)print*,'dsds',itotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (ds|ds)
iquad=ids_pairs_total*(ids_pairs_total+1)/2
ifirst_quad_cpu=0
ilast_quad_cpu=0
iall_on_cpu=0
a=float(iquad)/float(nprocs)
ieach=floor(a)
iquads_cpu=ieach
j=ieach
k=mod(iquad,nprocs)
if(k.ne.0)then
j=iquad-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
iquads_cpu(icount)=iquads_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
ifirst_quad_cpu(1)=1
ilast_quad_cpu(1)=ifirst_quad_cpu(1)+iquads_cpu(1)-1
iall_on_cpu(1)=ilast_quad_cpu(1)-ifirst_quad_cpu(1)+1
do i=2,nprocs
ifirst_quad_cpu(i)=ilast_quad_cpu(i-1)+1
ilast_quad_cpu(i)=ifirst_quad_cpu(i)+iquads_cpu(i)-1
iall_on_cpu(i)=ilast_quad_cpu(i)-ifirst_quad_cpu(i)+1
end do
if(myrank.eq.0)then
if(debug_ab2)then
print*,'ds|ds Distribution Among Processors'
print*,'Processor          First ds            Last ds         Total'
do i=1,nprocs
print*,i,'       ',ifirst_quad_cpu(i),'       ',ilast_quad_cpu(i),iall_on_cpu(i)
end do
end if
end if
call mpi_barrier(mpi_comm_world,ierr)

!!do ii=ifirst_quad_cpu(myrank+1),ilast_quad_cpu(myrank+1)

 ii=ifirst_quad_cpu(myrank+1)
 ifind=0
 do j=1,ids_pairs_total
 do i=1,j
 ifind=ifind+1
 if(ifind.eq.ii)then
 irowstart=i
 jcolstart=j
 goto 31
 end if
 end do
 end do
 31 continue

icurrent=0
do mmm=jcolstart,ids_pairs_total
kkk=irowstart
do nnn=kkk,mmm
icurrent=icurrent+1
if(icurrent.gt.iall_on_cpu(myrank+1))goto 150


iloop=nnn
kloop=mmm

idshell=ids_pairs(iloop,1)
jsshell=ids_pairs(iloop,2)

distsqab=(shell_coor(1,jsshell)-shell_coor(1,idshell))**2+(shell_coor(2,jsshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jsshell)-shell_coor(3,idshell))**2

!if(distsqab.gt.rdist)goto 2013

kdshell=ids_pairs(kloop,1)
lsshell=ids_pairs(kloop,2)
distsqcd=(shell_coor(1,kdshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,kdshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,kdshell)-shell_coor(3,lsshell))**2

!if(distsqcd.gt.rdist)goto 2013

dsds_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jsshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)

do la=1,num_contr(kdshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,kdshell),shell_exp(isig,lsshell),kdshell,lsshell,dkcd,zc_plus_zd)
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jsshell)*shell_coefs(la,kdshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)

call buildp(shell_exp(la,kdshell),shell_exp(isig,lsshell),kdshell,lsshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_dsds(0,zabcd,za_plus_zb,zc_plus_zd,idshell,kdshell,pvector,qvector,wvector,dsdsints)
dsds_buffer=dsds_buffer+dsdsints*dterm
end if
end do
end do
end do
end do

do l=1,3
do k=1,3
do i=1,3
dsds_buffer(i,i,k,l)=dsds_buffer(i,i,k,l)*dorbscale
dsds_buffer(k,l,i,i)=dsds_buffer(k,l,i,i)*dorbscale
end do
end do
end do

if(iloop.ne.kloop)then
 kcount=0
 do l=1,3
 do k=1,l
 kcount=kcount+1
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(27,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!if(dabs(dsds_buffer(i,j,k,l)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)

if(dabs(dsds_buffer(i,j,k,l)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(kcount,kdshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=dsds_buffer(i,j,k,l)
end if



!write(29,*)dsds_buffer(i,j,k,l),dsds_buffer(j,i,k,l),dsds_buffer(i,j,l,k)
 end do
 end do
 end do
 end do

else

 do kl=1,6
 do ij=1,kl
i=imapdl(ij)
j=imapdr(ij)
k=imapdl(kl)
l=imapdr(kl)
!write(27,*)index_basis(ij,idshell),index_basis(1,jsshell),index_basis(kl,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!if(dabs(dsds_buffer(i,j,k,l)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(1,jsshell),index_basis(kl,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)


if(dabs(dsds_buffer(i,j,k,l)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(kl,kdshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=dsds_buffer(i,j,k,l)
end if



end do
end do

end if



2013 continue 
end do
irowstart=1
end do
150 continue

if(myrank.eq.0)print*,'ddps',itotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (dd|ps)
kstart=ifirst_on_cpu_dd(myrank+1)-1

do iloop=1,idd_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=idd_pairs(kstart,1)
jdshell=idd_pairs(kstart,2)

distsqab=(shell_coor(1,jdshell)-shell_coor(1,idshell))**2+(shell_coor(2,jdshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jdshell)-shell_coor(3,idshell))**2

do kloop=1,ips_pairs_total
kpshell=ips_pairs(kloop,1)
lsshell=ips_pairs(kloop,2)
distsqcd=(shell_coor(1,kpshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,kpshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,kpshell)-shell_coor(3,lsshell))**2

ddps_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jdshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,dkab,za_plus_zb)
if(distsqab.gt.rdist)goto 2014
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

do la=1,num_contr(kpshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,dkcd,zc_plus_zd)
if(distsqcd.gt.rdist)goto 2014
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jdshell)*shell_coefs(la,kpshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)
call buildp(shell_exp(la,kpshell),shell_exp(isig,lsshell),kpshell,lsshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_ddps(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jdshell,kpshell,pvector,qvector,wvector,ddpsints)
ddps_buffer=ddps_buffer+ddpsints*dterm
end if
end do
end do
end do
end do

do m=1,3
do l=1,3
do k=1,3
do i=1,3
ddps_buffer(i,i,k,l,m)=ddps_buffer(i,i,k,l,m)*dorbscale
ddps_buffer(k,l,i,i,m)=ddps_buffer(k,l,i,i,m)*dorbscale
end do
end do
end do
end do


if(idshell.ne.jdshell)then
 do m=1,3
 kcount=0
 do l=1,3
 do k=1,l
 kcount=kcount+1
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(28,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(m,kpshell),index_basis(1,lsshell),&
!ddps_buffer(i,j,k,l,m)
!if(dabs(ddps_buffer(i,j,k,l,m)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(m,kpshell),index_basis(1,lsshell),&
!ddps_buffer(i,j,k,l,m)

if(dabs(ddps_buffer(i,j,k,l,m)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(kcount,jdshell)
bfuncs(3,itotal)=index_basis(m,kpshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ddps_buffer(i,j,k,l,m)
end if



 end do
 end do
 end do
 end do
 end do

else
 do m=1,3 
 do kl=1,6
 do ij=1,kl
i=imapdl(ij)
j=imapdr(ij)
k=imapdl(kl)
l=imapdr(kl)
!write(28,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(m,kpshell),index_basis(1,lsshell),&
!ddps_buffer(i,j,k,l,m)
!if(dabs(ddps_buffer(i,j,k,l,m)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(m,kpshell),index_basis(1,lsshell),&
!ddps_buffer(i,j,k,l,m)


if(dabs(ddps_buffer(i,j,k,l,m)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(m,kpshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ddps_buffer(i,j,k,l,m)
end if




 end do
 end do
 end do

end if

2014 continue
end do
end do



if(myrank.eq.0)print*,'dsdp'

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (ds|dp)
kstart=ifirst_on_cpu_ds(myrank+1)-1

do iloop=1,ids_pairs_cpu(myrank+1)
kstart=kstart+1
idshell=ids_pairs(kstart,1)
jsshell=ids_pairs(kstart,2)

distsqab=(shell_coor(1,jsshell)-shell_coor(1,idshell))**2+(shell_coor(2,jsshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jsshell)-shell_coor(3,idshell))**2
!if(distsqab.gt.rdist)goto 2015

do kloop=1,idp_pairs_total
kdshell=idp_pairs(kloop,1)
lpshell=idp_pairs(kloop,2)
distsqcd=(shell_coor(1,kdshell)-shell_coor(1,lpshell))**2+&
(shell_coor(2,kdshell)-shell_coor(2,lpshell))**2&
+(shell_coor(3,kdshell)-shell_coor(3,lpshell))**2
!if(distsqcd.gt.rdist)goto 2015
dsdp_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jsshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)

do la=1,num_contr(kdshell)
do isig=1,num_contr(lpshell)

call buildk(distsqcd,shell_exp(la,kdshell),shell_exp(isig,lpshell),kdshell,lpshell,dkcd,zc_plus_zd)
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jsshell)*shell_coefs(la,kdshell)*shell_coefs(isig,lpshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jsshell),idshell,jsshell,pvector)
call buildp(shell_exp(la,kdshell),shell_exp(isig,lpshell),kdshell,lpshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_dsdp(0,zabcd,za_plus_zb,zc_plus_zd,idshell,kdshell,lpshell,pvector,qvector,wvector,dsdpints)
dsdp_buffer=dsdp_buffer+dsdpints*dterm
end if
end do
end do
end do
end do

do m=1,3
do l=1,3
do k=1,3
do i=1,3
dsdp_buffer(i,i,k,l,m)=dsdp_buffer(i,i,k,l,m)*dorbscale
dsdp_buffer(k,l,i,i,m)=dsdp_buffer(k,l,i,i,m)*dorbscale
end do
end do
end do
end do


 do m=1,3
 kcount=0
 do l=1,3
 do k=1,l
 kcount=kcount+1
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1

!write(29,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(m,lpshell),&
!dsdp_buffer(i,j,k,l,m)
!if(dabs(dsdp_buffer(i,j,k,l,m)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(m,lpshell),&
!dsdp_buffer(i,j,k,l,m)

if(dabs(dsdp_buffer(i,j,k,l,m)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(1,jsshell)
bfuncs(3,itotal)=index_basis(kcount,kdshell)
bfuncs(4,itotal)=index_basis(m,lpshell)
twoeints(itotal)=dsdp_buffer(i,j,k,l,m)
end if 


 end do
 end do
 end do
 end do
 end do

2015 continue
end do
end do




if(myrank.eq.0)print*,'ddpp',itotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (dd|pp)
kstart=ifirst_on_cpu_dd(myrank+1)-1

do iloop=1,idd_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=idd_pairs(kstart,1)
jdshell=idd_pairs(kstart,2)

distsqab=(shell_coor(1,jdshell)-shell_coor(1,idshell))**2+(shell_coor(2,jdshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jdshell)-shell_coor(3,idshell))**2

!if(distsqab.gt.rdist)goto 2016

do kloop=1,ipp_pairs_total
kpshell=ipp_pairs(kloop,1)
lpshell=ipp_pairs(kloop,2)
distsqcd=(shell_coor(1,kpshell)-shell_coor(1,lpshell))**2+&
(shell_coor(2,kpshell)-shell_coor(2,lpshell))**2&
+(shell_coor(3,kpshell)-shell_coor(3,lpshell))**2

!if(distsqcd.gt.rdist)goto 2016

ddpp_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jdshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

do la=1,num_contr(kpshell)
do isig=1,num_contr(lpshell)

call buildk(distsqcd,shell_exp(la,kpshell),shell_exp(isig,lpshell),kpshell,lpshell,dkcd,zc_plus_zd)
dkcd_global=dkcd
kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jdshell)*shell_coefs(la,kpshell)*shell_coefs(isig,lpshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

call buildp(shell_exp(la,kpshell),shell_exp(isig,lpshell),kpshell,lpshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_ddpp(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jdshell,kpshell,lpshell,pvector,qvector,wvector,ddppints)
ddpp_buffer=ddpp_buffer+ddppints*dterm
end if
end do
end do
end do
end do

do n=1,3
do m=1,3
do l=1,3
do k=1,3
do i=1,3
ddpp_buffer(i,i,k,l,m,n)=ddpp_buffer(i,i,k,l,m,n)*dorbscale
ddpp_buffer(k,l,i,i,m,n)=ddpp_buffer(k,l,i,i,m,n)*dorbscale
end do
end do
end do
end do
end do

if(idshell.ne.jdshell)then

 if(kpshell.ne.lpshell)then
 do n=1,3
 do m=1,3
 kcount=0
 do l=1,3
 do k=1,l
 kcount=kcount+1
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
 !print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
 !dsds_buffer(i,j,k,l)
 !write(40,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(m,kpshell),index_basis(n,lpshell),&
 !ddpp_buffer(i,j,k,l,m,n)
! if(dabs(ddpp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(m,kpshell),index_basis(n,lpshell),&
! ddpp_buffer(i,j,k,l,m,n)

 if(dabs(ddpp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(kcount,jdshell)
bfuncs(3,itotal)=index_basis(m,kpshell)
bfuncs(4,itotal)=index_basis(n,lpshell)
twoeints(itotal)=ddpp_buffer(i,j,k,l,m,n)
end if




 end do
 end do
 end do
 end do
 end do
 end do

                                    else

 do n=1,3
 do m=1,n
 kcount=0
 do l=1,3
 do k=1,l
 kcount=kcount+1
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
 !print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
 !dsds_buffer(i,j,k,l)
 !write(40,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(m,kpshell),index_basis(n,lpshell),&
 !ddpp_buffer(i,j,k,l,m,n)
!  if(dabs(ddpp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(m,kpshell),index_basis(n,lpshell),&
! ddpp_buffer(i,j,k,l,m,n)

 if(dabs(ddpp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(kcount,jdshell)
bfuncs(3,itotal)=index_basis(m,kpshell)
bfuncs(4,itotal)=index_basis(n,lpshell)
twoeints(itotal)=ddpp_buffer(i,j,k,l,m,n)
end if



 end do
 end do
 end do
 end do
 end do
 end do
 end if

else !if idshell equals jdshell part

 if(kpshell.ne.lpshell)then
 do n=1,3
 do m=1,3
 do kl=1,6
 do ij=1,kl
 i=imapdl(ij)
 j=imapdr(ij)
 k=imapdl(kl)
 l=imapdr(kl)
 !print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
 !dsds_buffer(i,j,k,l)
 !write(40,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(m,kpshell),index_basis(n,lpshell),&
 !ddpp_buffer(i,j,k,l,m,n)
!  if(dabs(ddpp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(m,kpshell),index_basis(n,lpshell),&
! ddpp_buffer(i,j,k,l,m,n)

 if(dabs(ddpp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(m,kpshell)
bfuncs(4,itotal)=index_basis(n,lpshell)
twoeints(itotal)=ddpp_buffer(i,j,k,l,m,n)
end if



 end do
 end do
 end do
 end do
                  else
 do n=1,3
 do m=1,n
 do kl=1,6
 do ij=1,kl
 i=imapdl(ij)
 j=imapdr(ij)
 k=imapdl(kl)
 l=imapdr(kl)
 !print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
 !dsds_buffer(i,j,k,l)
 !write(40,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(m,kpshell),index_basis(n,lpshell),&
 !ddpp_buffer(i,j,k,l,m,n)
!  if(dabs(ddpp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(m,kpshell),index_basis(n,lpshell),&
! ddpp_buffer(i,j,k,l,m,n)

 if(dabs(ddpp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(m,kpshell)
bfuncs(4,itotal)=index_basis(n,lpshell)
twoeints(itotal)=ddpp_buffer(i,j,k,l,m,n)
end if



 end do
 end do
 end do
 end do
              end if
end if !for whole block



2016 continue
end do
end do





if(myrank.eq.0)print*,'dpdp',itotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (dp|dp)
iquad=idp_pairs_total*(idp_pairs_total+1)/2
ifirst_quad_cpu=0
ilast_quad_cpu=0
iall_on_cpu=0
a=float(iquad)/float(nprocs)
ieach=floor(a)
iquads_cpu=ieach
j=ieach
k=mod(iquad,nprocs)
if(k.ne.0)then
j=iquad-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
iquads_cpu(icount)=iquads_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
ifirst_quad_cpu(1)=1
ilast_quad_cpu(1)=ifirst_quad_cpu(1)+iquads_cpu(1)-1
iall_on_cpu(1)=ilast_quad_cpu(1)-ifirst_quad_cpu(1)+1
do i=2,nprocs
ifirst_quad_cpu(i)=ilast_quad_cpu(i-1)+1
ilast_quad_cpu(i)=ifirst_quad_cpu(i)+iquads_cpu(i)-1
iall_on_cpu(i)=ilast_quad_cpu(i)-ifirst_quad_cpu(i)+1
end do
if(myrank.eq.0)then
if(debug_ab2)then
print*,'dp|dp Distribution Among Processors'
print*,'Processor          First dp            Last dp         Total'
do i=1,nprocs
print*,i,'       ',ifirst_quad_cpu(i),'       ',ilast_quad_cpu(i),iall_on_cpu(i)
end do
end if
end if
call mpi_barrier(mpi_comm_world,ierr)

!!do ii=ifirst_quad_cpu(myrank+1),ilast_quad_cpu(myrank+1)

 ii=ifirst_quad_cpu(myrank+1)
 ifind=0
 do j=1,idp_pairs_total
 do i=1,j
 ifind=ifind+1
 if(ifind.eq.ii)then
 irowstart=i
 jcolstart=j
 goto 32
 end if
 end do
 end do
 32 continue

icurrent=0
do mmm=jcolstart,idp_pairs_total
kkk=irowstart
do nnn=kkk,mmm
icurrent=icurrent+1
if(icurrent.gt.iall_on_cpu(myrank+1))goto 151

iloop=nnn
kloop=mmm

 


idshell=idp_pairs(iloop,1)
jpshell=idp_pairs(iloop,2)

distsqab=(shell_coor(1,jpshell)-shell_coor(1,idshell))**2+(shell_coor(2,jpshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jpshell)-shell_coor(3,idshell))**2
!if(distsqab.gt.rdist)goto 2017

kdshell=idp_pairs(kloop,1)
lpshell=idp_pairs(kloop,2)
distsqcd=(shell_coor(1,kdshell)-shell_coor(1,lpshell))**2+&
(shell_coor(2,kdshell)-shell_coor(2,lpshell))**2&
+(shell_coor(3,kdshell)-shell_coor(3,lpshell))**2
!if(distsqcd.gt.rdist)goto 2017

dpdp_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jpshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,pvector)

do la=1,num_contr(kdshell)
do isig=1,num_contr(lpshell)

call buildk(distsqcd,shell_exp(la,kdshell),shell_exp(isig,lpshell),kdshell,lpshell,dkcd,zc_plus_zd)
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jpshell)*shell_coefs(la,kdshell)*shell_coefs(isig,lpshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jpshell),idshell,jpshell,pvector)

call buildp(shell_exp(la,kdshell),shell_exp(isig,lpshell),kdshell,lpshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_dpdp(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jpshell,kdshell,lpshell,pvector,qvector,wvector,dpdpints)
dpdp_buffer=dpdp_buffer+dpdpints*dterm
end if
end do
end do
end do
end do


do n=1,3
do m=1,3
do l=1,3
do k=1,3
do i=1,3
dpdp_buffer(i,i,k,l,m,n)=dpdp_buffer(i,i,k,l,m,n)*dorbscale
dpdp_buffer(l,m,k,i,i,n)=dpdp_buffer(l,m,k,i,i,n)*dorbscale
end do
end do
end do
end do
end do

if(iloop .ne. kloop)then

 do n=1,3
 kcount=0
 do l=1,3
 do k=1,l
 kcount=kcount+1
 do m=1,3
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(41,*)index_basis(icount,idshell),index_basis(m,jpshell),index_basis(kcount,kdshell),index_basis(n,lpshell),&
!dpdp_buffer(i,j,m,k,l,n)
!if(dabs(dpdp_buffer(i,j,m,k,l,n)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(m,jpshell),index_basis(kcount,kdshell),index_basis(n,lpshell),&
!dpdp_buffer(i,j,m,k,l,n)

if(dabs(dpdp_buffer(i,j,m,k,l,n)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(m,jpshell)
bfuncs(3,itotal)=index_basis(kcount,kdshell)
bfuncs(4,itotal)=index_basis(n,lpshell)
twoeints(itotal)=dpdp_buffer(i,j,m,k,l,n)
end if



 end do
 end do
 end do
 end do
 end do
 end do

else


 
 do lm=1,6
 do n=1,3
 do ij=1,lm
 kmax=n 
 if(ij.ne.lm)kmax=3
 do k=1,kmax

!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
 i=imapdl(ij)
 j=imapdr(ij)
 l=imapdl(lm)
 m=imapdr(lm)

!write(41,*)index_basis(ij,idshell),index_basis(k,jpshell),index_basis(lm,kdshell),index_basis(n,lpshell),&
!dpdp_buffer(i,j,k,l,m,n)
!if(dabs(dpdp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(k,jpshell),index_basis(lm,kdshell),index_basis(n,lpshell),&
!dpdp_buffer(i,j,k,l,m,n)


if(dabs(dpdp_buffer(i,j,k,l,m,n)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(k,jpshell)
bfuncs(3,itotal)=index_basis(lm,kdshell)
bfuncs(4,itotal)=index_basis(n,lpshell)
twoeints(itotal)=dpdp_buffer(i,j,k,l,m,n)
end if


 end do
 end do
 end do
 end do







end if

2017 continue
end do
irowstart=1
end do
151 continue



if(myrank.eq.0)print*,'ddds',itotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (dd|ds)

kstart=ifirst_on_cpu_dd(myrank+1)-1

do iloop=1,idd_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=idd_pairs(kstart,1)
jdshell=idd_pairs(kstart,2)

distsqab=(shell_coor(1,jdshell)-shell_coor(1,idshell))**2+(shell_coor(2,jdshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jdshell)-shell_coor(3,idshell))**2

do kloop=1,ids_pairs_total
kdshell=ids_pairs(kloop,1)
lsshell=ids_pairs(kloop,2)
distsqcd=(shell_coor(1,kdshell)-shell_coor(1,lsshell))**2+&
(shell_coor(2,kdshell)-shell_coor(2,lsshell))**2&
+(shell_coor(3,kdshell)-shell_coor(3,lsshell))**2

ddds_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jdshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,dkab,za_plus_zb)
!if(distsqab.gt.rdist)goto 2018
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

do la=1,num_contr(kdshell)
do isig=1,num_contr(lsshell)

call buildk(distsqcd,shell_exp(la,kdshell),shell_exp(isig,lsshell),kdshell,lsshell,dkcd,zc_plus_zd)
!if(distsqcd.gt.rdist)goto 2018
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jdshell)*shell_coefs(la,kdshell)*shell_coefs(isig,lsshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

call buildp(shell_exp(la,kdshell),shell_exp(isig,lsshell),kdshell,lsshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_ddds(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jdshell,kdshell,pvector,qvector,wvector,dddsints)
ddds_buffer=ddds_buffer+dddsints*dterm
end if
end do
end do
end do
end do

do n=1,3
do m=1,3
do l=1,3
do k=1,3
do i=1,3
ddds_buffer(i,i,k,l,m,n)=ddds_buffer(i,i,k,l,m,n)*dorbscale
ddds_buffer(k,l,i,i,m,n)=ddds_buffer(k,l,i,i,m,n)*dorbscale
ddds_buffer(k,l,m,n,i,i)=ddds_buffer(k,l,m,n,i,i)*dorbscale
end do
end do
end do
end do
end do


if(idshell.ne.jdshell)then
 lcount=0
 do m=1,3
 do n=1,m
 lcount=lcount+1
 kcount=0
 do l=1,3
 do k=1,l
 kcount=kcount+1
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(42,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(lcount,kdshell),index_basis(1,lsshell),&
!ddds_buffer(i,j,k,l,n,m)
!if(dabs(ddds_buffer(i,j,k,l,n,m)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(lcount,kdshell),index_basis(1,lsshell),&
!ddds_buffer(i,j,k,l,n,m)

if(dabs(ddds_buffer(i,j,k,l,n,m)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(kcount,jdshell)
bfuncs(3,itotal)=index_basis(lcount,kdshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ddds_buffer(i,j,k,l,n,m)
end if



 end do
 end do
 end do
 end do
 end do
 end do

else
 lcount=0
 do m=1,3
 do n=1,m
 lcount=lcount+1
 do kl=1,6
 do ij=1,kl
 i=imapdl(ij)
 j=imapdr(ij)
 k=imapdl(kl)
 l=imapdr(kl)

!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(42,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(lcount,kdshell),index_basis(1,lsshell),&
!ddds_buffer(i,j,k,l,n,m)!!! check this if they wrong. i switch m,n to n,m
!if(dabs(ddds_buffer(i,j,k,l,n,m)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(lcount,kdshell),index_basis(1,lsshell),&
!ddds_buffer(i,j,k,l,n,m)

if(dabs(ddds_buffer(i,j,k,l,n,m)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(lcount,kdshell)
bfuncs(4,itotal)=index_basis(1,lsshell)
twoeints(itotal)=ddds_buffer(i,j,k,l,n,m)
end if



 end do
 end do
 end do
 end do
end if


2018 continue
end do
end do

if(myrank.eq.0)print*,'dddp',itotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (dd|dp)
kstart=ifirst_on_cpu_dd(myrank+1)-1

do iloop=1,idd_pairs_cpu(myrank+1)
kstart=kstart+1

idshell=idd_pairs(kstart,1)
jdshell=idd_pairs(kstart,2)

distsqab=(shell_coor(1,jdshell)-shell_coor(1,idshell))**2+(shell_coor(2,jdshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jdshell)-shell_coor(3,idshell))**2

!if(distsqab.gt.rdist)goto 2019

do kloop=1,idp_pairs_total
kdshell=idp_pairs(kloop,1)
lpshell=idp_pairs(kloop,2)
distsqcd=(shell_coor(1,kdshell)-shell_coor(1,lpshell))**2+&
(shell_coor(2,kdshell)-shell_coor(2,lpshell))**2&
+(shell_coor(3,kdshell)-shell_coor(3,lpshell))**2

!if(distsqcd.gt.rdist)goto 2019


dddp_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jdshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

do la=1,num_contr(kdshell)
do isig=1,num_contr(lpshell)

call buildk(distsqcd,shell_exp(la,kdshell),shell_exp(isig,lpshell),kdshell,lpshell,dkcd,zc_plus_zd)
dkcd_global=dkcd

kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jdshell)*shell_coefs(la,kdshell)*shell_coefs(isig,lpshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

call buildp(shell_exp(la,kdshell),shell_exp(isig,lpshell),kdshell,lpshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_dddp(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jdshell,kdshell,lpshell,pvector,qvector,wvector,dddpints)
dddp_buffer=dddp_buffer+dddpints*dterm
end if
end do
end do
end do
end do

do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do i=1,3
dddp_buffer(i,i,k,l,m,n,io)=dddp_buffer(i,i,k,l,m,n,io)*dorbscale
dddp_buffer(k,l,i,i,m,n,io)=dddp_buffer(k,l,i,i,m,n,io)*dorbscale
dddp_buffer(k,l,m,n,i,i,io)=dddp_buffer(k,l,m,n,i,i,io)*dorbscale
end do
end do
end do
end do
end do
end do


if(idshell.ne.jdshell)then
 do io=1,3
 lcount=0
 do m=1,3
 do n=1,m
 lcount=lcount+1
 kcount=0
 do l=1,3
 do k=1,l
 kcount=kcount+1
 icount=0
 do j=1,3
 do i=1,j
  icount=icount+1
!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(43,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(lcount,kdshell),index_basis(io,lpshell),&
!dddp_buffer(i,j,k,l,n,m,io)
!if(dabs(dddp_buffer(i,j,k,l,n,m,io)).gt.cutoff_2e)write(3,*)index_basis(icount,idshell),index_basis(kcount,jdshell),index_basis(lcount,kdshell),index_basis(io,lpshell),&
!dddp_buffer(i,j,k,l,n,m,io)

if(dabs(dddp_buffer(i,j,k,l,n,m,io)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(icount,idshell)
bfuncs(2,itotal)=index_basis(kcount,jdshell)
bfuncs(3,itotal)=index_basis(lcount,kdshell)
bfuncs(4,itotal)=index_basis(io,lpshell)
twoeints(itotal)=dddp_buffer(i,j,k,l,n,m,io)
end if


 end do
 end do
 end do
 end do
 end do
 end do
 end do

else
 do io=1,3
 lcount=0
 do m=1,3
 do n=1,m
 lcount=lcount+1
 do kl=1,6
 do ij=1,kl

 i=imapdl(ij)
 j=imapdr(ij)
 k=imapdl(kl)
 l=imapdr(kl)

!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(43,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(lcount,kdshell),index_basis(io,lpshell),&
!dddp_buffer(i,j,k,l,n,m,io)
!if(dabs(dddp_buffer(i,j,k,l,n,m,io)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(lcount,kdshell),index_basis(io,lpshell),&
!dddp_buffer(i,j,k,l,n,m,io)

if(dabs(dddp_buffer(i,j,k,l,n,m,io)).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(lcount,kdshell)
bfuncs(4,itotal)=index_basis(io,lpshell)
twoeints(itotal)=dddp_buffer(i,j,k,l,n,m,io)
end if


 end do
 end do
 end do
 end do
 end do




end if

2019 continue
end do
end do


if(myrank.eq.0)print*,'dddd',itotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! (dd|dd)
iquad=idd_pairs_total*(idd_pairs_total+1)/2
ifirst_quad_cpu=0
ilast_quad_cpu=0
iall_on_cpu=0
a=float(iquad)/float(nprocs)
ieach=floor(a)
iquads_cpu=ieach
j=ieach
k=mod(iquad,nprocs)
if(k.ne.0)then
j=iquad-(nprocs)*ieach
icount=0
do ll=1,j
icount=icount+1
iquads_cpu(icount)=iquads_cpu(icount)+1
if(icount.eq.nprocs)icount=0
end do
end if
ifirst_quad_cpu(1)=1
ilast_quad_cpu(1)=ifirst_quad_cpu(1)+iquads_cpu(1)-1
iall_on_cpu(1)=ilast_quad_cpu(1)-ifirst_quad_cpu(1)+1
do i=2,nprocs
ifirst_quad_cpu(i)=ilast_quad_cpu(i-1)+1
ilast_quad_cpu(i)=ifirst_quad_cpu(i)+iquads_cpu(i)-1
iall_on_cpu(i)=ilast_quad_cpu(i)-ifirst_quad_cpu(i)+1
end do
if(myrank.eq.0)then
if(debug_ab2)then
print*,'dd|dd Distribution Among Processors'
print*,'Processor          First dd            Last dd         Total'
do i=1,nprocs
print*,i,'       ',ifirst_quad_cpu(i),'       ',ilast_quad_cpu(i),iall_on_cpu(i)
end do
end if
end if
call mpi_barrier(mpi_comm_world,ierr)

!!do ii=ifirst_quad_cpu(myrank+1),ilast_quad_cpu(myrank+1)

 ii=ifirst_quad_cpu(myrank+1)
 ifind=0
 do j=1,idd_pairs_total
 do i=1,j
 ifind=ifind+1
 if(ifind.eq.ii)then
 irowstart=i
 jcolstart=j
 goto 33
 end if
 end do
 end do
 33 continue

icurrent=0
do mmm=jcolstart,idd_pairs_total
kkk=irowstart
do nnn=kkk,mmm
icurrent=icurrent+1
if(icurrent.gt.iall_on_cpu(myrank+1))goto 152

iloop=nnn
kloop=mmm



idshell=idd_pairs(iloop,1)
jdshell=idd_pairs(iloop,2)

distsqab=(shell_coor(1,jdshell)-shell_coor(1,idshell))**2+(shell_coor(2,jdshell)&
-shell_coor(2,idshell))**2+(shell_coor(3,jdshell)-shell_coor(3,idshell))**2

!if(distsqab.gt.rdist)goto 2020

kdshell=idd_pairs(kloop,1)
ldshell=idd_pairs(kloop,2)
distsqcd=(shell_coor(1,kdshell)-shell_coor(1,ldshell))**2+&
(shell_coor(2,kdshell)-shell_coor(2,ldshell))**2&
+(shell_coor(3,kdshell)-shell_coor(3,ldshell))**2

!if(distsqcd.gt.rdist)goto 2020

dddd_buffer=0.
do mu=1,num_contr(idshell)
do nu=1,num_contr(jdshell)

call buildk(distsqab,shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,dkab,za_plus_zb)
dkab_global=dkab
!call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

do la=1,num_contr(kdshell)
do isig=1,num_contr(ldshell)

call buildk(distsqcd,shell_exp(la,kdshell),shell_exp(isig,ldshell),kdshell,ldshell,dkcd,zc_plus_zd)
dkcd_global=dkcd


kabkcd=dkab_global*dkcd_global
dterm=shell_coefs(mu,idshell)*shell_coefs(nu,jdshell)*shell_coefs(la,kdshell)*shell_coefs(isig,ldshell)
if(dabs(dterm*kabkcd).gt.thresh)then
call buildp(shell_exp(mu,idshell),shell_exp(nu,jdshell),idshell,jdshell,pvector)

call buildp(shell_exp(la,kdshell),shell_exp(isig,ldshell),kdshell,ldshell,qvector)
zabcd=za_plus_zb+zc_plus_zd
wvector=za_plus_zb*pvector+zc_plus_zd*qvector
wvector=wvector/zabcd

call twoe_dddd(0,zabcd,za_plus_zb,zc_plus_zd,idshell,jdshell,kdshell,ldshell,pvector,qvector,wvector,ddddints)
dddd_buffer=dddd_buffer+ddddints*dterm
end if
end do
end do
end do
end do

do ip=1,3
do io=1,3
do n=1,3
do m=1,3
do l=1,3
do k=1,3
do i=1,3
dddd_buffer(i,i,k,l,m,n,isquish(io,ip))=dddd_buffer(i,i,k,l,m,n,isquish(io,ip))*dorbscale
dddd_buffer(k,l,i,i,m,n,isquish(io,ip))=dddd_buffer(k,l,i,i,m,n,isquish(io,ip))*dorbscale
dddd_buffer(k,l,m,n,i,i,isquish(io,ip))=dddd_buffer(k,l,m,n,i,i,isquish(io,ip))*dorbscale
dddd_buffer(k,l,m,n,io,ip,isquish(i,i))=dddd_buffer(k,l,m,n,io,ip,isquish(i,i))*dorbscale
end do
end do
end do
end do
end do
end do
end do

if(iloop.ne.kloop)then

     do iop=1,6
     mnmax=6
     if(kdshell.eq.ldshell)mnmax=iop
     do mn=1,mnmax
     do kl=1,6
     ijmax=6
     if(idshell.eq.jdshell)ijmax=kl
     do ij=1,ijmax

  i=imapdl(ij)
  j=imapdr(ij)
  k=imapdl(kl)
  l=imapdr(kl)
  m=imapdl(mn)
  n=imapdr(mn)
  io=imapdl(iop)
  ip=imapdr(iop)
!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(44,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(mn,kdshell),index_basis(iop,ldshell),&
!dddd_buffer(i,j,k,l,m,n,io,ip)
!if(dabs(dddd_buffer(i,j,k,l,m,n,io,ip)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(mn,kdshell),index_basis(iop,ldshell),&
!dddd_buffer(i,j,k,l,m,n,io,ip)


if(dabs(dddd_buffer(i,j,k,l,m,n,isquish(io,ip))).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(mn,kdshell)
bfuncs(4,itotal)=index_basis(iop,ldshell)
twoeints(itotal)=dddd_buffer(i,j,k,l,m,n,isquish(io,ip))
!!write(90,*)bfuncs(1,itotal),bfuncs(2,itotal),bfuncs(3,itotal),bfuncs(4,itotal),dddd_buffer(i,j,k,l,m,n,isquish(io,ip))
end if



    end do
    end do
    end do
    end do

                   else ! if iloop equal kloop

   if(idshell.ne.jdshell)then !also mean kdshell not equal ldshell 

     do iop=1,6
     do mn=1,6
     do kl=1,iop
     ijmax=mn
     if(kl.ne.iop)ijmax=6
     do ij=1,ijmax
  i=imapdl(ij)
  j=imapdr(ij)
  k=imapdl(kl)
  l=imapdr(kl)
  m=imapdl(mn)
  n=imapdr(mn)
  io=imapdl(iop)
  ip=imapdr(iop)
!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(44,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(mn,kdshell),index_basis(iop,ldshell),&
!dddd_buffer(i,j,k,l,m,n,io,ip)
!if(dabs(dddd_buffer(i,j,k,l,m,n,io,ip)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(mn,kdshell),index_basis(iop,ldshell),&
!dddd_buffer(i,j,k,l,m,n,io,ip)

if(dabs(dddd_buffer(i,j,k,l,m,n,isquish(io,ip))).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(mn,kdshell)
bfuncs(4,itotal)=index_basis(iop,ldshell)
twoeints(itotal)=dddd_buffer(i,j,k,l,m,n,isquish(io,ip))
!!write(90,*)bfuncs(1,itotal),bfuncs(2,itotal),bfuncs(3,itotal),bfuncs(4,itotal),dddd_buffer(i,j,k,l,m,n,isquish(io,ip))
end if



     end do
     end do
     end do
     end do

                           else

 irow=1
 ipppp_write=0
 do iop=1,6
 do mn=1,6
 do kl=1,6
 do ij=1,6
CALL PACK(ij,kl,ileft)
CALL PACK(mn,iop,iright)
CALL PACK2G(ileft,iright,joffset)
do ii=1,irow
if(joffset.eq.ipppp_write(ii))goto 30
end do
ipppp_write(irow)=joffset
irow=irow+1
  i=imapdl(ij)
  j=imapdr(ij)
  k=imapdl(kl)
  l=imapdr(kl)
  m=imapdl(mn)
  n=imapdr(mn)
  io=imapdl(iop)
  ip=imapdr(iop)
!print*,'dsds',index_basis(icount,idshell),index_basis(1,jsshell),index_basis(kcount,kdshell),index_basis(1,lsshell),&
!dsds_buffer(i,j,k,l)
!write(44,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(mn,kdshell),index_basis(iop,ldshell),&
!dddd_buffer(i,j,k,l,m,n,io,ip)
!if(dabs(dddd_buffer(i,j,k,l,m,n,io,ip)).gt.cutoff_2e)write(3,*)index_basis(ij,idshell),index_basis(kl,jdshell),index_basis(mn,kdshell),index_basis(iop,ldshell),&
!dddd_buffer(i,j,k,l,m,n,io,ip)

if(dabs(dddd_buffer(i,j,k,l,m,n,isquish(io,ip))).gt.cutoff_2e)then
itotal=itotal+1
if(itotal.gt.jmaxints*1000)then
print*,'Two electron integral buffer exceeded on processor',myrank
stop
end if

bfuncs(1,itotal)=index_basis(ij,idshell)
bfuncs(2,itotal)=index_basis(kl,jdshell)
bfuncs(3,itotal)=index_basis(mn,kdshell)
bfuncs(4,itotal)=index_basis(iop,ldshell)
twoeints(itotal)=dddd_buffer(i,j,k,l,m,n,isquish(io,ip))
end if





30 continue
 end do
 end do
 end do
 end do





     end if


end if

2020 continue
end do
irowstart=1
end do

152 continue



!print*,'itotal here',itotal




close(1)
close(2)
close(3)








itotalints_all_cpu=0
call mpi_reduce(itotal,itotalints_all_cpu,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
if(myrank.eq.0)print*,'There are ',itotalints_all_cpu,'integrals'
call mpi_barrier(mpi_comm_world,ierr)
 
call cpusec(timeout)
if(myrank.eq.0)write(*,200)timeout-timein

200 format(' Done!  Integral computation time (sec) =  ',F10.3)

return


print*,''
print*,'Processing two-electron integrals before SCF....'
open(unit=3,file='two_eri')
open(unit=5,file='iiii',form='unformatted',access='sequential')
40 continue
inds=0
rbuffer=0
do i=1,600
read(3,*,iostat=io)inds(1,i),inds(2,i),inds(3,i),inds(4,i),rbuffer(i)
 if(io<0)exit
end do
write(5)rbuffer,inds
 if(io<0)then
print*,'Done!'
print*,''
close(5)
 return
 else
 goto 40
 end if






return
end

