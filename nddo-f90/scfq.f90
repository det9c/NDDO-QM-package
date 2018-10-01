subroutine scfq(q0,natoms,ff,stat1,stat2,gout,numint)
use tables
use control
use indices
use constants
implicit double precision(a-h,o-z)
interface

subroutine find
end subroutine find

subroutine geometry
end subroutine geometry
end interface
integer,intent(in)::natoms,numint      
double precision,intent(inout) ::ff
double precision, dimension(numint),intent(inout) ::gout
double precision,dimension(natoms,3),intent(in)::q0
logical::stat1,stat2
!size(x0)
q=q0
if(OPTIMIZE .and. .not. TS  .and. .not. MODES)then
write(*,*)'GEOMETRY AT CURRENT ITERATE'
     do i=1,ninput
write(*,301)zstore(i),q(i,1),opt(i,1),q(i,2)*180.0d0/pi,opt(i,2),q(i,3)*180.0d0/pi,opt(i,3),ref(i,1),ref(i,2),ref(i,3)
     end do
301 format(i3,f10.5,i3,f10.5,i3,f10.5,i3,i3,i3,i3)
end if

!TRIAL=status
!keep=stat2

call geometry
call buildb
call nodummy



if(allocated(ifirst))deallocate(ifirst)
if(allocated(ifirst2))deallocate(ifirst2)
if(allocated(ilast))deallocate(ilast)
if(allocated(ilast2))deallocate(ilast2)
call two_electron
call hcore



TRIAL=stat1
keep=stat2
!DIIS=.FALSE.
if(restricted)then
call rhf
else
call uhf
end if
TRIAL=.true.
keep=.false.


if(restricted)then
call find
else
call finduhf
end if
gout=gint
ff=etotal

!call find
return 
end


