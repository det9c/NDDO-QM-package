subroutine default(ikey)
! sets keyword defaults
use tables
use control
interface
subroutine atoi(input,j)
character,intent(in)::input
integer,intent(inout)::j
end subroutine atoi

subroutine char_to_float(word,value,length)
      character(200),intent(in)::word
      integer,intent(in)::length
      double precision,intent(inout)::value
end subroutine char_to_float


end interface
integer,intent(in)::ikey
integer::iexp
character(200)::ch
double precision::valuedp

scftol=1D-6
!METHOD='AM1'
abinitio=.true.
DEORTHO=.false.
GRADIENT=.false.
TRIAL=.false.
XYZ=.false.
OPTIMIZE=.false.
OPTTOL=1.0d0
parameters=.false.
core='standard'
CENTERS=.false.
newton=.false.
restricted=.true.
MULT=1
sys_charge=0.0d0
DFTHF=.false.
LAPACK=.FALSE.
SPARKLES=.FALSE.
densityin=.false.
densityout=.false.
FREQUENCY=.FALSE.
TS=.FALSE.
DIIS=.FALSE.
MODES=.FALSE.
CALCALL=.FALSE.
UPDATE=200
IRC=.FALSE.
IRCSTPS=50
EULER=.FALSE.
RK4=.FALSE.
FORWARD=.FALSE.
SAM=.FALSE.
DEL1=1D-3
STE1=.1
ERIBUFF=500
CUTOFF=15.0D0
cutoff_sq=cutoff*cutoff
debug=.false.
CGSCF=.false.
sdtocg=0.5
pbc=.false.


do i=1,ikey
ch=' '
ch=keyword(i)




  if(ch=='AM1')then
  METHOD='AM1'
  core='GAUSSIAN'
  abinitio=.false.
  elseif(ch=='PM3')then
  METHOD='PM3'
  core='GAUSSIAN'
   abinitio=.false.
  elseif(ch=='MNDO')then
  METHOD='MNDO'
   abinitio=.false.
  elseif(ch=='MNDOD')then
  METHOD='MNDOD'
   abinitio=.false.
  elseif(ch(1:4)=='CONV')then
  call atoi(ch(6:),iexp)
  scftol=10**(-float(iexp))

 
 elseif(ch=='DEORTHO')then
  DEORTHO=.true.
 elseif(ch=='GRADIENT')then
  GRADIENT=.true.
 elseif(ch=='TRIAL')then
  TRIAL=.true.
 elseif(ch=='XYZ')then
  XYZ=.true.
  elseif(ch=='OPTIMIZE')then
  OPTIMIZE=.true.

elseif(ch(1:6)=='OPTTOL')then
  call atoi(ch(8:),iexp)
  OPTTOL=10**(-float(iexp))
elseif(ch(1:6)=='CHARGE')then
  if(ch(8:8)=='-')then
     call atoi(ch(9:),iexp)
     SYS_CHARGE=-float(iexp)
  else
     call atoi(ch(8:),iexp)
     SYS_CHARGE=float(iexp)
     end if

elseif(ch(1:4)=='MULT')then
  call atoi(ch(6:),iexp)
  mult=int(iexp)

elseif(ch=='PARAMETERS')then
  parameters=.true.
elseif(ch=='GAUSSIAN')then
  core='GAUSSIAN'

elseif(ch=='3CENTER')then
  CENTERS=.true.
elseif(ch=='NEWTON')then
  newton=.true.
elseif(ch=='UHF')then
  restricted=.false.
elseif(ch=='HFDFT')then
  DFTHF=.true.
elseif(ch=='LAPACK')then
  LAPACK=.TRUE.
elseif(ch=='SPARKLES')then
  SPARKLES=.TRUE.
elseif(ch=='DENSITYIN')then
  densityin=.TRUE.
elseif(ch=='DENSITYOUT')then
  densityout=.TRUE.
elseif(ch=='FREQUENCY')then
  FREQUENCY=.TRUE.
elseif(ch=='TS')then
  TS=.TRUE.
elseif(ch=='DIIS')then
  DIIS=.TRUE.
elseif(ch=='MODES')then
 MODES=.TRUE.
elseif(ch=='CALCALL')then
 CALCALL=.TRUE.
elseif(ch(1:6)=='UPDATE')then
  call atoi(ch(8:),iexp)
  UPDATE=IEXP
elseif(ch=='IRC')then
 IRC=.TRUE.
elseif(ch=='EULER')then
 EULER=.TRUE.
elseif(ch=='RK4')then
 RK4=.TRUE.
elseif(ch=='FORWARD')then
 FORWARD=.TRUE.
elseif(ch=='SAM')then
 SAM=.TRUE.
elseif(ch(1:4)=='DEL1')then
  call atoi(ch(6:),iexp)
  DEL1=10**(-float(iexp))
elseif(ch(1:4)=='STE1')then
  call atoi(ch(6:),iexp)
  STE1=10**(-float(iexp))
elseif(ch(1:5)=='BASIS')then
 basis_set=ch(7:) 
elseif(ch(1:7)=='ERIBUFF')then
length=len_trim(ch)
call char_to_float(ch,valuedp,length)

ERIBUFF=int(valuedp)

elseif(ch(1:6)=='CUTOFF')then
length=len_trim(ch)
call char_to_float(ch,valuedp,length)
cutoff=valuedp
cutoff_sq=cutoff*cutoff
elseif(ch(1:5)=='DEBUG')then
debug=.true.
elseif(ch(1:5)=='CGSCF')then
CGSCF=.true.
elseif(ch(1:6)=='SDTOCG')then
length=len_trim(ch)
call char_to_float(ch,valuedp,length)
sdtocg=valuedp
elseif(ch(1:3)=='PBC')then
pbc=.true.
elseif(ch(1:6)=='LINDEP')then
length=len_trim(ch)
call char_to_float(ch,valuedp,length)
lindep_tol=valuedp


 end if




end do
end subroutine default



