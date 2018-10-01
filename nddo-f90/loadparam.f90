subroutine loadparam
use tables
character(6)::a
double precision::value
integer::zz,jparam
write(*,*)
write(*,*)
write(*,*)'+++++++++++++++++++++++++++++++++++++++++++++++++'
write(*,*)'Parameters to be read from external source...'
write(*,*)
write(*,*)'User has requested that the following parameters overwrite the defaults:'
itype=size(ntype,1)
open(unit=10,file='PARAMETERS')
! read parameter file and load values until there are no more lines to be read
do
read(10,*,iostat=io)a,zz,value
if(io<0)exit
write(*,20)a,zz,value

! locate row on which parameter is stored
jparam=2
do i=1,itype
if(ntype(i)==zz)then
jparam=1
exit
end if
end do
if(jparam==1)then
!now assign parameter to the proper row
if(a=='USS')then
uss(i)=value
else if(a=='UPP')then
upp(i)=value
else if(a=='UDD')then
udd(i)=value
else if(a=='ZS')then
zs(i)=value
else if(a=='ZP')then
zp(i)=value
else if(a=='ZD')then
zetad(i)=value
else if(a=='BETAS')then
betas(i)=value
else if(a=='BETAP')then
betap(i)=value
else if(a=='BETAD')then
betad(i)=value
else if(a=='GSS')then
gss(i)=value
else if(a=='GPP')then
gpp(i)=value
else if(a=='GSP')then
gsp(i)=value
else if(a=='GP2')then
gp2(i)=value
else if(a=='HSP')then
hsp(i)=value
else if(a=='ZSN')then
zsone(i)=value
else if(a=='ZPN')then
zpone(i)=value
else if(a=='ZDN')then
zdone(i)=value
else if(a=='ALPHA')then
alpha(i)=value
! PUT IN MOPAC FNXX PARAMETERS
else if(a=='A1' .or. a=='FN11')then
g1(1,i)=value
else if(a=='B1' .or. a=='FN21')then
g2(1,i)=value
else if(a=='C1' .or. a=='FN31')then
g3(1,i)=value
else if(a=='A2' .or. a=='FN12')then
g1(2,i)=value
else if(a=='B2'.or. a=='FN22')then
g2(2,i)=value
else if(a=='C2'.or. a=='FN32')then
g3(2,i)=value
else if(a=='A3'.or. a=='FN13')then
g1(3,i)=value
else if(a=='B3'.or. a=='FN23')then
g2(3,i)=value
else if(a=='C3'.or. a=='FN33')then
g3(3,i)=value
else if(a=='ATHREE')then
alpha3(i)=value
else if(a=='A1C')then
ac(1,i)=value
else if(a=='B1C')then
bc(1,i)=value
else if(a=='A2C')then
ac(2,i)=value
else if(a=='B2C')then
bc(2,i)=value
else if(a=='A4'.or. a=='FN14')then
g1(4,i)=value
else if(a=='B4'.or. a=='FN24')then
g2(4,i)=value
else if(a=='C4'.or. a=='FN34')then
g3(4,i)=value


end if
end if
end do





20 format(A10,I2,F10.5)
end subroutine loadparam

