subroutine fzero(uz,series,low)
use gaussian_basis
implicit double precision(a-h,o-z)
double precision,intent(in)::uz
double precision,intent(inout)::series
integer,intent(in)::low

dimension func(20),sfac(7)
pi=dacos(-1.0d0)
index=uz/.01D0
ir=int(index)
pk=table(ir+1)
tnew=float(ir)*.01D0
delta=tnew-uz
func(16)=pk
r1=2.0D0*tnew
r2=exp(-1.0D0*tnew)
do 13 j=15,1,-1
    func(j)=0.0D0
    func(j)=((r1*func(j+1))+r2)/(2.0D0*(j-1)+1.0D0)
 13            continue
!     explicitly putting in kfactorial for taylor series
               sfac(1)=1.0D0
               sfac(2)=1.0D0
               sfac(3)=2.0D0
               sfac(4)=6.0D0
               sfac(5)=24.0D0
               sfac(6)=120.0D0
               sfac(7)=720.0D0
               ip=-1
               tseries=0.0D0
               series=0.0D0
               lot=low+1
               ihigh=low+7
!               do 20 k=lot,ihigh
!                  p=p+1.0D0
!              tseries=tseries+(func(k)*(delta)**(p))/sfac(p+1)
! 20               continue

               do 20 k=lot,ihigh
                  ip=ip+1
                  rfac=1.
                  do mm=1,ip
                  rfac=rfac*delta   
                  end do
              tseries=tseries+(func(k)*rfac)/sfac(ip+1)
 20               continue



             series=tseries
                  return
                  end

