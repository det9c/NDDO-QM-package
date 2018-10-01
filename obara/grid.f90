subroutine gen_grid
use gaussian_basis
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
integer rnum2
pi=dacos(-1.0d0)
!     generation of interpolation table using explicit
!     formulas given by frank harris
!     m is given in harris
      w=15.0D0
      rnum=(2.0D0*w)+1.0D0
      rnum2=(2.0D0*w)-1.0D0
      total=1.0D0
!****************
!  this calculates (2m-1)!! [num2 above]
 do 155 l=rnum2,1,-2
        total=total*l
 155  continue
       total=total*1.0D0
!********************
!  this calculates (2m+2k+1)!!
         do 1000 i=1,100
            dumber=0.0
         dumber=(rnum+2.0D0*(i-1))
         fac=1.0D0
         do 11 j=dumber,1,-2
         fac=fac*j
 11      continue
         denom(i)=fac
 1000 continue

!*****************
         sum=0.0D0
         do 16 n=1,4000
            t=0.0
         t=(n-1.0D0)*.01D0
         sum=0.0
         do 12 k=1,100
            sum=sum+(((2.0D0*t)**(k-1))/denom(k))
 12         continue
            table(n)=sum*total*exp(-t)
!            print*,'for t =',t,'f is',table(n)
 16         continue
            return
            end
