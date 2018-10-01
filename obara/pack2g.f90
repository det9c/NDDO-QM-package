subroutine pack2g(indi,indj,ij)
use gaussian_basis      

!     this subroutine packs indices i and j into a single integer
!     note that i is STRICTLY <= j 


integer:: i,j,k,l,nlen
integer,intent(in)::indi,indj
integer,intent(inout)::ij
      i=indi
      j=indj
         if(indi.ne.indj)then
         i=min0(indi,indj)
         j=max0(indi,indj)
         end if
      ij=i+   ioffset2(j)            !!!j*(j-1)/2       !ioffset(j)


      return
      end
      
