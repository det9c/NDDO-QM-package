subroutine make_hess_sam
use constants
use tables
use scratch_array
use control
implicit double precision (a-h,o-z)
interface
subroutine scfopt(x0,y0,z0,gout,natoms,ff)
double precision,intent(inout) ::ff
double precision, dimension(natoms),intent(inout) ::x0,y0,z0
double precision, dimension(3*natoms),intent(inout) ::gout
end subroutine scfopt

subroutine mkvector(matrix,vector,idim1)
use indices
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::matrix
double precision,dimension(:),intent(inout)::vector
integer,intent(in)::idim1
end subroutine mkvector




end interface

double precision,allocatable,dimension(:,:)::hessian,gplus,gminus,hessold
double precision,allocatable,dimension(:)::vec
!double precision,allocatable,dimension(:)::x0,y0,z0



delta=1D-3
i3n=3*numat
iupper=i3n*(i3n+1)/2
allocate(vec(iupper))
allocate(hessian(i3n,i3n))
allocate(gplus(i3n,1))
allocate(gminus(i3n,1))

optimize=.false. ! set this to fool program into not printing xyz all the 
hessian=zero
icol=0

SAVE_TREE=.true.
do i=1,numat
!print*,'Working on components for atom',i
icol=icol+1
! compute d2E/[dx(i) dq(j)] where q=x,y,z
xold=x(i)
x(i)=x(i)+delta
call scfopt(x,y,z,gplus,numat,ff)
x(i)=x(i)-2.0d0*delta
call scfopt(x,y,z,gminus,numat,ff)
gplus=(gplus-gminus)/two/delta
x(i)=xold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)

icol=icol+1
! compute d2E/[dy(i) dq(j)] where q=x,y,z
yold=y(i)
y(i)=y(i)+delta
call scfopt(x,y,z,gplus,numat,ff)
y(i)=y(i)-2.0d0*delta
call scfopt(x,y,z,gminus,numat,ff)
gplus=(gplus-gminus)/two/delta
y(i)=yold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)



icol=icol+1
! compute d2E/[dz(i) dq(j)] where q=x,y,z
zold=z(i)
z(i)=z(i)+delta
call scfopt(x,y,z,gplus,numat,ff)
z(i)=z(i)-2.0d0*delta
call scfopt(x,y,z,gminus,numat,ff)
gplus=(gplus-gminus)/two/delta
z(i)=zold
hessian(1:i3n,icol:icol)=gplus(1:i3n,1:1)

end do


hessian=hessian*.529177d0*.529177d0/627.51d0 ! convert to atomic units

! symmetrize hessian
hessian=hessian+transpose(hessian)
hessian=hessian/two

!call matprt(hessian,i3n,i3n,i3n,i3n)
call mkvector(hessian,vec,i3n)

call dcopy(iupper,vec,1,rpd2xyz(1,ircpt),1)

deallocate(hessian)
deallocate(gplus)
deallocate(gminus)
deallocate(vec)
end subroutine make_hess_sam


