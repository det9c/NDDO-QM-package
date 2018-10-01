subroutine pulay_ab(oldfock,error,fock,iter,ldim)
use gaussian_basis
implicit double precision (a-h,o-z)
double precision,dimension(:,:),intent(in)::oldfock,error
double precision,dimension(:),intent(inout)::fock
integer,intent(in)::iter,ldim
double precision,dimension(ldim,ldim)::amatrix,ainv,atemp,unit
double precision,dimension(ldim)::bvec,cvec  ! solving Ab=c
integer,dimension(ldim)::ind
double precision,dimension(iveclength)::vec1
zero=0.

amatrix=1.0d0
bvec=zero
cvec=zero
nlen=num_gauss**2
do j=1,iter
cvec(j)=zero
do i=1,j
top=ddot(nlen,error(1,j),1,error(1,i),1)
amatrix(j,i)=top
amatrix(i,j)=top
end do
end do
cvec(ldim)=1.0d0
amatrix(ldim,ldim)=zero



call gaussj(amatrix,ldim,ldim,cvec,1,1)



fock=zero
do i=1,iter
call dcopy(iveclength,oldfock(1,i),1,vec1,1)
call dscal(iveclength,cvec(i),vec1,1)
call daxpy(iveclength,1.0d0,vec1,1,fock,1)
end do









return
end subroutine pulay_ab


