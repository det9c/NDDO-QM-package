subroutine outerp(n,alpha,v1,v2,mat)
implicit double precision (a-h,o-z)
double precision,dimension(n,1),intent(in)::v1,v2
double precision,dimension(n,n),intent(inout)::mat
double precision,intent(in)::alpha
integer,intent(in)::n

do i=1,n
do j=1,n
mat(j,i)=mat(j,i)+alpha*v1(j,1)*v2(i,1)
end do
end do

end subroutine outerp
