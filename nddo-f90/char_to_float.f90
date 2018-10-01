subroutine char_to_float(word,value,length)
!/////////////////////////////////////////////////////
!convert a character string to a floating point number
!/////////////////////////////////////////////////////
      implicit double precision (a-h,o-z)
      character(200),intent(in)::word
      integer,intent(in)::length
      double precision,intent(inout)::value

      nstart=scan(word,'=')
      nstart=nstart+1
      sign=1.0d0
      negative=scan(word,'-')
      if(negative/=0)then
      sign=-1.0d0
      nstart=nstart+1
      end if
      ndecimal=scan(word,'.')

      if(ndecimal.eq.0)then
       write(*,*)'Floating point input required'
       write(*,*)'Error:',word
       stop
      end if

! compute floating point total of keyword
      value=0.0d0
      nbefore=ndecimal-nstart
      iend=ndecimal-1
      do i=iend,nstart,-1
      j=ndecimal-i-1
      power=float(j)
      call atof(word(i:i),s)
      value=value+s*10**power
      end do

      istart=ndecimal+1
      do i=istart,length
      j=i-ndecimal
      power=float(j)*(-1.0d0)
      call atof(word(i:i),s)
      value=value+s*10**power
      end do

      value=value*sign

return
end subroutine char_to_float

subroutine atof(input,s)
character,intent(in)::input
double precision,intent(inout)::s
if(input=='0')then
s=0.
elseif(input=='1')then
s=1.
elseif(input=='2')then
s=2.
elseif(input=='3')then
s=3.
elseif(input=='4')then
s=4.
elseif(input=='5')then
s=5.
elseif(input=='6')then
s=6.
elseif(input=='7')then
s=7.
elseif(input=='8')then
s=8.
elseif(input=='9')then
s=9.
end if
end subroutine atof
