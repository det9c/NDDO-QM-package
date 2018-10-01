subroutine overlap_ss(dist,r,t,ssout)
implicit double precision (a-h,o-z)
pi=dacos(-1.0d0)
rplust=r+t
ssout=((pi/rplust)**(1.5D0))*dexp(-r*t*dist/rplust)
return
end
