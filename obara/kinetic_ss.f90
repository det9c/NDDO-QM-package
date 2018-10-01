subroutine kinetic_ss(dist,r,t,ss,dkssout)
use gaussian_basis
implicit double precision (a-h,o-z)
pi=dacos(-1.0d0)
rplust=r+t
s=r*t/rplust
zeta_global=s
dkssout=s*(3.0d0-2*s*dist)*ss

return
end
