       subroutine bfgsq(natoms,nopt,force,q,internal,icount2,numint
     $,intadd)
      integer          nmax, mmax,iwrite,natoms,nopt
      parameter        (nmax=1024, mmax=17)
c        nmax is the dimension of the largest problem to be solved.
c        mmax is the maximum number of limited memory corrections.
 
c     Declare the variables needed by the code.
c       A description of all these variables is given at the end of 
c       the driver.
 
      character*60     task, csave
      logical          lsave(4),stat1,stat2
      integer          n, m, iprint,ii,icount2,
     +                 nbd(nopt), iwa(3*nopt), isave(44)
      double precision f, factr, pgtol, force,delta,fnew,fold,qold,
     +               x(nopt), l(nopt), u(nopt), g(nopt), dsave(29), 
     +                 wa(2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax)
     $,rnorm

c     Declare a few additional variables for this sample problem.

      integer          i,k,NCYCLE,kk
      double precision x0(natoms),y0(natoms),z0(natoms),q(natoms,3)
     $,gint(numint)
      integer internal(2,nopt),numint,intadd(numint,2)
c     We wish to have output at every iteration.
      iwrite=2
      iprint = 1

c     We specify the tolerances in the stopping criteria.

c      factr=1.0d+7
c      pgtol=1.0d-5
      pgtol=0.0D0
c     We specify the dimension n of the sample problem and the number
c        m of limited memory corrections stored.  (n and m should not
c        exceed the limits nmax and mmax respectively.)
c
c     the dimension of our problem is the 3N cartesian (3N-6 if I were smart)
c     coordinates
c
c      n=3*natoms
      n=nopt
      m=20
c      delta=.000572957d0
c      delta=12.0d0
        delta=.00001d0
C
C     DEFINE STARTING POINT AND THE VARIABLE BOUNDS
C
c     input coordinates are on file xyz as follows:
c     atomic# x y z -> for each atom in ANGSTROMS!
c
      do 100 i=1,nopt
         x(i)=q(  internal(1,i), internal(2,i)  )
 100  continue
      close(10)

c     set bounds on variables (none for optimization) and convert to bohr
       

      do 10 i=1,n
         nbd(i)=0
         if(i.le.icount2)then
         nbd(i)=1
         l(i)=0.5d0
         end if
  10  continue


c 
      task = 'START'
      NCYCLE=0
      print*,''
      print*,'///////////////////////////////////////////////////'
      print*,'             Geometry optimization in internals'
      print*,'    Optimization considered complete when the      '
      print*,'    gradient norm is less than ',force
      print*,'///////////////////////////////////////////////////'

      
      print*,'There are ',nopt,' optimizable degrees of freedom'

c        ------- the beginning of the loop ----------
 
 111  continue
      
c     This is the call to the L-BFGS-B code.
 
      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
     +            csave,lsave,isave,dsave)
 
      if (task(1:2) .eq. 'FG') then
c        the minimization routine has returned to request the
c        function f and gradient g values at the current x.
c
c
c
         NCYCLE=NCYCLE+1
         if(NCYCLE.GT.100)STOP
          print*,''
         PRINT*,'-----------------------------------------------'
         PRINT*,'GEOMETRY OPTIMIZATION CYCLE NUMBER: ',NCYCLE
         do 510 i=1,nopt
        q(  internal(1,i), internal(2,i)  ) = x(i)
 510     continue
          stat1=.false.
          stat2=.true.
        call scfq(q,natoms,f,stat1,stat2,gint,numint)
        f=f*627.51d0/27.21d0
      do i=1,nopt
           do j=1,numint
           if(intadd(j,1).eq.internal(1,i).and.
     $     intadd(j,2).eq.internal(2,i))goto 55
           end do
 55        g(i)=gint(j)
      end do
         
          
 

c             stop






c
c     calculate ||F|| of current point to monitor optimization
c
         rnorm=0.0D0
         do 80 i=1,n
         rnorm=rnorm+g(i)**2
 80      continue
         rnorm=dsqrt(rnorm)
         print*,''
         print*,'****************************************************'
         print*,'*       AT THE END OF CYCLE NUMBER ',NCYCLE,':      '
c         PRINT*,'* GRADIENT NORM OF OPTIMIZED COORDINATES (Kcal/(mol
c!     $*Ang)= ',RNORM
         WRITE(*,150)RNORM
 150  FORMAT('GRADIENT NORM OF OPTIMIZED COORDINATES (Kcal/(mol*Ang) ='
     $,F10.5)
         PRINT*,'*            ENERGY = ',F,' Kcal/mol                '
         print*,'*****************************************************'
         


          if(rnorm.lt.force)then
c          if(rnorm.lt.5.0d0)then
         write(*,*)'OPTIMIZATION CONVERGED AFTER',ncycle,' CYCLES'
        write(*,*)''
c put the optimal coordinates in the zmat
        do 500 i=1,nopt
        q(  internal(1,i), internal(2,i)  ) = x(i)
500     continue
      
        
         return
         end if

         
         
         
         





c          go back to the minimization routine.
         goto 111
      endif
c
      if (task(1:5) .eq. 'NEW_X') then
         goto 111
       end if

c        the minimization routine has returned with a new iterate,
c         and we have opted to continue the iteration.

c           ---------- the end of the loop -------------
 
c     If task is neither FG nor NEW_X we terminate execution.
c
c
c    if we get here then optimization failed
         write(*,*)'*********OPTIMIZATION FAILED********'

      return
 
      end
      
      
         
      




c======================= The end of driver1 ============================

c     --------------------------------------------------------------
c             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     n is an INTEGER variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an INTEGER variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended, and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended. 
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit, it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound, l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound, u(i) need not be defined.
c
c     nbd is an INTEGER array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds, 
c              3 if x(i) has only an upper bound.
c
c     f is a DOUBLE PRECISION variable.  If the routine setulb returns
c       with task(1:2)= 'FG', then f must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb
c       returns with taskb(1:2)= 'FG', then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy; 
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length 
c       (2mmax + 4)nmax + 12mmax^2 + 12mmax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an INTEGER  array of length 3nmax used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry, it must be set to 'START'.
c       On a return with task(1:2)='FG', the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X', an iteration of the
c         algorithm has concluded, and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration. 
c       When
c         task(1:4)='CONV', the termination test in L-BFGS-B has been 
c           satisfied;
c         task(1:4)='ABNO', the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR', the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run 
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a LOGICAL working array of dimension 4.
c       On exit with task = 'NEW_X', the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an INTEGER working array of dimension 44.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in isave
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X', it contains information that
c         the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in dsave
c
c     --------------------------------------------------------------
c           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     << An example of subroutine 'timer' for AIX Version 3.2 >>
c
c     subroutine timer(ttime)
c     double precision ttime
c     integer itemp, integer mclock
c     itemp = mclock()
c     ttime = dble(itemp)*1.0d-2
c     return
c     end
c-----------------------------------------------------------------------
