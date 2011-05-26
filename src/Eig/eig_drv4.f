      subroutine eig_drv4 (n,nev,ncv,sigma,d,v,work,ierr)
c
c $Id$
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: eig_drv4
c
c\Description:
c     Driver to solve a generalized eigenvalue problem using the
c     reverse communication interface in ARPACK, shift and invert mode.
c     This subroutine is based on the test driver DSDRV4 from the
c     symmetric examples of the ARPACK distribution.
c
c     ... Suppose we want to solve A*x = lambda*M*x in shift-invert mode,
c         where A and M are stiffness and mass matrix, respectively.
c
c     ... OP = (inv[A - sigma*M])*M  and  B = M.
c
c     ... Use mode 3 of DSAUPD.
c
c\Usage:
c  call eig_drv4 ( N, NEV, NCV, SIGMA, D, V, WORK, IERR )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c
c  SIGMA   Double precision.  (INPUT)
c          The shift value.
c
c  D       Double precision  NCV by 2 array.  (OUTPUT)
c          The first column of D contain the eigenvalues.
c
c  V       Double precision  N by NCV array.  (OUTPUT)
c          The NCV columns of V contain the Lanczos basis vectors.
c
c  WORK    Double precision  work array.  (OUTPUT/WORKSPACE)
c
c
c  IERR    Integer.  (OUTPUT)
c          Error flag.
c
c\EndDoc
c-----------------------------------------------------------------------
c\BeginLib
c
c\Routines called:
c     dsaupd  ARPACK reverse communication interface routine.
c     dseupd  ARPACK routine that returns Ritz values and vectors.
c     daxpy   Level 1 BLAS that computes y <- alpha*x+y.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     eig_sol Linear solver routine that computes inv[A - sigma*M]*x
c     eig_av  Matrix vector multiplication routine that computes A*x.
c     eig_mv  Matrix vector multiplication routine that computes M*x.
c
c\Author
c     Knut Morten Okstad
c     Dept. of Applied Mathematics
c     SINTEF ICT
c     Trondheim, Norway
c
c\EndLib
c-----------------------------------------------------------------------
c
      implicit none
c
c     %-----------%
c     | Arguments |
c     %-----------%
C
      integer          n, nev, ncv, ierr
      Double precision sigma, d(ncv,2), v(n,ncv), work(*)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer          iparam(11), ipntr(11)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        bmat*1, which*2
      integer          ipresid, ipselec, ipworkd, ipworkl, lworkl,
     &                 ido, j, nconv, maxitr, ishfts, mode
      Double precision tol
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision dnrm2
      external         dnrm2, daxpy, eig_sol, eig_av, eig_mv
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic        abs
c
c     %-----------------------%
c     | Executable statements |
c     %-----------------------%
c
      print *, ' '
      print *, ' Entering EIG_DRV4',n,nev,ncv,sigma
c
c     %----------------------------------------------------%
c     | The number N is the dimension of the matrix.  A    |
c     | generalized eigenvalue problem is solved (BMAT =   |
c     | 'G'.) NEV is the number of eigenvalues (closest to |
c     | the shift SIGMA) to be approximated.  Since the    |
c     | shift-invert mode is used, WHICH is set to 'LM'.   |
c     | The user can modify NEV, NCV, SIGMA to solve       |
c     | problems of different sizes, and to get different  |
c     | parts of the spectrum.                             |
c     %----------------------------------------------------%
c
      bmat = 'G'
      which = 'LM'
c
c     %--------------------------------------------------%
c     | The work array WORKL is used in DSAUPD as        |
c     | workspace.  Its dimension LWORKL is set as       |
c     | illustrated below.  The parameter TOL determines |
c     | the stopping criterion.  If TOL<=0, machine      |
c     | precision is used.  The variable IDO is used for |
c     | reverse communication and is initially set to 0. |
c     | Setting INFO=0 indicates that a random vector is |
c     | generated in DSAUPD to start the Arnoldi         |
c     | iteration.                                       |
c     %--------------------------------------------------%
c
      lworkl = ncv*(ncv+8)
      tol = 0.0D0
      ido = 0
      ierr = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 3 specified in the      |
c     | documentation of DSAUPD is used (IPARAM(7) = 3).  |
c     | All these options may be changed by the user.     |
c     | For details, see the documentation in DSAUPD.     |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 3
c
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode
c
      ipworkd = 1
      ipworkl = ipworkd + 3*n
      ipresid = ipworkl + lworkl
      ipselec = ipresid + n
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     %-------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DSAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
         call dsaupd (ido, bmat, n, which, nev, tol, work(ipresid),
     &                ncv, v, n, iparam, ipntr, work(ipworkd),
     &                work(ipworkl), lworkl, ierr)

         if (ido .eq. -1) then
c
c           %--------------------------------------------%
c           | Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x  |
c           | to force the starting vector into the      |
c           | range of OP.  The user should supply       |
c           | his/her own matrix vector multiplication   |
c           | routine and a linear system solver here.   |
c           | The matrix vector multiplication routine   |
c           | takes workd(ipntr(1)) as the input vector. |
c           | The final result is returned to            |
c           | workd(ipntr(2)).                           |
c           %--------------------------------------------%
c
            call eig_mv (work(ipntr(1)), work(ipntr(2)))
            call eig_sol (work(ipntr(2)), work(ipntr(2)), ierr)
            if (ierr .ne. 0) then
               print *, ' '
               print *, ' Error with eig_sol in EIG_DRV4'
               print *, ' '
               return
            end if
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         else if (ido .eq. 1) then
c
c           %-----------------------------------------%
c           | Perform y <-- OP*x = inv[A-sigma*M]*M*x |
c           | M*x has been saved in workd(ipntr(3)).  |
c           | the user only needs the linear system   |
c           | solver here that takes workd(ipntr(3))  |
c           | as input, and returns the result to     |
c           | workd(ipntr(2)).                        |
c           %-----------------------------------------%
c
            call eig_sol (work(ipntr(3)), work(ipntr(2)), ierr)
            if (ierr .ne. 0) then
               print *, ' '
               print *, ' Error with eig_sol in EIG_DRV4'
               print *, ' '
               return
            end if
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         else if (ido .eq. 2) then
c
c           %-----------------------------------------%
c           |          Perform  y <--- M*x            |
c           | Need the matrix vector multiplication   |
c           | routine here that takes workd(ipntr(1)) |
c           | as the input and returns the result to  |
c           | workd(ipntr(2)).                        |
c           %-----------------------------------------%
c
            call eig_mv (work(ipntr(1)), work(ipntr(2)))
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if
c
c     %-----------------------------------------%
c     | Either we have convergence, or there is |
c     | an error.                               |
c     %-----------------------------------------%
c
      if (ierr .lt. 0) then
c
c        %--------------------------%
c        | Error message, check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with dsaupd, info = ',ierr
         print *, ' Check the documentation of dsaupd'
         print *, ' '
         return
c
      else if (ierr .eq. 1) then
         print *, ' '
         print *, ' Maximum number of iterations reached.'
         print *, ' '
      else if (ierr .eq. 3) then
         print *, ' '
         print *, ' No shifts could be applied during implicit',
     &        ' Arnoldi update, try increasing NCV.'
         print *, ' '
      end if
c
c     %-------------------------------------------%
c     | No fatal errors occurred.                 |
c     | Post-Process using DSEUPD.                |
c     |                                           |
c     | Computed eigenvalues may be extracted.    |
c     |                                           |
c     | Eigenvectors may also be computed now if  |
c     | desired.  (indicated by rvec = .true.)    |
c     %-------------------------------------------%
c
      call dseupd (.true., 'All', work(ipselec), d, v, n, sigma,
     &             bmat, n, which, nev, tol, work(ipresid),
     &             ncv, v, n, iparam, ipntr,
     &             work(ipworkd), work(ipworkl), lworkl, ierr)
c
c     %----------------------------------------------%
c     | Eigenvalues are returned in the first column |
c     | of the two dimensional array D and the       |
c     | corresponding eigenvectors are returned in   |
c     | the first NEV columns of the two dimensional |
c     | array V if requested.  Otherwise, an         |
c     | orthogonal basis for the invariant subspace  |
c     | corresponding to the eigenvalues in D is     |
c     | returned in V.                               |
c     %----------------------------------------------%
c
      if (ierr .ne. 0) then
c
c        %------------------------------------%
c        | Error condition:                   |
c        | Check the documentation of DSEUPD. |
c        %------------------------------------%
c
         print *, ' '
         print *, ' Error with dseupd, info = ', ierr
         print *, ' Check the documentation of dseupd'
         print *, ' '
         return
c
      end if
c
      nconv = iparam(5)
      do 30 j = 1, nconv
c
c        %---------------------------%
c        | Compute the residual norm |
c        |                           |
c        |  || A*x - lambda*M*x ||   |
c        |                           |
c        | for the NCONV accurately  |
c        | computed eigenvalues and  |
c        | eigenvectors.  (iparam(5) |
c        | indicates how many are    |
c        | accurate to the requested |
c        | tolerance)                |
c        %---------------------------%
c
         call eig_av (v(1,j), work)
         call eig_mv (v(1,j), work(n+1))
         call daxpy (n, -d(j,1), work(n+1), 1, work, 1)
         d(j,2) = dnrm2(n, work, 1)
         d(j,2) = d(j,2) / abs(d(j,1))
c
 30   continue
c
      call dmout (6, nconv, 2, d, ncv, -6,
     &            'Ritz values and relative residuals')
c
c
c     %------------------------------------------%
c     | Print additional convergence information |
c     %------------------------------------------%
c
      print *, ' '
      print *, ' EIG_DRV4'
      print *, ' ========'
      print *, ' '
      print *, ' Size of the matrix is ', n
      print *, ' The number of Ritz values requested is ', nev
      print *, ' The number of Arnoldi vectors generated',
     &     ' (NCV) is ', ncv
      print *, ' What portion of the spectrum: ', which
      print *, ' The number of converged Ritz values is ', nconv
      print *, ' The number of Implicit Arnoldi update',
     &     ' iterations taken is ', iparam(3)
      print *, ' The number of OP*x is ', iparam(9)
      print *, ' The convergence criterion is ', tol
      print *, ' '
c
c     %----------------------------%
c     | Done with program eig_drv4 |
c     %----------------------------%
c
      end
