      subroutine DNQSOL(DNQFJ, N, X, FVEC, XTOL, IOPT, W, IDIMW)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 2001-05-25 DNQSOL Krogh Minor change for making .f90 version.
c>> 2000-12-01 DNQSOL Krogh  Removed unused parameter P001.
c>> 1996-05-16 DNQSOL Krogh  Changes to use .C. and C%%.
c>> 1996-03-30 DNQSOL Krogh  Added external stmts. SIN => VSIN, etc.
c>> 1994-11-02 DNQSOL Krogh  Changes to use M77CON
c>> 1992-04-27 DNQSOL CLL Deleted unreferenced stmt label.
c>> 1992-04-07 CAO Extra comma in Print removed (error from VAX compile)
c>> 1992-01-15 CLL
c>> 1991-12-18 CLL & FTK  Adding treatment of slow convergence to 0.
c>> 1991-12-05 CLL & FTK  Adding Option vector interface.
c>> 1990-04-20 CLL@JPL Adapting code from Minpack for MATH77
c
c  Solves a system of N nonlinear equations in N unknowns.
c  DNQSOL is the the user-interface subroutine.  It calls DNQSL1 which
c  contains the top-level logic of the solution algorithm.
c  DNQSOL & DNQSL1 also need:
c         Other subroutines that are in this file:
c            DNQFDJ, DNQDOG, DNQQFM, DNQQRF, DNQUPD.
c         Other subprograms from the MATH77 library: DNRM2, DERV1,
c            [D/R]1MACH (Fortan 77 only), IERV1, & IERM1.
C         A user-provided subroutine: DNQFJ.
c
c  Most of these subprograms are derived from MINPACK-1.
c  MINPACK-1, 1980, was developed by Jorge J. More',
c  Burton S. Garbow, and Kenneth E. Hillstrom, Argonne Nat'l Lab.
c  The MINPACK-1 code was obtained as FILE05 from MINPACK/EX from
c  Netlib, downloaded to JPL on Tue Feb  6 12:17:45 EST 1990.
c
c     Old Name         New Name
c     --------         --------
c     HYBRJ1, HYBRD1   DNQSOL (Completely redesigned.)
c     HYBRJ, HYBRD     DNQSL1 (Algorithm and code changes.)
c     DOGLEG           DNQDOG
c     ENORM            DNRM2 in BLAS and MATH77
c     FDJAC1           DNQFDJ
c     QFORM            DNQQFM
c     QRFAC            DNQQRF
c     R1MPYQ           DNQAQ
c     R1UPDT           DNQUPD
c     [D/S]PMPAR       [D/R]1MACH in file amach.f (Fortran 77 only)
c     FCN              DNQFJ
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c                  Arguments for DNQSOL
c
c           call DNQSOL(DNQFJ, N, X, FVEC, XTOL, IOPT, W, IDIMW)
c
c  DNQFJ    Name of user-supplied subroutine.
c
c  N        [in]     Problem size
c  X(N)     [inout]  Initial and final x-vector.
c  FVEC(N)  [out]    Final F values.
c  XTOL     [in]     Rel. Conv. tolerance on weighted X
c  IOPT()   [inout]  First 3 elements contain output values.
c           IOPT(1) = INFO.  Output status.
c           IOPT(2) = NFEV.  No. of F evals used.
c           IOPT(3) = NJEV.  No. of evals of Jacobian.
c
c                    Ramaining elements in IOPT() select options.
c
c     Option   No. of       Affected variables       Affected variables
c     Number   arguments    in DNQSOL.               in DNQSL1.
c        1         0        HAVEJ                    HAVEJ
c        2         1        DMODE, HAVED, W(4:3+N)   HAVED, DIAG(1:N)
c        3         1        NPRINT                   NPRINT
c        4         1        MAXFEV                   MAXFEV
c        5         2        ML, MU                   ML, MU
c        6         0        W(1)                     EPSFCN
c        7         0        W(2)                     FACTOR
c        8         0        TRACE                    TRACE
c
c Functionality of options, listed by option numbers in square brackets.
c       [1]  If set, user is not computing a Jacobian.
c            This subr sets HAVEJ = .false.
c       [2]  Arg: DMODE = 1, 2, or 3.
c            1. This subr sets DIAG() to all ones and HAVED = .true.
c            2. User has set DIAG().  This subr sets HAVED = .true.
c            3. This subr sets HAVED = .false. so DNQSL1 will set
c               DIAG() dynamically.
c       [3]  Arg: NPRINT   Print control.
c       [4]  Arg: MAXFEV   Limit on no. of F evals.
c       [5]  Args: ML & MU  Band structure.
c       [6]  If set means EPSFCN has been set in W(1).
c       [7]  If set means FACTOR has been set in W(2).
c       [8]  If set, this subr sets TRACE = .true., else this subr
c            sets TRACE = .false.  When TRACE is .true., DNQSL1 prints
c            detailed intermediate results.
c
c  W()  [inout]  W(1) and W(2) may be used to pass EPSFCN and FACTOR
c       to the subroutine.  W(3) contains TOLTST on return.
c       W( 4 : 3+(15*N + 3*N**2)/2 ) is used as work space.
c
c       EPSFCN     W(1)     Error in F evals.  Used in computing
c                              approx derivs.
c       FACTOR     W(2)     Algorithm parameter.
c       TOLTST     W(3)     Output.  Final value of quantity compared
c                              with XTOL for convergence test.
c       DIAG(N)    W(4:N+3) Scaling values.  May be input or
c                              computed.  See option 2.
c       WA1(N)     W()      Work space of length N.
c       WA2(N)     W()      Work space of length N.
c       WA3(N)     W()      Work space of length N.
c       WA4(N)     W()      Work space of length N.
c       GNSTEP(N)  W()      Work space of length N.
c       QTF(N)     W()      Wrk space.  At end has (Q**t)*F.
c       FJAC(N,N)  W() Work space for Jacobian.  At end has Q of
c                      QR factorization.
c       R( (N + N**2)/2 )  W()  Wrk spc.  At end has Packed R of
c                      QR factorization.
c  IDIMW  [in]  Dimension of W().  Require IDIMW .ge. 3+(15*N+3*N**2)/2
c     ------------------------------------------------------------------
c--D replaces "?": ?NQSOL,?NQSL1,?ERV1,?NQFJ,?NQDOG,?NRM2,?NQFDJ
c--&               ?NQQFM,?NQQRF,?NQAQ,?NQUPD
c     Also uses IERM1, IERV1
c     ------------------------------------------------------------------
      external D1MACH, DNQFJ
      integer N, IOPT(*), IDIMW
      double precision X(N), FVEC(N), XTOL, W(IDIMW)
c
      integer IWTOLT, IWDIAG, IWA1, IWA2, IWA3, IWA4, IWGNST
      integer IWQTF, IWFJAC, IWR
      parameter(IWTOLT = 3, IWDIAG = 4 )
      double precision D1MACH, EPSFCN, EPSMCH, FAC1, FACTOR
      integer DMODE, J, JABS, K, NI, NPRINT, MAXF1, MAXFEV, ML, MU
      logical JPOS, HAVEJ, HAVED, TRACE
      parameter(FAC1 = 0.75d0,  MAXF1 = 200)
      save EPSMCH
      data EPSMCH / 0.0d0 /
c     ------------------------------------------------------------------
c
      if(EPSMCH .eq. 0.0d0) EPSMCH = D1MACH(4)
      NI = N
      IOPT(1) = 1
      if (NI .le. 0) then
         call IERM1('DNQSOL',IOPT(1),0,'Require N > 0','N',NI,'.')
         go to 900
      endif
      if (IDIMW .lt. 3 + (NI*(15+3*NI))/2) then
         call IERM1('DNQSOL',IOPT(1),0,'Require IDIMW .ge. NEED',
     *              'IDIMW',IDIMW,',')
         call IERV1('NEED =', 3 + (NI*(15+3*NI))/2,'.')
         go to 900
      endif
c                                   Set default values.
      HAVEJ = .true.
      DMODE = 1
      NPRINT = 0
      MAXFEV = MAXF1 * (NI + 1)
      ML = NI - 1
      MU = ML
      EPSFCN = EPSMCH
      FACTOR = FAC1
      TRACE = .false.
c
c                  Loop on K beginning with K = 4 and
c                  terminating when an option code, J, is zero.
      K = 4
   20 continue
      J = IOPT(K)
      JABS = abs(J)
      JPOS = J .gt. 0
      go to (40, 31, 32, 33, 34, 35, 36, 37, 38), JABS+1
c
c          ANSI Standard Fortran 77 drops thru to here if JABS is
c          larger than 7.  This is an error condition.
c
         call IERM1('DNQSOL',IOPT(1),0,'IOPT(K) must be in [-7..7]',
     *                  'K',K,',')
         call IERV1('IOPT(K)',J,'.')
         go to 900
c
   31 HAVEJ = .not. JPOS
      K = K+1
      go to 20
c                     Option 2.  Argument = 1, 2, or 3. Default = 1.
c                     1. This subr sets DIAG() to all ones.
c                     2. User has set DIAG().
c                     3. Subr DNQSL1 sets DIAG() dynamically.

   32 if( JPOS .and. IOPT(K+1) .eq. 2) then
         DMODE = 2
      elseif( JPOS .and. IOPT(K+1) .eq. 3) then
         DMODE = 3
      elseif(.not. JPOS .or. IOPT(K+1) .eq. 1) then
         DMODE = 1
      else
c                               Error.
         call IERM1('DNQSOL',IOPT(1),0,'Bad argument for Option 2.',
     *        'Argument',IOPT(K+1),'.')
         go to 900
      endif
      K = K+2
      go to 20
   33 if(JPOS) then
         NPRINT = IOPT(K+1)
      else
         NPRINT = 0
      endif
      K = K+2
      go to 20
   34 if(JPOS) then
         MAXFEV = IOPT(K+1)
      else
         MAXFEV = MAXF1 * (NI + 1)
      endif
      K = K+2
      go to 20
   35 if(JPOS) then
         ML = IOPT(K+1)
         MU = IOPT(K+2)
      else
         ML = NI+1
         MU = ML
      endif
      K = K+3
      go to 20
   36 if(JPOS) then
         EPSFCN = W(1)
      else
         EPSFCN = EPSMCH
      endif
      K = K+1
      go to 20
   37 If(JPOS) then
         FACTOR = W(2)
      else
         FACTOR = FAC1
      endif
      K = K+1
      go to 20
   38 If(JPOS) then
         TRACE = .true.
      else
         TRACE = .false.
      endif
      K = K+1
      go to 20
c                                                 End loop on K
   40 continue
c
c                     Option 2.  DMODE = 1, 2, or 3.
c                     1. This subr sets DIAG() to all ones.
c                     2. User has set DIAG().
c                     3. Subr DNQSL1 sets DIAG() dynamically.

      if(DMODE .eq. 1) then
         HAVED = .true.
         do 50 K = IWDIAG, IWDIAG+NI-1
            W(K) = 1.0d0
   50    continue
      else
         HAVED = DMODE .eq. 2
      endif
c
      IWA1 = IWDIAG + NI
      IWA2 = IWA1 + NI
      IWA3 = IWA2 + NI
      IWA4 = IWA3 + NI
      IWGNST = IWA4 + NI
      IWQTF = IWGNST + NI
      IWFJAC = IWQTF + NI
      IWR = IWFJAC + NI*NI
c     IWNEXT = IWR + (N * (N+1)) / 2    Next available loc in W().
c
      call DNQSL1(DNQFJ, NI, X, FVEC, XTOL,
     1   IOPT(1), IOPT(2), IOPT(3),
     2   NPRINT, HAVEJ, MAXFEV, HAVED, ML, MU,
     3   EPSFCN, FACTOR, TRACE, W(IWTOLT), W(IWDIAG),
     4   W(IWA1), W(IWA2), W(IWA3), W(IWA4), W(IWGNST), W(IWQTF),
     5   W(IWFJAC), W(IWR))
      return
c                             Error return
  900 continue
      IOPT(2) = 0
      IOPT(3) = 0
      W(3) = 0.0d0
      return
      end
c     ==================================================================
      subroutine DNQSL1(DNQFJ, N, X, FVEC, XTOL,
     *                 INFO, NFEV, NJEV,
     *                 NPRINT, HAVEJ, MAXFEV, HAVED, ML, MU,
     *                 EPSFCN, FACTOR, TRACE, TOLTST, DIAG,
     *                 WA1, WA2, WA3, WA4, GNSTEP, QTF, FJAC, R)
c>> 1991-12-04 CLL
c>> 1991-12-02 CLL
c>> 1991-06-18 CLL@JPL Adapting code from Minpack for MATH77

c     26 arguments.
c     Dimension of R() must be (N + N**2)/2.
c     Total space occupied by EPSFCN, FACTOR, and TOLTST through R is
c     3 + (15*N + 3*N**2)/2

      external DNQFJ
      integer N, MAXFEV, NPRINT, INFO, NFEV, NJEV, ML, MU
      logical HAVEJ, HAVED, TRACE
      double precision XTOL, EPSFCN, FACTOR, TOLTST
      double precision X(N), FVEC(N), FJAC(N,N), DIAG(N), R(*)
      double precision QTF(N), WA1(N), WA2(N), WA3(N), WA4(N)
      double precision GNSTEP(N)
C     **********
C
C     SUBROUTINE DNQSL1
C
C     THE PURPOSE OF DNQSL1 IS TO FIND A ZERO OF A SYSTEM OF
C     N NONLINEAR FUNCTIONS IN N VARIABLES BY A MODIFICATION
C     OF THE POWELL HYBRID METHOD. THE USER MUST PROVIDE A
C     SUBROUTINE WHICH CALCULATES THE FUNCTIONS and THE JACOBIAN.
C
C     ------------------------------------------------------------------
c                         Arguments
c
c   DNQFJ is THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
c     CALCULATES THE FUNCTIONS and THE JACOBIAN. DNQFJ MUST
c     BE DECLARED IN AN EXTERNAL STATEMENT IN THE USER
c     CALLING PROGRAM.  DNQFJ will not be called with IFLAG = 2
c     if HAVEJ is .false.  DNQFJ will not be called with IFLAG = 0
c     if NPRINT is <= 0.
c     DNQFJ is specified as follows:
C
c     subroutine DNQFJ(N, X, FVEC, FJAC, IFLAG)
c     integer N, IFLAG
c     double precision X(N), FVEC(N), FJAC(N,N)
c     ----------
c     if IFLAG = 0, Print X() and FVEC() and return.
c     IF IFLAG = 1 CALCULATE THE FUNCTIONS AT X AND
c     RETURN THIS VECTOR IN FVEC. DO NOT ALTER FJAC.
c     IF IFLAG = 2 CALCULATE THE JACOBIAN AT X AND
c     RETURN THIS MATRIX IN FJAC. DO NOT ALTER FVEC.
c     Set IFLAG to a negative value to force an immediate
c     termination of the solution procedure.  Otherwise do not
c     alter IFLAG.
c     ---------
c     RETURN
c     END
C
c   N is A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
c     OF FUNCTIONS and VARIABLES.
C
c   X is AN ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN
c     AN INITIAL ESTIMATE OF THE SOLUTION VECTOR. ON OUTPUT X
c     CONTAINS THE FINAL ESTIMATE OF THE SOLUTION VECTOR.
C
c   FVEC is AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
c     THE FUNCTIONS EVALUATED AT THE OUTPUT X.
C
c   XTOL is A NONNEGATIVE INPUT VARIABLE. TERMINATION
c     OCCURS WHEN THE RELATIVE ERROR BETWEEN TWO CONSECUTIVE
c     ITERATES is AT MOST XTOL.
C
c   INFO [integer,out]  If the user has terminated execution by setting
c     IFLAG negative in DNQFJ, INFO is set to IFLAG.
c     Otherwise, INFO is set as follows:
C
c     INFO = 0   Successful termination.  Radius of trust region has
c                been reduced to at most max(XTOL, machine precision).
C
c     INFO = 1   IMPROPER INPUT PARAMETERS.
C
c     INFO = 2   Number of calls to DNQFJ for function evaluations has
c                reached MAXFEV.
C
c     INFO = 3   XTOL is TOO SMALL. NO FURTHER IMPROVEMENT IN
c                THE APPROXIMATE SOLUTION X is POSSIBLE.
C
c     INFO = 4   Iteration is not making good progress, as
c                measured by the improvement through the last
c                five Jacobian evaluations.
C
c     INFO = 5   Iteration is not making good progress, as
c                measured by the improvement through last
c                ten function evaluations.
C
c   NFEV [out,integer]  The number of calls to DNQFJ with IFLAG = 1.
C
c   NJEV [out,integer] The number of evaluations of the Jacobian matrix.
c     If HAVEJ is .true. this will be the number of calls to DNQFJ with
c     IFLAG = 2.  Otherwise it is the number of times the Jacobian has
c     been approximately computed by differencing.
C
c   NPRINT [in, integer]  Enables controlled printing of iterates if it
c     is positive. In this case, DNQFJ is called with IFLAG = 0 at the
c     beginning of the first iteration and every NPRINTth time a new X
c     vector is accepted as an improvement, and at termination.
c     On these calls the new best X and FVEC are made available for
c     printing. FVEC and FJAC should not be altered.
c     If NPRINT is not positive, no special calls to DNQFJ with
c     IFLAG = 0 will be made.
C
c   HAVEJ [in, logical]  True means the user subroutine DNQFJ contains
c     code for computing the Jacobian matrix, and false means it does
c     not.
c
c   MAXFEV is A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
c     OCCURS WHEN THE NUMBER OF CALLS TO DNQFJ WITH IFLAG = 1
c     HAS REACHED MAXFEV.
C
c   HAVED  = true means initial values of DIAG() are given by the
c     calling program.  False means this subroutine must compute
c     initial values of DIAG().  It will set DIAG(j) = the euclidean
c     norm of column j, unless this is zero, in which case it will
c     set DIAG(j) = 0.0.
C
c   ML and MU specify the band structure, if any, of the Jacobian
c     matrix.  All nonzero elements of the Jacobian matrix lie
c     within the first ML subdiagonals, the main diagonal, and the
c     first MU superdiagonals.
c     ML and MU are only used when HAVEJ is .false. and are only useful
c     if ML+MU+1 < N.  In this case they are used to
c     reduce the number of function evaluations in estimating
c     derivatives.  If the Jacobian has no band structure set
c     ML = MU = N-1.
C
c   EPSFCN is an input variable used in determining a suitable
c     step length for the forward-difference approximation. This
c     approximation assumes that the relative errors in the
c     functions are of the order of max(EPSFCN, Machine precision).
C
c   FACTOR is a positive input variable used in determining the
c     initial step bound.  This bound is set to the product of
c     FACTOR and the euclidean norm of DIAG*X if nonzero, or else
c     to FACTOR itself.  In most cases FACTOR should lie in the
c     interval (0.1, 10.0).  Default: FACTOR = 0.75.
C
c   TRACE [in, logical]  If true, this subr will print detailed
c     intermediate output.  Otherwise it will not.
c
c   TOLTST  [out]  Final value of quantity that is compared with
c     XTOL for convergence test.
c
c   DIAG is an array of length N. If HAVED = false,
c     DIAG is internally set. If HAVED = true, DIAG()
c     MUST CONTAIN POSITIVE ENTRIES THAT SERVE AS
c     MULTIPLICATIVE SCALE FACTORS FOR THE VARIABLES.
C
c   WA1, WA2, WA3, and WA4 are work arrays of length N.
c
c   GNSTEP()  [scratch]  Work array of length N to save the
c     Gauss-Newton step vector computed in DNQDOG.
c
c   QTF is AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
c     THE VECTOR (Q TRANSPOSE)*FVEC.
C
c   FJAC is AN OUTPUT N BY N ARRAY WHICH CONTAINS THE
c     ORTHOGONAL MATRIX Q PRODUCED BY THE QR FACTORIZATION
c     OF THE FINAL APPROXIMATE JACOBIAN.
C
c   R is AN OUTPUT ARRAY OF LENGTH LR WHICH CONTAINS THE packed
c     UPPER TRIANGULAR MATRIX PRODUCED BY THE QR FACTORIZATION
c     OF THE FINAL APPROXIMATE JACOBIAN, STORED ROWWISE.
C     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C     SUBPROGRAMS CALLED
C
C       USER-SUPPLIED ...... DNQFJ
C
C       MINPACK-SUPPLIED ... DNQDOG,D1MACH,DNRM2,DNQFDJ,
C                            DNQQFM,DNQQRF,DNQAQ,DNQUPD
C
C       FORTRAN-SUPPLIED ... abs,max,min,mod
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE'
c     Argonne Reports: ANL-80-68 and ANL-80-74, 1980.
C
c     1991-12-09 CLL at JPL.  Replacing integer argument MODE that had
c     values 2 or 1 with logical argument HAVED related to MODE by
c     HAVED = MODE .eq. 2.  Thus the user must set HAVED = .true. when
c     supplying the DIAG() values, and .false. otherwise.
C     **********
c     ------------------------------------------------------------------
c                Description of some of the local variables.
c
c  DELTA [flpt]  Diameter of trust region.
c  HLIM0 [flpt]  Upper limit on DELTA when working with a computed
c     Jacobian.
c  HLIM1 [flpt]  Upper limit on DELTA when working with an updated
c     Jacobian.
c  JACT [integer]  Can have values of COMPUT, UPDATE, or KEEP.
c     Initially set to COMPUT.  At the beginning of the main loop
c     we either compute a new Jacobian, update the Jacobian, or keep
c     the old Jacobian, depending on the setting of JACT.
c  JACT0 [integer]  Saves the value of JACT at the beginning of the
c     main loop.  As JACT gets changed in the loop, JACT0 is still
c     available as a record of what it was at the beginning of the loop.
c  JEVAL [logical]  Set true whenever the Jacobian is computed, and set
c     false when it is updated.
c  NBEST [integer]  Counter, incremented each time an x-vector is
c     accepted as being a better approximation to the solution.  Used
c     in connection with NPRINT to trigger calles to DNQFJ for printing.
c  NCFAIL [integer]  Counts consecutive "failed" steps since the last
c     computation of the Jacobian.  NCFAIL is set to 0 when the Jacobian
c     is computed or when a step "succeeds" in the sense that
c     RATIO .ge. 0.1.  It is incremented when RATIO .lt. 0.1.
c  NLOOP [integer]  Counter for main iteration loop.
c  NUPDAT [integer]  Counts number of consecutive times the Jacobian
c     matrix is updated.
c  TRYZER [logical]  Initially set to true.  While true, the algorithm
c     will monitor X's to see if they seem to be all approaching zero.
c     If so will try setting them all to zero.  If this gives an exactly
c     zero function vector then we are finished.  If not, we set TRYZER
c     to false and restore X to its previous value (even if the function
c     value at X = 0 was an improvement) and omit any further testing
c     for X's approaching zero.  (We tryed accepting the X reached by
c     this exceptional step if the function value was an improvement,
c     but in one test case this caused the algorithm to end at a local
c     nonzero minimum rather than finding a zero.)
c     ------------------------------------------------------------------
      external D1MACH, DNRM2
      integer COMPUT, I, IFLAG, IWA(1), J, JACT, JACT0
      integer KEEP, L, LDFJAC, LR
      integer MSUM, NBEST, NCFAIL, NCSUC, NEXTPR
      integer NLOOP, NSLOW1, NSLOW2, NUMNWT, NUPDAT, UPDATE
      logical JEVAL, NEWX, NEWTOK, SING, TRYZER
      double precision D1MACH,DNRM2
      double precision ACTRED,DELTA,EPSMCH,FNORM,FNORM1, HLIM0, HLIM1
      double precision ONE,PNORM, PRERED,P1,P5,P0001,RATIO
      double precision SUM,TEMP,XNORM, ZERO
      parameter(COMPUT = 1, UPDATE = 2, KEEP = 3)
      parameter(ONE = 1.0d0, P1 = 0.1d0, P5 = 0.5d0)
      parameter(P0001 = 0.0001d0, ZERO = 0.0d0)
      save EPSMCH
      data EPSMCH /0.0d0 /
c     ------------------------------------------------------------------
C                     Set EPSMCH to the machine precision.
C
      if(EPSMCH .eq. 0.0d0) EPSMCH = D1MACH(4)
C
C                  Initialize values of output arguments.
C
      INFO = 1
      NFEV = 0
      NJEV = 0
      TOLTST = 0.0d0
      TRYZER = .true.
C
C        CHECK THE INPUT PARAMETERS FOR ERRORS.
C        We assume the condition N > 0 has already been checked in
c        the user-interface subroutine that called this one.

      IF ( XTOL .lt. ZERO .or. MAXFEV .le. 0
     *   .or. FACTOR .le. ZERO ) then
           call IERM1('DNQSL1',INFO,0,
     *      'Require MAXFEV > 0, XTOL .gt. 0.0, FACTOR > 0.0',
     *      'MAXFEV',MAXFEV,',')
           call DERV1('XTOL',XTOL,',')
           call DERV1('FACTOR',FACTOR,'.')
           go to 300
      endif
      if( .not. HAVEJ .and. (ML .lt. 0 .or. MU .lt. 0)) then
           call IERM1('DNQSL1',INFO,0,
     *      'With HAVEJ false, require ML .ge. 0 and MU .ge. 0',
     *              'ML',ML,',')
         call IERV1('MU',MU,',')
         go to 300
      endif
c                            HAVED = true means the user has set DIAG().
      IF ( HAVED ) then
         DO 10 J = 1, N
            IF (DIAG(J) .le. ZERO) then
                 call IERM1('DNQSL1',INFO,0,
     *         'With HAVED = .true., require all DIAG(J) > 0.0',
     *                    'J',J,',')
                 call DERV1('DIAG(J)',DIAG(J),'.')
               go to 300
            endif
   10    CONTINUE
      endif
c                               Initialize algorithm variables.
      INFO = 0
      JACT   = COMPUT
      LDFJAC = N
      LR     = (N*(N+1)) / 2
      MSUM   = min(ML + MU + 1, N)
      NBEST  = 1
      NCSUC  = 0
      NEXTPR = 1
      NLOOP  = 0
      NSLOW1 = 0
      NSLOW2 = 0
      NUMNWT = 0
C
C                  Evaluate the function at the starting point.
C                  Calculate and test its norm.
C
      IFLAG = 1
C%%     (*dnqfj)( n, x, fvec, fjac, &iflag );
      CALL DNQFJ(N, X, FVEC, FJAC, IFLAG)
      NFEV = 1
      IF (IFLAG .lt. 0) GO TO 300
      FNORM = DNRM2(N,FVEC,1)
      if(TRACE) then
         print'(1x,i5,a/(6x,5g15.6))',NLOOP,
     *        ' Initial X:',(X(J),J=1,N)
         print'(1x,5x,a,g15.6)',
     *      ' Initial FNORM:',FNORM
      endif
      if(FNORM .eq. 0.0d0) then
         go to 300
      endif
C
C                                     Beginning of main loop.
C
   30 continue
         NLOOP = NLOOP + 1
         JACT0 = JACT
C
C              Compute, Update, or Keep Jacobian, depending on JACT.
C
      if (JACT .eq. COMPUT) then
         JEVAL = .TRUE.
         NUPDAT = 0
         NCFAIL = 0
C
C        CALCULATE THE JACOBIAN MATRIX.
C
         if(TRACE) print'(1x,i5,a)',NLOOP,
     *       ' Computing new Jacobian matrix.'
         NJEV = NJEV + 1
         if(HAVEJ) then
            IFLAG = 2
C%%           (*dnqfj)( n, x, fvec, fjac, &iflag );
            CALL DNQFJ(N, X, FVEC, FJAC, IFLAG)
         else
            CALL DNQFDJ(DNQFJ,N,X,FVEC,FJAC,LDFJAC,
     *                  IFLAG,ML,MU,EPSFCN,WA1, WA2)
            NFEV = NFEV + MSUM
         endif
         IF (IFLAG .lt. 0) GO TO 300
C
C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
C
         CALL DNQQRF(N,N,FJAC,LDFJAC, .false., IWA,1,WA1,WA2,WA3)
C
C        On the first iteration and if HAVED is .false., scale according
C        to the norms of the columns of the initial Jacobian.
C        Also on the first iteration calculate the norm of the scaled X
C        and initialize the trust region diameter, DELTA.
C
         IF (NLOOP .eq. 1) then
            IF ( .not. HAVED ) then
               DO 40 J = 1, N
                  DIAG(J) = WA2(J)
                  IF (WA2(J) .eq. ZERO) DIAG(J) = ONE
   40          CONTINUE
            endif
C
            DO 60 J = 1, N
               WA3(J) = DIAG(J)*X(J)
   60       CONTINUE
            XNORM = DNRM2(N,WA3,1)
            DELTA = FACTOR*XNORM
            IF (DELTA .eq. ZERO) DELTA = FACTOR
         endif
C
C        FORM (Q TRANSPOSE)*FVEC and STORE IN QTF.
C
         DO 80 I = 1, N
            QTF(I) = FVEC(I)
   80       CONTINUE
         DO 120 J = 1, N
            IF (FJAC(J,J) .eq. ZERO) GO TO 110
            SUM = ZERO
            DO 90 I = J, N
               SUM = SUM + FJAC(I,J)*QTF(I)
   90          CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 100 I = J, N
               QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  100          CONTINUE
  110       CONTINUE
  120       CONTINUE
C
C        COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
c        The diagonal elts come from WA1().  The strictly upper
c        triangular elts come from FJAC(,).  The upper triangular matrix
c        will be stored, packed by rows, in R().
C
         SING = .FALSE.
         DO 150 J = 1, N
            L = J
            DO 130 I = 1, J-1
               R(L) = FJAC(I,J)
               L = L + N - I
  130          CONTINUE
            R(L) = WA1(J)
            IF (WA1(J) .eq. ZERO) SING = .true.
  150       CONTINUE
C
C        ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
C
         CALL DNQQFM(N,N,FJAC,LDFJAC,WA1)
C
C        RESCALE IF NECESSARY.
C
         if ( .not. HAVED ) then
            DO 160 J = 1, N
               DIAG(J) = max(DIAG(J),WA2(J))
  160       CONTINUE
         endif
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elseif(JACT .eq. UPDATE) then

C
C           CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
C           and UPDATE QTF IF NECESSARY.
C
         if(TRACE) print'(1x,i5,a)',NLOOP,
     *       ' Updating Jacobian matrix.'
            NUPDAT = NUPDAT + 1
            JEVAL = .FALSE.
            DO 280 J = 1, N
               SUM = ZERO
               DO 270 I = 1, N
                  SUM = SUM + FJAC(I,J)*WA4(I)
  270             CONTINUE
               WA2(J) = (SUM - WA3(J))/PNORM
               WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
               IF (RATIO .ge. P0001) QTF(J) = SUM
  280          CONTINUE
C
C           COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
C
            CALL DNQUPD(N,N,R,LR,WA1,WA2,WA3,SING)
            CALL DNQAQ(N,N,FJAC,LDFJAC,WA2,WA3)
            CALL DNQAQ(1,N,QTF,1,WA2,WA3)
      else
            if(TRACE) print'(1x,i5,a)',NLOOP,
     *       ' Keeping Jacobian matrix unchanged.'
      endif
C
C           Now have a new or updated or retained Jacobian matrix.
C
C           IF REQUESTED, CALL DNQFJ TO ENABLE PRINTING OF ITERATES.
C
            if (NPRINT .gt. 0) then
               if (NBEST .eq. NEXTPR) then
                  IFLAG = 0
C%%               (*dnqfj)( n, x, fvec, fjac, &iflag );
                  CALL DNQFJ(N, X, FVEC, FJAC, IFLAG)
                  IF (IFLAG .lt. 0) GO TO 300
                  NEXTPR = NEXTPR + NPRINT
               endif
            endif
C
C           Determine the direction P, using a dogleg method, and
c           returning -P in WA1().
C
            CALL DNQDOG(N,R,LR,DIAG,QTF,DELTA,WA1,NEWTOK,WA2,WA3,
     *                  JACT0 .eq. KEEP, GNSTEP)
c
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if(TRYZER) then
c                                  NUMNWT counts number of consecutive
c                                  full Newton steps.
            if(NEWTOK) then
               NUMNWT = NUMNWT + 1
            else
               NUMNWT = 0
            endif
c
c              Test for convergence of some x components toward 0.
c              If this seems to be happening try setting such
c              components to 0.
c
            if(NUMNWT .ge. 5 .and. NCSUC .ge. 4) then
               NUMNWT = 0
               do 204 J = 1,N
                  WA2(J) = X(J) - WA1(J)
                  if(abs(WA2(J)) .le. 0.75d0 * abs(X(J)) ) then
                     WA2(J) = 0.0d0
                  else
                     go to 203
                  endif
  204          continue
               if(TRACE) print'(1x,i5,a)',NLOOP,
     *            ' Trial setting of X() to zero.'
C
C              EVALUATE THE FUNCTION AT WA2() and CALCULATE ITS NORM.
C
                  IFLAG = 1
c%%                  (*dnqfj)( n, wa2, wa4, fjac, &iflag );
                  CALL DNQFJ(N, WA2, WA4, FJAC, IFLAG)
                  NFEV = NFEV + 1
                  IF (IFLAG .lt. 0) GO TO 300
                  FNORM1 = DNRM2(N,WA4,1)
                  if(TRACE) print'(1x,i5,a,g15.6)',NLOOP,
     *             ' FNORM1 =      ', FNORM1
                  if(FNORM1 .eq. 0.0d0) then
c
C                    Accept new point as final solution.
c                    Update X() and FVEC() and go to termination.
C
                     INFO = 0
                     TOLTST = 0.0d0
                     do 201 J = 1, N
                        X(J) = WA2(J)
                        FVEC(J) = WA4(J)
  201                continue
                     if(TRACE) print'(1x,i5,a,(6x,5g15.6))',NLOOP,
     *                  ' Accepting X = all zeros.'
                     go to 300
                  else
                     TRYZER = .false.
                  endif
            endif
c                      The following "endif" matches "if(TRYZER)then"
            endif
c     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C           STORE THE DIRECTION P and X + P. CALCULATE THE NORM OF P.
C
  203       continue
            DO 200 J = 1, N
               WA1(J) = -WA1(J)
               WA2(J) = X(J) + WA1(J)
               WA3(J) = DIAG(J)*WA1(J)
  200       continue

            PNORM = DNRM2(N,WA3,1)
            if(TRACE) then
               print'(1x,i5,a,/1x,5x,2g15.6)',NLOOP,
     *         '    DELTA          PNORM',DELTA,PNORM
               print'(6x,a/(6x,5g15.6))',' Trial X:',(WA2(J),J=1,N)
            endif
C
C           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
C
            IF (NLOOP .eq. 1) then
               DELTA = min(DELTA,PNORM)
               HLIM0 = DELTA
               HLIM1 = DELTA
            endif
C
C           EVALUATE THE FUNCTION AT X + P and CALCULATE ITS NORM.
C
            IFLAG = 1
c%%           (*dnqfj)( n, wa2, wa4, fjac, &iflag );
            CALL DNQFJ(N, WA2, WA4, FJAC, IFLAG)
            NFEV = NFEV + 1
            IF (IFLAG .lt. 0) GO TO 300
            FNORM1 = DNRM2(N,WA4,1)
C
C           COMPUTE THE SCALED ACTUAL REDUCTION.
C
            ACTRED = -ONE
            IF (FNORM1 .lt. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C           COMPUTE THE SCALED PREDICTED REDUCTION.
C
            L = 1
            DO 220 I = 1, N
               SUM = ZERO
               DO 210 J = I, N
                  SUM = SUM + R(L)*WA1(J)
                  L = L + 1
  210             CONTINUE
               WA3(I) = QTF(I) + SUM
  220          CONTINUE
            TEMP = DNRM2(N,WA3,1)
            PRERED = ZERO
            IF (TEMP .lt. FNORM) PRERED = ONE - (TEMP/FNORM)**2
C
C           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
C           REDUCTION.
C
            RATIO = ZERO
            IF (PRERED .gt. ZERO) RATIO = ACTRED/PRERED
            if(TRACE) print'(1x,i5,a,/1x,5x,4g15.6)',NLOOP,
     *      '    FNORM1        ACTRED          PRERED        RATIO',
     *      FNORM1,ACTRED,PRERED,RATIO
C
c           Analyze RATIO, NCSUC and JEVAL to decide on accepting or
c           rejecting the new X, and assigning new values to
c           NCSUC, JACT, and DELTA.
c
            if( RATIO .lt. 0.0000D0) then
               NCSUC = 0
               NCFAIL = NCFAIL + 1
               NEWX = .false.
               if(JEVAL) HLIM0 = min(HLIM0, 0.707107d0 * PNORM)
               HLIM1 = min(HLIM1, 0.707107d0 * PNORM)
               if( JEVAL .or. (NCFAIL .le. 1 .and. NUPDAT .le. 2)) then
                  JACT = KEEP
                  DELTA = 0.5d0 * PNORM
               else
                  JACT = COMPUT
                  DELTA =  HLIM0
               endif
            elseif( RATIO .lt. 0.1D0) then
               NCSUC = 0
               NCFAIL = NCFAIL + 1
               NEWX = .true.
               if(NCFAIL .le. 1 .and. NUPDAT .le. 2) then
                  JACT = UPDATE
                  DELTA = 0.5d0 * PNORM
               else
                  JACT = COMPUT
                  DELTA = HLIM0
               endif
            else
c                                     Here we have RATIO .ge. 0.1
               NCSUC = NCSUC + 1
               NCFAIL = 0
               NEWX = .true.
               JACT = UPDATE
               if(RATIO .lt. 0.5d0) then
                  if(NCSUC .ge. 5)
     *               HLIM1 = max(HLIM1, 1.414214d0 * PNORM)
                  if(NCSUC .ge. 2)
     *               DELTA = min(HLIM1, max(DELTA, 1.414214d0 * PNORM))
               elseif(RATIO .lt. 0.9d0) then
                  if(JACT0 .eq. COMPUT)
     *               HLIM0 = max(HLIM0, 1.414214d0 * PNORM)
                  if(NCSUC .ge. 4)
     *               HLIM1 = max(HLIM1, 1.414214d0 * PNORM)
                  if(NCSUC .ge. 2)
     *               DELTA = min(HLIM1, max(DELTA, 1.414214d0 * PNORM))
               elseif(RATIO .lt. 1.1d0) then
                  if(JACT0 .eq. COMPUT)
     *               HLIM0 = max(HLIM0, 2.0d0 * PNORM)
                  if(NCSUC .eq. 1) then
                     DELTA = 1.414214d0 * PNORM
                  else
                     DELTA =  2.0d0 * PNORM
                  endif
                  HLIM1 = max(HLIM1, DELTA)
               endif
            endif
            HLIM0 = max(HLIM0, HLIM1)
            if(TRACE) print'(1x,i5,a,a,/1x,5x,3i8,3g13.4)',NLOOP,
     *      '     NCSUC  NCFAIL  NUPDAT',
     *      ' DELTA        HLIM0        HLIM1',
     *             NCSUC, NCFAIL,  NUPDAT,  DELTA,HLIM0,   HLIM1
C
            if(NEWX) then
c                               Accept new X, FVEC, and their norms.
               DO 250 J = 1, N
                  X(J) = WA2(J)
                  WA2(J) = DIAG(J)*X(J)
                  FVEC(J) = WA4(J)
  250          CONTINUE
               XNORM = DNRM2(N,WA2,1)
               FNORM = FNORM1
               NBEST = NBEST + 1
               if(TRACE) print'(1x,i5,a,g15.6)',NLOOP,
     *         ' Accepting new X with XNORM =  ',XNORM
            endif
C
C                        DETERMINE THE PROGRESS OF THE ITERATION.
C
            if( ACTRED .ge. 0.001d0) then
               NSLOW1 = 0
            else
               NSLOW1 = NSLOW1 + 1
            endif
            if( ACTRED .ge. 0.1d0) then
               NSLOW2 = 0
            elseif( JACT0 .eq. COMPUT) then
               NSLOW2 = NSLOW2 + 1
            endif
            if(TRACE) print'(1x,i5,a,/1x,5x,2(i11,4x))',NLOOP,
     *      '     NSLOW1         NSLOW2',
     *      NSLOW1,       NSLOW2
C
C                           TEST FOR CONVERGENCE.
C
            IF (DELTA .le. XTOL*XNORM .or. FNORM .eq. ZERO) then
                INFO = 0
               if(TRACE) print'(1x,i5,a,/1x,5x,i14,g15.6)',NLOOP,
     *         '          INFO   XNORM', INFO, XNORM
                go to 295
            endif
C
C                    TESTS FOR TERMINATION and STRINGENT TOLERANCES.
C
            IF (NFEV .ge. MAXFEV) INFO = 2
            IF (P1*max(P1*DELTA,PNORM) .le. EPSMCH*XNORM) INFO = 3
            IF (NSLOW2 .eq. 5) INFO = 4
            IF (NSLOW1 .eq. 10) INFO = 5
            IF (INFO .ne. 0) then
               if(TRACE) print'(1x,i5,a,/1x,5x,i14,g15.6)',NLOOP,
     *         '          INFO   XNORM', INFO, XNORM
               call IERM1('DNQSL1',INFO, 0,'Unsuccessful termination.',
     *         'INFO',INFO,'.')
               go to 295
            endif
      go to 30
C                                              End of main loop.
C
c                   Come to following stmt when INFO has been set to
c                   2, 3, 4, or 5, or to 0 due to successful XTOL test.
  295 continue
               if(XNORM .ne. 0.0d0) then
                  TOLTST = DELTA / XNORM
               else
                  TOLTST = DELTA
               endif
c
c                    Jump to following statement with IFLAG negative
c                    or INFO = 1 or INFO  = 0 due to FNORM being zero.
c                    Here we have TOLTST = 0.0.
  300 continue
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
      IF (IFLAG .lt. 0) INFO = IFLAG
      if(TRACE) print'(1x,i5,a,i3)',NLOOP,
     *      ' Quitting with INFO = ',INFO
      IFLAG = 0
c%%   if (nprint > 0) (*dnqfj)( n, x, fvec, fjac, &iflag );
      IF (NPRINT .gt. 0) CALL DNQFJ(N,X,FVEC,FJAC, IFLAG)
      if(INFO .lt. 0) then
         call IERM1('DNQSL1',INFO, 0,
     *      'Quitting because user code set IFLAG negative.',
     *      'IFLAG',INFO,'.')
      endif
      return
C
C     Last line of subroutine DNQSL1.
C
      END
c     ==================================================================
      subroutine DNQFDJ(DNQFJ,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
     *                  WA1,WA2)
c>> 1991-12-04 CLL  Changed arg list of user supplied subroutine.
c>> 1991-06-18 CLL@JPL Adapting code from Minpack for MATH77
      external DNQFJ
      integer N,LDFJAC,IFLAG,ML,MU
      double precision EPSFCN
      double precision X(N),FVEC(N),FJAC(LDFJAC,N),WA1(N),WA2(N)
C     **********
C
C     SUBROUTINE DNQFDJ
C
C     THIS SUBROUTINE COMPUTES A FORWARD-DIFFERENCE APPROXIMATION
C     TO THE N BY N JACOBIAN MATRIX ASSOCIATED WITH A SPECIFIED
C     PROBLEM OF N FUNCTIONS IN N VARIABLES. IF THE JACOBIAN HAS
C     A BANDED FORM, THEN FUNCTION EVALUATIONS ARE SAVED BY ONLY
C     APPROXIMATING THE NONZERO TERMS.
C
C     THE SUBROUTINE STATEMENT IS
C
C     subroutine DNQFDJ(DNQFJ,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
C                         WA1,WA2)
C
C     WHERE
C
C       DNQFJ IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
C         CALCULATES THE FUNCTIONS. DNQFJ MUST BE DECLARED
C         IN AN EXTERNAL STATEMENT IN THE USER CALLING
C         PROGRAM, and SHOULD BE WRITTEN AS FOLLOWS.
C
C         subroutine DNQFJ(N,X,FVEC,IFLAG)
C         integer N,IFLAG
C         double precision X(N),FVEC(N)
C         ----------
C         CALCULATE THE FUNCTIONS AT X AND
C         RETURN THIS VECTOR IN FVEC.
C         ----------
C         RETURN
C         END
C
C         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY DNQFJ UNLESS
C         THE USER WANTS TO TERMINATE EXECUTION OF DNQFDJ.
C         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF FUNCTIONS and VARIABLES.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       FVEC IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
C         FUNCTIONS EVALUATED AT X.
C
C       FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE
C         APPROXIMATION TO THE JACOBIAN MATRIX EVALUATED AT X.
C
C       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
C
C       IFLAG IS AN INTEGER VARIABLE WHICH CAN BE USED TO TERMINATE
C         THE EXECUTION OF DNQFDJ. SEE DESCRIPTION OF DNQFJ.
C
C       ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
C         THE NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE
C         JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
C         ML TO AT LEAST N - 1.
C
C       MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
C         THE NUMBER OF SUPERDIAGONALS WITHIN THE BAND OF THE
C         JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
C         MU TO AT LEAST N - 1.
C
C       EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
C         STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
C         APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
C         FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
C         THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
C         ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
C         PRECISION.
C
C       WA1 and WA2 ARE WORK ARRAYS OF LENGTH N. IF ML + MU + 1 IS AT
C         LEAST N, THEN THE JACOBIAN IS CONSIDERED DENSE, and WA2 IS
C         NOT REFERENCED.
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... D1MACH
C
C       FORTRAN-SUPPLIED ... abs,max,sqrt
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
c     ------------------------------------------------------------------
      external D1MACH
      integer I,J,K,MSUM
      double precision  EPS,EPSMCH,H,TEMP,ZERO
c++ CODE for ~.C. is active
      double precision DUMMY(1,1)
c++ CODE for .C. & (.N. == 'S') is inactive
c%%     float *dummy;
c++ CODE for .C. & (.N. == 'D') is inactive
c%%     double *dummy;
C++ End
      double precision D1MACH
      parameter(ZERO = 0.0d0)
C
C     EPSMCH IS THE MACHINE PRECISION.
C
      EPSMCH = D1MACH(4)
C
      EPS = sqrt(max(EPSFCN,EPSMCH))
      IFLAG = 1
      MSUM = ML + MU + 1
      IF (MSUM .lt. N) GO TO 40
C
C        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
C
         DO 20 J = 1, N
            TEMP = X(J)
            H = EPS*abs(TEMP)
            IF (H .eq. ZERO) H = EPS
            X(J) = TEMP + H
c%%           (*dnqfj)( n, x, wa1, dummy, iflag );
            CALL DNQFJ(N, X, WA1, DUMMY, IFLAG)
            IF (IFLAG .lt. 0) GO TO 30
            X(J) = TEMP
            DO 10 I = 1, N
               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
         GO TO 110
   40 CONTINUE
C
C        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
C
         DO 90 K = 1, MSUM
            DO 60 J = K, N, MSUM
               WA2(J) = X(J)
               H = EPS*abs(WA2(J))
               IF (H .eq. ZERO) H = EPS
               X(J) = WA2(J) + H
   60          CONTINUE
c%%           (*dnqfj)( n, x, wa1, dummy, iflag );
            CALL DNQFJ(N, X, WA1, DUMMY, IFLAG)
            IF (IFLAG .lt. 0) GO TO 100
            DO 80 J = K, N, MSUM
               X(J) = WA2(J)
               H = EPS*abs(WA2(J))
               IF (H .eq. ZERO) H = EPS
               DO 70 I = 1, N
                  FJAC(I,J) = ZERO
                  IF (I .ge. J - MU .and. I .le. J + ML)
     *               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
      RETURN
C
C     Last line of subroutine DNQFDJ.
C
      END
c     ==================================================================
      subroutine DNQAQ(M,N,A,LDA,V,W)
      integer M,N,LDA
      double precision A(LDA,N),V(N),W(N)
C     **********
C
C     SUBROUTINE DNQAQ
C
C     GIVEN AN M BY N MATRIX A, THIS SUBROUTINE COMPUTES A*Q WHERE
C     Q IS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     and GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE WHICH
C     ELIMINATE ELEMENTS IN THE I-TH and N-TH PLANES, RESPECTIVELY.
C     Q ITSELF IS NOT GIVEN, RATHER THE INFORMATION TO RECOVER THE
C     GV, GW ROTATIONS IS SUPPLIED.
C
C     THE SUBROUTINE STATEMENT IS
C
C       subroutine DNQAQ(M,N,A,LDA,V,W)
C
C     WHERE
C
C       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF ROWS OF A.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF COLUMNS OF A.
C
C       A IS AN M BY N ARRAY. ON INPUT A MUST CONTAIN THE MATRIX
C         TO BE POSTMULTIPLIED BY THE ORTHOGONAL MATRIX Q
C         DESCRIBED ABOVE. ON OUTPUT A*Q HAS REPLACED A.
C
C       LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.
C
C       V IS AN INPUT ARRAY OF LENGTH N. V(I) MUST CONTAIN THE
C         INFORMATION NECESSARY TO RECOVER THE GIVENS ROTATION GV(I)
C         DESCRIBED ABOVE.
C
C       W IS AN INPUT ARRAY OF LENGTH N. W(I) MUST CONTAIN THE
C         INFORMATION NECESSARY TO RECOVER THE GIVENS ROTATION GW(I)
C         DESCRIBED ABOVE.
C
C     SUBROUTINES CALLED
C
C       FORTRAN-SUPPLIED ... abs,sqrt
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
c     ------------------------------------------------------------------
      integer I,J,NMJ,NM1
      double precision VCOS,ONE,VSIN,TEMP
      parameter(ONE = 1.0d0)
C
C     APPLY THE FIRST SET OF GIVENS ROTATIONS TO A.
C
      NM1 = N - 1
      IF (NM1 .lt. 1) GO TO 50
      DO 20 NMJ = 1, NM1
         J = N - NMJ
         IF (abs(V(J)) .gt. ONE) VCOS = ONE/V(J)
         IF (abs(V(J)) .gt. ONE) VSIN = sqrt(ONE-VCOS**2)
         IF (abs(V(J)) .le. ONE) VSIN = V(J)
         IF (abs(V(J)) .le. ONE) VCOS = sqrt(ONE-VSIN**2)
         DO 10 I = 1, M
            TEMP = VCOS*A(I,J) - VSIN*A(I,N)
            A(I,N) = VSIN*A(I,J) + VCOS*A(I,N)
            A(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
C
C     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
C
      DO 40 J = 1, NM1
         IF (abs(W(J)) .gt. ONE) VCOS = ONE/W(J)
         IF (abs(W(J)) .gt. ONE) VSIN = sqrt(ONE-VCOS**2)
         IF (abs(W(J)) .le. ONE) VSIN = W(J)
         IF (abs(W(J)) .le. ONE) VCOS = sqrt(ONE-VSIN**2)
         DO 30 I = 1, M
            TEMP = VCOS*A(I,J) + VSIN*A(I,N)
            A(I,N) = -VSIN*A(I,J) + VCOS*A(I,N)
            A(I,J) = TEMP
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
C
C     Last line of subroutine DNQAQ.
C
      END
c     ==================================================================
      subroutine DNQDOG(N,R,LR,DIAG,QTB,DELTA,X,NEWTOK,WA1,WA2,
     *                  SAMEJ, GNSTEP)
c>> 1992-01-03 CLL
      integer N,LR
      logical SAMEJ, NEWTOK
      double precision DELTA, GNSTEP(N)
      double precision R(LR),DIAG(N),QTB(N),X(N),WA1(N),WA2(N)
C     **********
C
C     subroutine DNQDOG
C
c     Given an M by N matrix A, an N by N nonsingular diagonal
c     matrix D, an M-vector B, and a positive number DELTA, the
c     problem is to determine the convex combination X of the
c     gauss-newton and scaled gradient directions that minimizes
c     (A*X - B) in the least squares sense, subject to the
c     restriction that the euclidean norm of D*X be at most DELTA.
c
c     This subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     QR factorization of A. that is, if A = Q*R, where Q has
c     orthogonal columns and R is an upper triangular matrix,
c     then DNQDOG needs the full upper triangle of R and
c     the first N components of (Q transpose)*B.
c
c     The subroutine statement is
C
C       subroutine DNQDOG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
C
C     where
C
c  N is A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.
C
c  R() [in]  An ARRAY OF LENGTH LR WHICH MUST CONTAIN THE UPPER
c    TRIANGULAR MATRIX R STORED BY ROWS.
C
c  LR is A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
c    (N*(N+1))/2.
C
c  DIAG() [in]  An ARRAY OF LENGTH N WHICH MUST CONTAIN THE
c    DIAGONAL ELEMENTS OF THE MATRIX D.
C
c  QTB() [in]  An ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
c    N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.
C
c  DELTA is a POSITIVE INPUT VARIABLE WHICH SPECIFIES AN UPPER
c    BOUND ON THE EUCLIDEAN NORM OF D*X.
C
c  X() [out]  An ARRAY OF LENGTH N WHICH CONTAINS THE DESIRED
c    CONVEX COMBINATION OF THE GAUSS-NEWTON DIRECTION and THE
c    SCALED GRADIENT DIRECTION.
c
c  NEWTOK [logical, out]  True means the full Newton step was
c    used.  False means a modified, shorter, step was used.
c
c  WA1() and WA2() are work arrays of length N.
c
c  SAMEJ [logical, in]  True means we have the same Jacobian matrix as
c    on the previous call to this subr.  The Gauss-Newton vector in
c    GNSTEP() can be reused.
c
c  GNSTEP() [inout]  On return holds the Gauss-Newton vector.  On entry
c    with SAMEJ = .true., contains the GN vector from the previous call.
C     ------------------------------------------------------------------
C       SUBPROGRAMS CALLED
C
c       MINPACK-SUPPLIED ... D1MACH,DNRM2
C
C       FORTRAN-SUPPLIED ... abs,max,min,sqrt
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
c     ------------------------------------------------------------------
      external D1MACH, DNRM2
      integer I,J,JJ,JP1,K,L, NP1
      double precision ALPHA,BNORM,EPSMCH,GNORM,ONE,QNORM,SGNORM,SUM
      double precision TEMP,ZERO
      double precision D1MACH,DNRM2
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
      save EPSMCH
      data EPSMCH / 0.0d0 /
c     ------------------------------------------------------------------
C                   Set EPSMCH to the machine precision.
C
      if(EPSMCH .eq. 0.0d0) EPSMCH = D1MACH(4)
      if(.not. SAMEJ) then
C
C     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.
C
      NP1 = N+1
      JJ = (N*(N + 1))/2 + 1
      DO 50 K = 1, N
         J = NP1 - K
         JP1 = J + 1
         JJ = JJ - K
         L = JJ + 1
         SUM = ZERO
            DO 10 I = JP1, N
               SUM = SUM + R(L)*GNSTEP(I)
               L = L + 1
   10       CONTINUE
         TEMP = R(JJ)
         IF (TEMP .eq. ZERO) then
            L = J
            DO 30 I = 1, J-1
               TEMP = max(TEMP,abs(R(L)))
               L = L + N - I
   30       CONTINUE
            TEMP = EPSMCH*TEMP
         endif
         if (TEMP .eq. ZERO) then
            GNSTEP(J) = 0.0d0
         else
            GNSTEP(J) = (QTB(J) - SUM)/TEMP
         endif
   50    CONTINUE
      endif
C
C     TEST WHETHER THE GAUSS-NEWTON DIRECTION is ACCEPTABLE.
C
      DO 60 J = 1, N
         WA1(J) = ZERO
         WA2(J) = DIAG(J)*GNSTEP(J)
   60    CONTINUE
      QNORM = DNRM2(N,WA2,1)
      NEWTOK = QNORM .le. DELTA
      if (NEWTOK) then
         do 65 J = 1,N
            X(J) = GNSTEP(J)
   65    continue
         go to 140
      endif
C
C     THE GAUSS-NEWTON DIRECTION is NOT ACCEPTABLE.
C     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
C
      L = 1
      DO 80 J = 1, N
         TEMP = QTB(J)
         DO 70 I = J, N
            WA1(I) = WA1(I) + R(L)*TEMP
            L = L + 1
   70       CONTINUE
         WA1(J) = WA1(J)/DIAG(J)
   80    CONTINUE
C
C     CALCULATE THE NORM OF THE SCALED GRADIENT and TEST FOR
C     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT is ZERO.
C
      GNORM = DNRM2(N,WA1,1)
      if (GNORM .eq. ZERO) then
         ALPHA = DELTA/QNORM
         do 85 J = 1, N
            X(J) = ALPHA*GNSTEP(J)
   85    continue
         go to 140
      endif
C
C     CALCULATE THE POINT ALONG THE SCALED GRADIENT
C     AT WHICH THE QUADRATIC is MINIMIZED.
C
      DO 90 J = 1, N
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   90    CONTINUE
      L = 1
      DO 110 J = 1, N
         SUM = ZERO
         DO 100 I = J, N
            SUM = SUM + R(L)*WA1(I)
            L = L + 1
  100       CONTINUE
         WA2(J) = SUM
  110    CONTINUE
      TEMP = DNRM2(N,WA2,1)
      SGNORM = (GNORM/TEMP)/TEMP
C
C     TEST WHETHER THE SCALED GRADIENT DIRECTION is ACCEPTABLE.
C
      ALPHA = ZERO
      if (SGNORM .lt. DELTA) then
C
C     THE SCALED GRADIENT DIRECTION is NOT ACCEPTABLE.
C     FINALLY, CALCULATE THE POINT ALONG THE dogleg
C     AT WHICH THE QUADRATIC is MINIMIZED.
C
      BNORM = DNRM2(N,QTB,1)
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
      TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2
     *       + sqrt((TEMP-(DELTA/QNORM))**2
     *               +(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
      ALPHA = ((DELTA/QNORM)*(ONE - (SGNORM/DELTA)**2))/TEMP
      endif
C
C     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
C     DIRECTION and THE SCALED GRADIENT DIRECTION.
C
      TEMP = (ONE - ALPHA)*min(SGNORM,DELTA)
      DO 130 J = 1, N
         X(J) = TEMP*WA1(J) + ALPHA*GNSTEP(J)
  130    CONTINUE
  140 CONTINUE
      RETURN
C
C     Last line of subroutine DNQDOG.
C
      END
      subroutine DNQQFM(M,N,Q,LDQ,WA)
      integer M,N,LDQ
      double precision Q(LDQ,M),WA(M)
C     **********
C
C     SUBROUTINE DNQQFM
C
C     THIS SUBROUTINE PROCEEDS FROM THE COMPUTED QR FACTORIZATION OF
C     AN M BY N MATRIX A TO ACCUMULATE THE M BY M ORTHOGONAL MATRIX
C     Q FROM ITS FACTORED FORM.
C
C     THE SUBROUTINE STATEMENT IS
C
C       subroutine DNQQFM(M,N,Q,LDQ,WA)
C
C     WHERE
C
C       M is A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF ROWS OF A and THE ORDER OF Q.
C
C       N is A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF COLUMNS OF A.
C
C       Q is AN M BY M ARRAY. ON INPUT THE FULL LOWER TRAPEZOID IN
C         THE FIRST MIN(M,N) COLUMNS OF Q CONTAINS THE FACTORED FORM.
C         ON OUTPUT Q HAS BEEN ACCUMULATED INTO A SQUARE MATRIX.
C
C       LDQ is A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY Q.
C
C       WA is A WORK ARRAY OF LENGTH M.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... min
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
c     ------------------------------------------------------------------
      integer I,J,JM1,K,L,MINMN,NP1
      double precision ONE,SUM,TEMP,ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
C
C     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.
C
      MINMN = min(M,N)
      IF (MINMN .lt. 2) GO TO 30
      DO 20 J = 2, MINMN
         JM1 = J - 1
         DO 10 I = 1, JM1
            Q(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
C     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
C
      NP1 = N + 1
      IF (M .lt. NP1) GO TO 60
      DO 50 J = NP1, M
         DO 40 I = 1, M
            Q(I,J) = ZERO
   40       CONTINUE
         Q(J,J) = ONE
   50    CONTINUE
   60 CONTINUE
C
C     ACCUMULATE Q FROM ITS FACTORED FORM.
C
      DO 120 L = 1, MINMN
         K = MINMN - L + 1
         DO 70 I = K, M
            WA(I) = Q(I,K)
            Q(I,K) = ZERO
   70       CONTINUE
         Q(K,K) = ONE
         IF (WA(K) .eq. ZERO) GO TO 110
         DO 100 J = K, M
            SUM = ZERO
            DO 80 I = K, M
               SUM = SUM + Q(I,J)*WA(I)
   80          CONTINUE
            TEMP = SUM/WA(K)
            DO 90 I = K, M
               Q(I,J) = Q(I,J) - TEMP*WA(I)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      RETURN
C
C     Last line of subroutine DNQQFM.
C
      END
      subroutine DNQQRF(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
      integer M,N,LDA,LIPVT
      integer IPVT(LIPVT)
      logical PIVOT
      double precision A(LDA,N),RDIAG(N),ACNORM(N),WA(N)
C     **********
C
C     SUBROUTINE DNQQRF
C
C     THIS SUBROUTINE USES HOUSEHOLDER TRANSFORMATIONS WITH COLUMN
C     PIVOTING (OPTIONAL) TO COMPUTE A QR FACTORIZATION OF THE
C     M BY N MATRIX A. THAT IS, DNQQRF DETERMINES AN ORTHOGONAL
C     MATRIX Q, A PERMUTATION MATRIX P, and AN UPPER TRAPEZOIDAL
C     MATRIX R WITH DIAGONAL ELEMENTS OF NONINCREASING MAGNITUDE,
C     SUCH THAT A*P = Q*R. THE HOUSEHOLDER TRANSFORMATION FOR
C     COLUMN K, K = 1,2,...,MIN(M,N), is OF THE FORM
C
C                           T
C           I - (1/U(K))*U*U
C
C     WHERE U HAS ZEROS IN THE FIRST K-1 POSITIONS. THE FORM OF
C     THIS TRANSFORMATION and THE METHOD OF PIVOTING FIRST
C     APPEARED IN THE CORRESPONDING LINPACK SUBROUTINE.
C
C     THE SUBROUTINE STATEMENT IS
C
C       subroutine DNQQRF(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
C
C     WHERE
C
C       M is A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF ROWS OF A.
C
C       N is A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF COLUMNS OF A.
C
C       A is AN M BY N ARRAY. ON INPUT A CONTAINS THE MATRIX FOR
C         WHICH THE QR FACTORIZATION is TO BE COMPUTED. ON OUTPUT
C         THE STRICT UPPER TRAPEZOIDAL PART OF A CONTAINS THE STRICT
C         UPPER TRAPEZOIDAL PART OF R, and THE LOWER TRAPEZOIDAL
C         PART OF A CONTAINS A FACTORED FORM OF Q (THE NON-TRIVIAL
C         ELEMENTS OF THE U VECTORS DESCRIBED ABOVE).
C
C       LDA is A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.
C
C       PIVOT is A LOGICAL INPUT VARIABLE. IF PIVOT is SET TRUE,
C         THEN COLUMN PIVOTING is ENFORCED. IF PIVOT is SET FALSE,
C         THEN NO COLUMN PIVOTING is DONE.
C
C       IPVT is AN INTEGER OUTPUT ARRAY OF LENGTH LIPVT. IPVT
C         DEFINES THE PERMUTATION MATRIX P SUCH THAT A*P = Q*R.
C         COLUMN J OF P is COLUMN IPVT(J) OF THE IDENTITY MATRIX.
C         IF PIVOT is FALSE, IPVT is NOT REFERENCED.
C
C       LIPVT is A POSITIVE INTEGER INPUT VARIABLE. IF PIVOT is FALSE,
C         THEN LIPVT MAY BE AS SMALL AS 1. IF PIVOT is TRUE, THEN
C         LIPVT MUST BE AT LEAST N.
C
C       RDIAG is AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         DIAGONAL ELEMENTS OF R.
C
C       ACNORM is AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         NORMS OF THE CORRESPONDING COLUMNS OF THE INPUT MATRIX A.
C         IF THIS INFORMATION is NOT NEEDED, THEN ACNORM CAN COINCIDE
C         WITH RDIAG.
C
C       WA is A WORK ARRAY OF LENGTH N. IF PIVOT is FALSE, THEN WA
C         CAN COINCIDE WITH RDIAG.
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... D1MACH,DNRM2
C
C       FORTRAN-SUPPLIED ... max,sqrt,min
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
c     ------------------------------------------------------------------
      external D1MACH, DNRM2
      integer I,J,JP1,K,KMAX,MINMN
      double precision AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
      double precision D1MACH,DNRM2
      parameter(ONE = 1.0d0, P05 = 0.05d0, ZERO = 0.0d0)
C
C     EPSMCH is THE MACHINE PRECISION.
C
      EPSMCH = D1MACH(4)
C
C     COMPUTE THE INITIAL COLUMN NORMS and INITIALIZE SEVERAL ARRAYS.
C
      DO 10 J = 1, N
         ACNORM(J) = DNRM2(M,A(1,J),1)
         RDIAG(J) = ACNORM(J)
         WA(J) = RDIAG(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
C
C     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
C
      MINMN = min(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
C
C        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
C
         KMAX = J
         DO 20 K = J, N
            IF (RDIAG(K) .gt. RDIAG(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX .eq. J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         RDIAG(KMAX) = RDIAG(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
C
C        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
C        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
C
         AJNORM = DNRM2(M-J+1,A(J,J),1)
         IF (AJNORM .eq. ZERO) GO TO 100
         IF (A(J,J) .lt. ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
C
C        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
C        and UPDATE THE NORMS.
C
         JP1 = J + 1
         IF (N .lt. JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .or. RDIAG(K) .eq. ZERO) GO TO 80
            TEMP = A(J,K)/RDIAG(K)
            RDIAG(K) = RDIAG(K)*sqrt(max(ZERO,ONE-TEMP**2))
            IF (P05*(RDIAG(K)/WA(K))**2 .gt. EPSMCH) GO TO 80
            RDIAG(K) = DNRM2(M-J,A(JP1,K),1)
            WA(K) = RDIAG(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         RDIAG(J) = -AJNORM
  110    CONTINUE
      RETURN
C
C     Last line of subroutine DNQQRF.
C
      END
      subroutine DNQUPD(M,N,S,LS,U,V,W,SING)
      integer M,N,LS
      logical SING
      double precision S(LS),U(M),V(N),W(M)
C     **********
C
C     SUBROUTINE DNQUPD
C
C     GIVEN AN M BY N LOWER TRAPEZOIDAL MATRIX S, AN M-VECTOR U,
C     and AN N-VECTOR V, THE PROBLEM is TO DETERMINE AN
C     ORTHOGONAL MATRIX Q SUCH THAT
C
C                   T
C           (S + U*V )*Q
C
C     is AGAIN LOWER TRAPEZOIDAL.
C
C     THIS SUBROUTINE DETERMINES Q AS THE PRODUCT OF 2*(N - 1)
C     TRANSFORMATIONS
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     WHERE GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE
C     WHICH ELIMINATE ELEMENTS IN THE I-TH and N-TH PLANES,
C     RESPECTIVELY. Q ITSELF is NOT ACCUMULATED, RATHER THE
C     INFORMATION TO RECOVER THE GV, GW ROTATIONS is RETURNED.
C
C     THE SUBROUTINE STATEMENT IS
C
C       subroutine DNQUPD(M,N,S,LS,U,V,W,SING)
C
C     WHERE
C
C       M is A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF ROWS OF S.
C
C       N is A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF COLUMNS OF S. N MUST NOT EXCEED M.
C
C       S is AN ARRAY OF LENGTH LS. ON INPUT S MUST CONTAIN THE LOWER
C         TRAPEZOIDAL MATRIX S STORED BY COLUMNS. ON OUTPUT S CONTAINS
C         THE LOWER TRAPEZOIDAL MATRIX PRODUCED AS DESCRIBED ABOVE.
C
C       LS is A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
C         (N*(2*M-N+1))/2.
C
C       U is AN INPUT ARRAY OF LENGTH M WHICH MUST CONTAIN THE
C         VECTOR U.
C
C       V is AN ARRAY OF LENGTH N. ON INPUT V MUST CONTAIN THE VECTOR
C         V. ON OUTPUT V(I) CONTAINS THE INFORMATION NECESSARY TO
C         RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.
C
C       W is AN OUTPUT ARRAY OF LENGTH M. W(I) CONTAINS INFORMATION
C         NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED
C         ABOVE.
C
C       SING is A LOGICAL OUTPUT VARIABLE. SING is SET TRUE IF ANY
C         OF THE DIAGONAL ELEMENTS OF THE OUTPUT S ARE ZERO. OTHERWISE
C         SING is SET FALSE.
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... D1MACH
C
C       FORTRAN-SUPPLIED ... abs,sqrt
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE,
C     JOHN L. NAZARETH
C
C     **********
c     ------------------------------------------------------------------
      external D1MACH
      integer I,J,JJ,L,NMJ,NM1
      double precision VCOS,COTAN,GIANT,ONE,P5,P25,VSIN,VTAN,TAU,TEMP,
     *                 ZERO
      double precision D1MACH
      parameter(ONE = 1.0d0, P5 = 0.5d0, P25 = 0.25d0, ZERO = 0.0d0)
      save GIANT
      data GIANT / 0.0d0 /
C
C     GIANT is THE LARGEST MAGNITUDE.
C
      if(GIANT .eq. 0.0d0) GIANT = D1MACH(2)
C
C     INITIALIZE THE DIAGONAL ELEMENT POINTER.
C
      JJ = (N*(2*M - N + 1))/2 - (M - N)
C
C     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.
C
      L = JJ
      DO 10 I = N, M
         W(I) = S(L)
         L = L + 1
   10    CONTINUE
C
C     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
C     IN SUCH A WAY THAT A SPIKE is INTRODUCED INTO W.
C
      NM1 = N - 1
      IF (NM1 .lt. 1) GO TO 70
      DO 60 NMJ = 1, NM1
         J = N - NMJ
         JJ = JJ - (M - J + 1)
         W(J) = ZERO
         IF (V(J) .eq. ZERO) GO TO 50
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF V.
C
         IF (abs(V(N)) .ge. abs(V(J))) GO TO 20
            COTAN = V(N)/V(J)
            VSIN = P5/sqrt(P25+P25*COTAN**2)
            VCOS = VSIN*COTAN
            TAU = ONE
            IF (abs(VCOS)*GIANT .gt. ONE) TAU = ONE/VCOS
            GO TO 30
   20    CONTINUE
            VTAN = V(J)/V(N)
            VCOS = P5/sqrt(P25+P25*VTAN**2)
            VSIN = VCOS*VTAN
            TAU = VSIN
   30    CONTINUE
C
C        APPLY THE TRANSFORMATION TO V and STORE THE INFORMATION
C        NECESSARY TO RECOVER THE GIVENS ROTATION.
C
         V(N) = VSIN*V(J) + VCOS*V(N)
         V(J) = TAU
C
C        APPLY THE TRANSFORMATION TO S and EXTEND THE SPIKE IN W.
C
         L = JJ
         DO 40 I = J, M
            TEMP = VCOS*S(L) - VSIN*W(I)
            W(I) = VSIN*S(L) + VCOS*W(I)
            S(L) = TEMP
            L = L + 1
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.
C
      DO 80 I = 1, M
         W(I) = W(I) + V(N)*U(I)
   80    CONTINUE
C
C     ELIMINATE THE SPIKE.
C
      SING = .FALSE.
      IF (NM1 .lt. 1) GO TO 140
      DO 130 J = 1, NM1
         IF (W(J) .eq. ZERO) GO TO 120
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF THE SPIKE.
C
         IF (abs(S(JJ)) .ge. abs(W(J))) GO TO 90
            COTAN = S(JJ)/W(J)
            VSIN = P5/sqrt(P25+P25*COTAN**2)
            VCOS = VSIN*COTAN
            TAU = ONE
            IF (abs(VCOS)*GIANT .gt. ONE) TAU = ONE/VCOS
            GO TO 100
   90    CONTINUE
            VTAN = W(J)/S(JJ)
            VCOS = P5/sqrt(P25+P25*VTAN**2)
            VSIN = VCOS*VTAN
            TAU = VSIN
  100    CONTINUE
C
C        APPLY THE TRANSFORMATION TO S and REDUCE THE SPIKE IN W.
C
         L = JJ
         DO 110 I = J, M
            TEMP = VCOS*S(L) + VSIN*W(I)
            W(I) = -VSIN*S(L) + VCOS*W(I)
            S(L) = TEMP
            L = L + 1
  110       CONTINUE
C
C        STORE THE INFORMATION NECESSARY TO RECOVER THE
C        GIVENS ROTATION.
C
         W(J) = TAU
  120    CONTINUE
C
C        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
C
         IF (S(JJ) .eq. ZERO) SING = .TRUE.
         JJ = JJ + (M - J + 1)
  130    CONTINUE
  140 CONTINUE
C
C     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.
C
      L = JJ
      DO 150 I = N, M
         S(L) = W(I)
         L = L + 1
  150    CONTINUE
      IF (S(JJ) .eq. ZERO) SING = .TRUE.
      RETURN
C
C     Last line of subroutine DNQUPD.
C
      END
