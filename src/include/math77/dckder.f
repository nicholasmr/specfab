      subroutine DCKDER(MODE, M, N, X, FVEC, FJAC, LDFJAC,
     *                  TEST, IMAX, JMAX, TSTMAX)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 1996-03-30 DCKDER Krogh  Added external statement.
c>> 1994-10-20 DCKDER Krogh  Changes to use M77CON
c>> 1992-01-16 DCKDER CLL@JPL
c>> 1991-12-16 CLL@JPL
C  DCKDER checks a Jacobian matrix of first partial derivatives for
c  consistency with central finite differences of function values.
c
c     Usage:
c
c        Set X() to a trial N-vector.
c        Compute FJAC(,) as the MxN Jacobian matrix of partials of
c           FVEC with respect to X, evaluated at X().
c        MODE = 1
c     10 call DCKDER(MODE,...)
c        if(MODE .eq. 2) then
c           Compute FVEC() as an M-vector of function values
c              evaluated at X().
c           go to 10
c        endif
c        Here the process is completed.
c        Results are in IMAX, JMAX, TSTMAX, and TEST(,).
c     ------------------------------------------------------------------
C                         Subroutine Arguments
c
C  MODE [inout]  On initial entry set MODE = 1.
c       DCKDER will return a number of times with MODE = 2.  The calling
c       code should compute FVEC() as a function of X() and call DCKDER
c       again, not altering MODE.  When DCKDER is finished it returns
c       with MODE = 3.
c  M [in]  Number of terms in FVEC() and number of rows of data
c          in FJAC(,).
c  N [in]  Number of terms in X() and number of cloumns of data
c          in FJAC(,).
c  X() [inout]  Initially must contain the x-vector around which the
c          testing will be done.  Contains perturbed x-vectors on
c          returns with MODE = 2.  On final return with MODE = 3,
c          contains the original x-vector exactly restored.
c  FVEC() [in]  The user stores function values in FVEC().
c  FJAC(,) [in]  The user stores the Jacobian matrix in FVEC(,).
c  LDFJAC [in]  Declared first dimension of the arrays FJAC(,) and
c               TEST(,).  Require LDFJAC .ge. M.
c  TEST(,) [out]  Array with same dimensions as FJAC().
c     For each i = 1, ..., M, and j = 1, ..., N, DCKDER sets TEST(i,j)
c     to the signed difference: FJAC(i,j) minus a central finite-
c     difference approximation to the first partial derivative of
c     f sub i with respect to x sub j.
c  IMAX, JMAX [out]  I and J indices of the value in TEST(,) of
c               largest magnitude.
c  TSTMAX [out]  Magnitude of the value in TEST(,) of largest magnitude,
c               i.e., TSTMAX = abs(TEST(IMAX, JMAX))
c     ------------------------------------------------------------------
c  D1MACH(1) is the underflow limit.
c  D1MACH(4) is the machine precision.
c     ------------------------------------------------------------------
c--D replaces "?": ?CKDER
c     ------------------------------------------------------------------
      external D1MACH
      double precision D1MACH
      integer I, IMAX, J, JMAX, LDFJAC, MODE, M, N
      double precision ALPHA, DELX, EPS, FAC, FJAC(LDFJAC,N), FVEC(M)
      double precision SMALL, TSTIJ, TEST(LDFJAC,N), TSTMAX, X(N), XJ
      logical PLUS
      save ALPHA, DELX, J, PLUS, SMALL, XJ
c     ------------------------------------------------------------------
      go to (10, 20), MODE
c                               ANSI Standard Fortran 77 comes here
c                               if MODE is not 1 or 2.
c
      call IERM1('DCKDER',1,0,'Require MODE = 1 or 2.','MODE',MODE,'.')
      return
c
   10 continue
      EPS = D1MACH(4)
      ALPHA = (3.0d0 * EPS)**0.333333d0
      SMALL = max(1.0d5 * D1MACH(1) / ALPHA, EPS*EPS)
      TSTMAX = 0.0d0
      IMAX = 0
      JMAX = 0
      MODE = 2
      J = 0
   15 continue
         J = J+1
         XJ = X(J)
         if(abs(XJ) .gt. SMALL) then
            DELX = ALPHA * XJ
         elseif(abs(XJ) .gt. 0.0d0) then
            DELX = ALPHA * SMALL
         else
            DELX = ALPHA
         endif
         X(J) = XJ + DELX
         PLUS = .true.
         return
c
c           Here we return to the user's code for computation of
c           FVEC() using X().  Execution will resume here.
c
   20    continue
         if(PLUS) then
c                         Using col N of TEST() to save FVEC()
c                         evaluated using XJ + DELX.
            do 30 I = 1,M
               TEST(I,N) = FVEC(I)
   30       continue
            X(J) = XJ - DELX
            PLUS = .false.
            return
c                  Returning again to the user's code for another
c                  evaluation of FVEC() using the perturbed X().
         endif
c
         FAC = 0.5d0 / DELX
         do 50 I = 1,M
            TSTIJ = FJAC(I,J) - FAC*(TEST(I,N) - FVEC(I))
            TEST(I,J) = TSTIJ
            if(abs(TSTIJ) .gt. TSTMAX) then
               TSTMAX = abs(TSTIJ)
               IMAX = I
               JMAX = J
            endif
   50    continue
c                                     Restore X(J)
         X(J) = XJ
      if(J .lt. N) go to 15
      MODE = 3
      return
      end
