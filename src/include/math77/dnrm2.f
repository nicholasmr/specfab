      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 1998-05-11 DNRM2  Krogh   Minor changes for conversion to C.
c>> 1996-08-29 DNRM2  Krogh   Coded an entirely different algorithm
c>> ....
C>> 1985-08-02 DNRM2  Lawson  Initial code.
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN X() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C New algorithm avoids underflow as well as overflow, and avoids divides
c as much as possible.  F. Krogh, August 29, 1996.
c
c ************************* Variable Definitions ***********************
c
c Let b be the base for floating point, b**E be the smallest integer
c power of b that overflows, b**e be the smallest integer power of b
c that does not underflow, and d be the number of base b digits in a
c a floating point number.  The descriptions below use b, e, E, and d.
c
c ABIG = b ** ((E+5) / 4)).    This is used as a factor to get the
c   final answer when combining SUMX from big XA, with SUM.
c ASMALL = b ** ((e+d-E+5) / 4)).    This is used as a factor to get the
c   final answer when combining SUMX from small XA, with SUM.
c DBIG = b ** (2*((E+5)/4)).  This is used as a factor to get the
c   final answer when only SUMX formed from big XA's is needed.
c DSMALL = b ** (2*((e+d-E+5)/4)).  This is used as a factor to get the
c   final answer when only SUMX formed from big XA's is needed.
c FBIG = b ** (-2*((E+5)/4)).  This is used as a multiplier when
c   accumulating SUMX for big XA.
c FSMALL = b ** (-2*((e+d-E+5)/4)).  This is used as a multiplier when
c   accumulating SUMX for small XA.
c I      Temporary index.
c ID     Number of base I1MACH(10) digits in Floating point number.
c IEMN   The minimum floating point exponent.
c IEMX   The maximum floating point exponent.
c INCX   Input, the increment (>0) between elements of X.
c N      Input, the number of elements in X.
c NN     N * INCX = last index processed.
c SUM    Place where sum of XA**2 is accumlated.
c SUMX   Place where scaled sum of XA**2 is accumlated.  In first loop
c   this is for XA**2 that would likely underflow, and in second loop
c   this is for XA**2 that would likely overflow.
c TBIG = b ** ((E - d - 1) / 2).  If XA > TBIG, number is "big".  Note
c   that for such XA's, (FBIG * XA) ** 2 should not underflow and the
c   accumlation overflows only if the final result would.
c TSMALL = b ** ((e+1) / 2).  If XA <= TSMALL, number is "small".  Note
c   that for such XA's, (FSMALL * XA)**2 should not underflow and the
c   accumlation should not overlflow.
c X      Input, we are getting the L2 norm of X.
c XA     Contains base of floating point numbers when getting saved
c   parameters.  Later contains abs(X(I)).
C     ------------------------------------------------------------------
c--D replaces "?": ?NRM2
C     ------------------------------------------------------------------
      integer N, INCX
      double precision X(*)
      external I1MACH
      integer I, ID, IEMN, IEMX, I1MACH, NN
      double precision SUM, SUMX, XA
      double precision ABIG,ASMALL,DBIG,DSMALL,FBIG,FSMALL,TBIG,TSMALL
      save ABIG,ASMALL,DBIG,DSMALL,FBIG,FSMALL,TBIG,TSMALL
c Values of these saved parameters for D.P. IEEE arithmetic are:
c   ABIG= .2315841784746324E+78     DBIG= .5363123171977043E+155
c   FBIG= .1864585182800050E-154    TBIG= .9989595361011182E+146
c ASMALL= .4887898181599363E-149  DSMALL= .2389154863368240E-298
c FSMALL= .4185580496821357E+299  TSMALL= .2983336292480080E-153
c Values of these saved parameters for S.P. IEEE arithmetic are:
c   ABIG= .8589935E+10      DBIG= .7378698E+20
c   FBIG= .1355253E-19      TBIG= .2251800E+16
c ASMALL= .1387779E-16    DSMALL= .1925930E-33
c FSMALL= .5192297E+34    TSMALL= .2168404E-18
c
      data ABIG / 0.D0 /
C     ------------------------------------------------------------------
c
c
      if (ABIG .eq. 0.D0) then
C++ Code for (.N. == 'D') is active
         IEMX = I1MACH(16)
         IEMN = I1MACH(15)
         ID = I1MACH(14)
C++ Code for (.N. == 'S') is inactive
C         IEMX = I1MACH(13)
C         IEMN = I1MACH(12)
C         ID = I1MACH(11)
C++ END
         XA = dble(I1MACH(10))
         ABIG = XA ** ((IEMX+5)/4)
         DBIG = ABIG ** 2
         FBIG = 1.D0 / DBIG
         TBIG = XA ** ((IEMX - ID - 1) / 2)
         ASMALL = XA ** ((IEMN + ID - IEMX + 5) / 4)
         DSMALL = ASMALL ** 2
         FSMALL = 1.D0 / DSMALL
         TSMALL = XA ** ((IEMN + 1) / 2)
      end if
      SUM = 0.D0
      if (N .gt. 0) then
         NN = N * INCX
         SUMX = 0.D0
c                      Loop when no big number yet encountered.
         do 100 I = 1, NN, INCX
            XA = abs(X(I))
            if (XA .lt. TSMALL) then
               SUMX = SUMX + (FSMALL * XA) ** 2
            else
               if (XA .gt. TBIG) go to 200
               SUM = SUM + XA**2
            end if
  100    continue
         if (SUM .ne. 0.D0) then
            if (SUMX .ge. 1.D0) then
               if (SUM .lt. 1.D0) then
                  SUM = ASMALL * sqrt(FSMALL*SUM + DSMALL*SUMX)
                  go to 400
               end if
            end if
            SUM = sqrt(SUM)
         else
            SUM = DSMALL * sqrt(SUMX)
         end if
         go to 400
c
  200    SUMX = 0.D0
c                      Loop when we have at least one big number.
         do 300 I = I, NN, INCX
            XA = abs(X(I))
            if (XA .gt. TSMALL) then
               if (XA .gt. TBIG) then
                  SUMX = SUMX + (FBIG * XA) ** 2
               else
                  SUM = SUM + XA**2
               end if
            end if
  300    continue
         if ((SUMX .le. 1.D10) .and. (SUM .ge. 1.D-10)) then
            SUM = ABIG * sqrt(FBIG*SUM + DBIG*SUMX)
         else
            SUM = DBIG * sqrt(SUMX)
         end if
      end if
  400 continue
      DNRM2 = SUM
      return
      end
