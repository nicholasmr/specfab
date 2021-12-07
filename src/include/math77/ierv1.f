      SUBROUTINE IERV1(LABEL,VALUE,FLAG)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 1995-11-15 IERV1 Krogh  Moved format up for C conversion.
C>> 1985-09-20 IERV1  Lawson  Initial code.
C
C     ------------------------------------------------------------
C     SUBROUTINE ARGUMENTS
C     --------------------
C     LABEL     An identifing name to be printed with VALUE.
C
C     VALUE     A integer to be printed.
C
C     FLAG      See write up for FLAG in ERMSG.
C
C     ------------------------------------------------------------
C
      COMMON/M77ERR/IDELTA,IALPHA
      INTEGER IDELTA,IALPHA,VALUE
      CHARACTER*(*) LABEL
      CHARACTER*1 FLAG
      SAVE /M77ERR/
 1002 FORMAT(3X,A,' = ',I5)
C
      IF (IALPHA.GE.-1) THEN
        WRITE (*,1002) LABEL,VALUE
        IF (FLAG .EQ. '.') CALL ERFIN
      ENDIF
      RETURN
C
      END
