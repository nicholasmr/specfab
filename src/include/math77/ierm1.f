      subroutine IERM1(SUBNAM,INDIC,LEVEL,MSG,LABEL,VALUE,FLAG)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
C>> 1990-01-18 CLL Added Integer stmt for VALUE.  Typed all variables.
C>> 1985-08-02 IERM1  Lawson  Initial code.
C
      integer INDIC, LEVEL, VALUE
      character*(*) SUBNAM,MSG,LABEL
      character*1 FLAG
      call ERMSG(SUBNAM,INDIC,LEVEL,MSG,',')
      call IERV1(LABEL,VALUE,FLAG)
C
      return
      end
