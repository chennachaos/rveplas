CDOC BEGIN_SUBROUTINE TANGEN
CDOC Sets prescribed displacements for arc-length tangential solution.
CDOC
CDOC This routine sets the array of prescribed displacements as needed
CDOC for the tangential solution of the Arc-Length Method.
CHST
CHST M.E.Honnor, April 1987
CHST
      SUBROUTINE TANGEN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
C
      DATA R0/0.0D0/
C***********************************************************************
C SETS UP PRESCRIBED DISPLACEMENTS FOR THE TANGENTIAL SOLUTION FOR
C THE ARC-LENGTH METHOD
C
C REFERENCE: Item (iii), Box 4.4
C***********************************************************************
      DO 20 IVFIX=1,NVFIX
        NLOCA=(NOFIX(IVFIX)-1)*NDOFN
        DO 10 IDOFN=1,NDOFN
          NPOS=NLOCA+IDOFN
          FIXED(NPOS,1)=PRESC(IVFIX,IDOFN)
          FIXED(NPOS,2)=R0
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE TANGEN
