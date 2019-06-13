CDOC BEGIN_SUBROUTINE JACDET
CDOC Jacobian determinant and cartesian derivatives for 2-D and 3-D 
CDOC isoparametric elements.
CDOC
CDOC This routine computes the Jacobian of the isoparametric mapping for
CDOC 2-D and 3-D isoparametric elements and returns the Jacobian 
CDOC determinant and the cartesian derivatives of the shape functions.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION CARTD  <  Matrix of cartesian derivatives of the
CDOC C                          element shape functions.
CDOC DOUBLE_PRECISION DERIV  >  Matrix of isoparametric derivatives of
CDOC C                          the element shape functions.
CDOC DOUBLE_PRECISION DETJAC <  Isoparametric mapping Jacobian
CDOC C                          determinant.
CDOC DOUBLE_PRECISION ELCOD  >  Array of element nodal coordinates.
CDOC INTEGER          IELEM  >  Element number (used only if warning
CDOC C                          messages are issued).
CDOC INTEGER          MDIME  >  Dimensioning parameter for arguments
CDOC C                          CARTD, DERIV and ELCOD.
CDOC INTEGER          NDIME  >  Number of spatial dimensions.
CDOC INTEGER          NNODE  >  Number of nodes of the element.
CDOC END_PARAMETERS
CDOC
CHST
CHST D. de Bortoli, March 2015: Added 3-D case, changing name from 
CHST                            'JACOB2' to 'JACDET'
CHST
      SUBROUTINE JACDET
     1(   CARTD      ,DERIV      ,DETJAC     ,ELCOD      ,IELEM      ,
     2    MDIME      ,NDIME      ,NNODE      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    CARTD(MDIME,*)     ,DERIV(MDIME,*)     ,ELCOD(MDIME,*)
      DIMENSION
     1    XJACI(NDIME,NDIME)  ,XJACM(NDIME,NDIME)
      DATA R0   /
     1     0.0D0/
C***********************************************************************
C EVALUATES THE JACOBIAN MATRIX, ITS DETERMINANT AND THE CARTESIAN
C DERIVATIVES OF THE SHAPE FUNCTIONS OF 2-D AND 3-D ISOPARAMETRIC 
C ELEMENTS
C
C REFERENCE: Section 4.1.2
C            Expression (4.33)
C***********************************************************************
 1000 FORMAT(//'  Warning from subroutine JACDET:'//
     1 10X,'Negative jacobian determinant ',G12.4,' Element number ',I5)
 1010 FORMAT(//'  Warning from subroutine JACDET:'//
     1 10X,'Zero jacobian determinant ',12X,' Element number ',I5/
     2 10X,'Jacobian matrix not inverted and cartesian derivatives ',/,
     3 10X,'of shape functions not computed')
C
C Check if number of spatial dimensions is valid
      IF((NDIME.NE.2).AND.(NDIME.NE.3))CALL ERRPRT('EI0078')
C
C Evaluate jacobian matrix XJACM
C ------------------------------
      DO 30 IDIME=1,NDIME
        DO 20 JDIME=1,NDIME
          XJACM(IDIME,JDIME)=R0
          DO 10 INODE=1,NNODE
            XJACM(IDIME,JDIME)=XJACM(IDIME,JDIME)+DERIV(IDIME,INODE)*
     1                         ELCOD(JDIME,INODE)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C Determinant of jacobian matrix
      IF(NDIME.EQ.2)THEN
        DETJAC=XJACM(1,1)*XJACM(2,2)-XJACM(1,2)*XJACM(2,1)
      ELSEIF(NDIME.EQ.3)THEN
C In 3-D, INVMT3 calculates both the inverse and the determinant of the
C jacobian matrix
        CALL INVMT3(XJACM,XJACI,DETJAC)
      ENDIF
      IF(DETJAC.LT.R0)THEN
        WRITE(*,1000)DETJAC,IELEM
        WRITE(16,1000)DETJAC,IELEM
      ELSEIF(DETJAC.EQ.R0)THEN
        WRITE(*,1010)IELEM
        WRITE(16,1010)IELEM
        GOTO 999
      ENDIF
C Inverse of jacobian matrix (2-D case)
      IF(NDIME.EQ.2)THEN
        XJACI(1,1)=XJACM(2,2)/DETJAC
        XJACI(2,2)=XJACM(1,1)/DETJAC
        XJACI(1,2)=-XJACM(1,2)/DETJAC
        XJACI(2,1)=-XJACM(2,1)/DETJAC
      ENDIF
C
C Evaluate cartesian derivatives of shape functions
C -------------------------------------------------
C
      DO 60 IDIME=1,NDIME
        DO 50 INODE=1,NNODE
          CARTD(IDIME,INODE)=R0
          DO 40 JDIME=1,NDIME
            CARTD(IDIME,INODE)=CARTD(IDIME,INODE)+XJACI(IDIME,JDIME)*
     1                         DERIV(JDIME,INODE)
   40     CONTINUE
   50   CONTINUE
   60 CONTINUE
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE JACDET
