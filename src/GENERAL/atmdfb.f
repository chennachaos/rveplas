CDOC BEGIN_SUBROUTINE ATMDFB
CDOC Computes the additional tangent modulus, q, for F-Bar elements
CDOC
CDOC This routine computes the additional tangent modulus, q, required
CDOC to assemble the tangent stiffness matrix for F-bar elements.
CDOC This routine contains the plane strain, axisymmetric and three-
CDOC dimensional implementations.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AMATX  >  Matrix of components of the current
CDOC C                          spatial tangent modulus, a.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION QMATX  <  The additional tangent modulus required
CDOC C                          to assemble the tangent stiffness matrix
CDOC C                          for F-bar elements.
CDOC DOUBLE_PRECISION STRES  >  Array of current Cauchy stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, September 1996: Initial coding
CHST   D. de Bortoli,     April 2015: Added MADIM as input argument;
CHST                                  Added 3-D case
CHST    
CHST      
      SUBROUTINE ATMDFB
     1(   AMATX      ,MADIM      ,NTYPE      ,QMATX      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION
     1    AMATX(MADIM,MADIM)    ,QMATX(MADIM,MADIM)   ,STRES(*)
      DATA
     1    R0   ,RP5  ,R1   ,R2   ,R3   /
     2    0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
C***********************************************************************
C COMPUTE THE ADDITIONAL TANGENT MODULUS "q" REQUIRED BY F-BAR
C ELEMENTS:
C                       1                 2
C               q  :=  --- a:(I (x) I) - --- [sigma] (x) I
C                       3                 3
C FOR THE AXISYMMETRIC AND THREE-DIMENSIONAL CASES, AND
C                       1                 1
C               q  :=  --- a:(I (x) I) - --- [sigma] (x) I
C                       2                 2
C FOR PLANE STRAIN.
C
C REFERENCE: Expressions (15.11) and (15.22)
C***********************************************************************
C
      IF(NTYPE.EQ.2)THEN
C Plane strain
        A=RP5
        B=-RP5
        QMATX(1,1)=A*(AMATX(1,1)+AMATX(1,4))+B*STRES(1)
        QMATX(2,1)=A*(AMATX(2,1)+AMATX(2,4))+B*STRES(3)
        QMATX(3,1)=A*(AMATX(3,1)+AMATX(3,4))+B*STRES(3)
        QMATX(4,1)=A*(AMATX(4,1)+AMATX(4,4))+B*STRES(2)
        QMATX(1,2)=R0
        QMATX(2,2)=R0
        QMATX(3,2)=R0
        QMATX(4,2)=R0
        QMATX(1,3)=R0
        QMATX(2,3)=R0
        QMATX(3,3)=R0
        QMATX(4,3)=R0
        QMATX(1,4)=QMATX(1,1)
        QMATX(2,4)=QMATX(2,1)
        QMATX(3,4)=QMATX(3,1)
        QMATX(4,4)=QMATX(4,1)
      ELSEIF(NTYPE.EQ.3)THEN
C Axisymmetric
        A=R1/R3
        B=-R2/R3
        QMATX(1,1)=A*(AMATX(1,1)+AMATX(1,4)+AMATX(1,5))+B*STRES(1)
        QMATX(2,1)=A*(AMATX(2,1)+AMATX(2,4)+AMATX(2,5))+B*STRES(3)
        QMATX(3,1)=A*(AMATX(3,1)+AMATX(3,4)+AMATX(3,5))+B*STRES(3)
        QMATX(4,1)=A*(AMATX(4,1)+AMATX(4,4)+AMATX(4,5))+B*STRES(2)
        QMATX(5,1)=A*(AMATX(5,1)+AMATX(5,4)+AMATX(5,5))+B*STRES(4)
        QMATX(1,2)=R0
        QMATX(2,2)=R0
        QMATX(3,2)=R0
        QMATX(4,2)=R0
        QMATX(5,2)=R0
        QMATX(1,3)=R0
        QMATX(2,3)=R0
        QMATX(3,3)=R0
        QMATX(4,3)=R0
        QMATX(5,3)=R0
        QMATX(1,4)=QMATX(1,1)
        QMATX(2,4)=QMATX(2,1)
        QMATX(3,4)=QMATX(3,1)
        QMATX(4,4)=QMATX(4,1)
        QMATX(5,4)=QMATX(5,1)
        QMATX(1,5)=QMATX(1,1)
        QMATX(2,5)=QMATX(2,1)
        QMATX(3,5)=QMATX(3,1)
        QMATX(4,5)=QMATX(4,1)
        QMATX(5,5)=QMATX(5,1)
      ELSEIF(NTYPE.EQ.4)THEN
C Three-dimensional
        A=R1/R3
        B=-R2/R3
        QMATX(1,1)=A*(AMATX(1,1)+AMATX(1,5)+AMATX(1,9))+B*STRES(1)
        QMATX(2,1)=A*(AMATX(2,1)+AMATX(2,5)+AMATX(2,9))+B*STRES(4)
        QMATX(3,1)=A*(AMATX(3,1)+AMATX(3,5)+AMATX(3,9))+B*STRES(6)
        QMATX(4,1)=A*(AMATX(4,1)+AMATX(4,5)+AMATX(4,9))+B*STRES(4)
        QMATX(5,1)=A*(AMATX(5,1)+AMATX(5,5)+AMATX(5,9))+B*STRES(2)
        QMATX(6,1)=A*(AMATX(6,1)+AMATX(6,5)+AMATX(6,9))+B*STRES(5)
        QMATX(7,1)=A*(AMATX(7,1)+AMATX(7,5)+AMATX(7,9))+B*STRES(6)
        QMATX(8,1)=A*(AMATX(8,1)+AMATX(8,5)+AMATX(8,9))+B*STRES(5)
        QMATX(9,1)=A*(AMATX(9,1)+AMATX(9,5)+AMATX(9,9))+B*STRES(3)
C
        QMATX(:,2:4)=R0
C
        QMATX(:,5)=QMATX(:,1)
C
        QMATX(:,6:8)=R0
C
        QMATX(:,9)=QMATX(:,1)
      ELSE
C Invalid analysis type
        CALL ERRPRT('EI0073')
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE ATMDFB
