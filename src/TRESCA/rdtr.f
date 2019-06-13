CDOC BEGIN_SUBROUTINE RDTR
CDOC Reads data for the Tresca elasto-plastic material model.
CDOC
CDOC This routine reads from the data file and echos to the results file
CDOC the material parameters necessary for the Tresca elasto-plastic
CDOC model with piece-wise linear isotropic hardening.
CDOC It also sets the array of real properties and some components of
CDOC the array of integer material properties.
CDOC These arrays are used by subroutines SUTR and CTTR.
CDOC This routine also checks whether the number of integer and real
CDOC material properties, the number of real state variables and the
CDOC number of logical algorithmic variables required by the present
CDOC model is compatible with the dimension of the corresponding global
CDOC arrays.
CDOC It also sets the unsymmetric tangent stiffness flag.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IPROPS <  Array of integer material properties.
CDOC INTEGER          MIPROP >  Dimension of the global array of integer
CDOC C                          material properties.
CDOC INTEGER          MLALGV >  Dimension of the global array of logical
CDOC C                          algorithmic variables.
CDOC INTEGER          MRPROP >  Dimension of the global array of real
CDOC C                          material variables.
CDOC INTEGER          MRSTAV >  Dimension of the global array of real
CDOC C                          state variables.
CDOC INTEGER          NLARGE >  Large deformation analysis flag.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RPROPS <  Array of real material properties.
CDOC LOGICAL          UNSYM  <  Unsymmetric tangent stiffness flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June 1996: Initial coding
CHST
CHST E.de Souza Neto, April 2011: I/O error message added
CHST
      SUBROUTINE RDTR
     1(   IPROPS     ,MIPROP     ,MLALGV     ,MRALGV     ,MRPROP     ,
     2    MRSTAV     ,NLARGE     ,NTYPE      ,RPROPS     ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UNSYM
      PARAMETER( IPHARD=4  ,NIPROP=3  ,NLALGV=4  ,NRALGV=2  ,NRSTAV=5 )
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)
      DATA R0 / 0.0D0 /
C***********************************************************************
C READ AND ECHO MATERIAL PROPERTIES FOR TRESCA TYPE ELASTO-PLASTIC
C MATERIAL WITH NON-LINEAR ISOTROPIC HARDENING
C
C REFERENCE: Section 8.1
C***********************************************************************
 1000 FORMAT(' Elasto-plastic with TRESCA yield criterion'/)
 1100 FORMAT(
     1' Mass density ...................................... =',G15.6/
     2' Young''s modulus ................................... =',G15.6/
     3' Poisson''s ratio ................................... =',G15.6)
 1200 FORMAT(/
     1' Number of points on hardening curve ............... =',I3//
     2'           Epstn        uniaxial yield stress '/)
 1300 FORMAT(2(5X,G15.6))
C
C Model not yet implemented for large strains with plane stress
      IF(NLARGE.EQ.1.AND.NTYPE.EQ.1)CALL ERRPRT('ED0198')
C Set unsymmetric tangent stiffness flag
      UNSYM=.FALSE.
C Read and echo some of the real properties
      WRITE(16,1000)
      READ(15,*,ERR=100,END=100)DENSE
      READ(15,*,ERR=100,END=100)YOUNG,POISS
      WRITE(16,1100)DENSE,YOUNG,POISS
      IF(YOUNG.LT.R0)CALL ERRPRT('ED0107')
C Hardening curve
      READ(15,*,ERR=100,END=100)NHARD
      WRITE(16,1200)NHARD
      IF(NHARD.LT.2) CALL ERRPRT('ED0108')
C Check dimensions of IPROPS
      IF(MIPROP.LT.NIPROP)CALL ERRPRT('ED0109')
      IPROPS(3)=NHARD
C Check dimensions of RPROPS
      NRPROP=IPHARD+NHARD*2-1
      IF(NRPROP.GT.MRPROP)CALL ERRPRT('ED0110')
C Store real properties in RPROPS
      RPROPS(1)=DENSE
      RPROPS(2)=YOUNG
      RPROPS(3)=POISS
C
C Read and set hardening curve
C
      DO 10 IHARD=1,NHARD
        READ(15,*,ERR=100,END=100)RPROPS(IPHARD+IHARD*2-2),
     1                            RPROPS(IPHARD+IHARD*2-1)
        WRITE(16,1300)RPROPS(IPHARD+IHARD*2-2),
     1                RPROPS(IPHARD+IHARD*2-1)
   10 CONTINUE
C
C Check dimension of RSTAVA, LALGVA and RALGVA
C
      IF(NRSTAV.GT.MRSTAV)CALL ERRPRT('ED0111')
      IF(NLALGV.GT.MLALGV)CALL ERRPRT('ED0112')
      IF(NRALGV.GT.MRALGV)CALL ERRPRT('ED0113')
C
      GOTO 200
C Issue error message and abort program execution in case of I/O error
  100 CALL ERRPRT('ED0210')
C
  200 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RDTR
