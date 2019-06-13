CDOC BEGIN_SUBROUTINE RDPDSC
CDOC Read data for planar double-slip single crystal model.
CDOC
CDOC This routine reads from the data file and echos to the results file
CDOC the material parameters necessary for the large strain
CDOC planar double-slip single crystal elasto-plastic model with
CDOC piece-wise linear isotropic Taylor hardening.
CDOC It also sets the array of real properties and some components of
CDOC the array of integer material properties.
CDOC These arrays are used by subroutines SUPDSC and CSTPDS.
CDOC This routine also checks whether the number of integer and real
CDOC material properties, the number of real state variables and the
CDOC number of logical and real algorithmic variables required by the
CDOC present model is compatible with the dimension of the corresponding
CDOC global arrays.
CDOC It also sets the unsymmetric tangent stiffness flag.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IPROPS <  Array of integer material properties.
CDOC INTEGER          MIPROP >  Dimension of the global array of integer
CDOC C                          material properties.
CDOC INTEGER          MLALGV >  Dimension of the global array of logical
CDOC C                          algorithmic variables.
CDOC INTEGER          MRALGV >  Dimension of the global array of real
CDOC C                          algorithmic variables.
CDOC INTEGER          MRPROP >  Dimension of the global array of real
CDOC C                          material variables.
CDOC INTEGER          MRSTAV >  Dimension of the global array of real
CDOC C                          state variables.
CDOC INTEGER          NLARGE >  Large strain analysis integer flag.
CDOC C                          Used for checking only.
CDOC C                          Present model is compatible only with
CDOC C                          large strain analysis.
CDOC INTEGER          NTYPE  >  Stress state type integer flag.
CDOC C                          Used for checking only.
CDOC C                          Present model is compatible only with
CDOC C                          plane strain analysis (NTYPE=2).
CDOC DOUBLE_PRECISION RPROPS <  Array of real material properties.
CDOC LOGICAL          UNSYM  <  Unsymmetric tangent stiffness flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, November 1998: Initial coding
CHST
CHST E.de Souza Neto, April 2011: I/O error message added
CHST
      SUBROUTINE RDPDSC
     1(   IPROPS     ,MIPROP     ,MLALGV     ,MRALGV     ,MRPROP     ,
     2    MRSTAV     ,NLARGE     ,NTYPE      ,RPROPS     ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   IPHARD=6   ,NIPROP=3   ,NLALGV=6   ,NRALGV=4   ,NRSTAV=5   )
C Arguments
      LOGICAL UNSYM
      DIMENSION
     1    IPROPS(MIPROP)     ,RPROPS(MRPROP)
C Local data
      DATA R1   ,R180   /
     1     1.0D0,180.0D0/
C***********************************************************************
C READ AND ECHO MATERIAL PROPERTIES FOR THE LARGE STRAIN PLANAR
C DOUBLE-SLIP SINGLE CRYSTAL ELASTO-PLASTIC MODEL WITH NON-LINEAR
C (PIECEWISE LINEAR) ISOTROPIC TAYLOR HARDENING
C***********************************************************************
 1000 FORMAT(' Large strain planar double-slip SINGLE CRYSTAL'/)
 1100 FORMAT(
     1' Mass density ...................................... =',G15.6/
     2' Shear modulus ..................................... =',G15.6/
     3' Bulk modulus ...................................... =',G15.6/
     4' Initial orientation of FIRST SLIP SYSTEM relative'          /
     5' to X-axis (degrees, counterclockwise-positive) .... =',G15.6/
     6' Initial orientation of SECOND SLIP SYSTEM relative'         /
     7' to first syst. (degrees, counterclockwise-positive) =',G15.6)
 1200 FORMAT(/
     1' Number of points on Taylor hardening curve ........ =',I3//
     2'       Accum. slip    Resolved Schmid yield stress '/)
 1300 FORMAT(2(5X,G15.6))
C
C Set unsymmetric tangent stiffness flag
      UNSYM=.TRUE.
C Check that stress state type and large strain flags are compatible
C with the present model
      IF(NTYPE.NE.2)CALL ERRPRT('ED0154')
      IF(NLARGE.NE.1)CALL ERRPRT('ED0155')
C
C Read and echo some of the real properties
      WRITE(16,1000)
      READ(15,*,ERR=100,END=100)DENSE
      READ(15,*,ERR=100,END=100)GMODU,BULK,THETA,BETA
      WRITE(16,1100)DENSE,GMODU,BULK,THETA,BETA
C number of points on hardening curve
      READ(15,*,ERR=100,END=100)NHARD
      WRITE(16,1200)NHARD
      IF(NHARD.LT.2) CALL ERRPRT('ED0148')
C check dimensions of IPROPS
      IF(NIPROP.GT.MIPROP)CALL ERRPRT('ED0149')
      IPROPS(3)=NHARD
C check dimensions of RPROPS
      NRPROP=IPHARD+NHARD*2-1
      IF(NRPROP.GT.MRPROP)CALL ERRPRT('ED0150')
C convert angles into radians
      RADEG=ACOS(-R1)/R180
      THETA=RADEG*THETA
      BETA=RADEG*BETA
      RPROPS(1)=DENSE
      RPROPS(2)=GMODU
      RPROPS(3)=BULK
      RPROPS(4)=THETA
      RPROPS(5)=BETA
C Read and set hardening curve
      DO 10 IHARD=1,NHARD
        READ(15,*,ERR=100,END=100)RPROPS(IPHARD+IHARD*2-2),
     1                            RPROPS(IPHARD+IHARD*2-1)
        WRITE(16,1300)RPROPS(IPHARD+IHARD*2-2),
     1                RPROPS(IPHARD+IHARD*2-1)
   10 CONTINUE
C Check dimension of RSTAVA, RALGVA and LALGVA
      IF(NRSTAV.GT.MRSTAV)CALL ERRPRT('ED0151')
      IF(NRALGV.GT.MRALGV)CALL ERRPRT('ED0152')
      IF(NLALGV.GT.MLALGV)CALL ERRPRT('ED0153')
C
      GOTO 200
C Issue error message and abort program execution in case of I/O error
  100 CALL ERRPRT('ED0200')
C
  200 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RPDSC
