CDOC BEGIN_SUBROUTINE RDVMMX
CDOC Read data for the von Mises plasticity model with mixed hardening.
CDOC
CDOC This routine reads from the data file and echos to the results file
CDOC the material parameters necessary for
CDOC the von Mises elasto-plastic model with piece-wise linear
CDOC mixed isotropic/kinematic hardening.
CDOC It also sets the array of real properties and some components of
CDOC the array of integer material properties.
CDOC These arrays are used by subroutines SUVMMX and CTVMMX.
CDOC This routine also checks whether the number of integer and real
CDOC material properties, the number of real state variables and the
CDOC number of logical algorithmic variables required by the present
CDOC model
CDOC is compatible with the dimension of the corresponding global
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
CDOC INTEGER          NLARGE >  Large strain analysis flag.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RPROPS <  Array of real material properties.
CDOC LOGICAL          UNSYM  <  Unsymmetric tangent stiffness flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1999: Initial coding
CHST
CHST E.de Souza Neto, April 2011: I/O error message added
CHST
      SUBROUTINE RDVMMX
     1(   IPROPS     ,MIPROP     ,MLALGV     ,MRPROP     ,MRSTAV     ,
     2    NLARGE     ,NTYPE      ,RPROPS     ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UNSYM
      PARAMETER( IPHARD=4  ,NLALGV=2  ,NRSTAV=9 )
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)
      DATA R0   /0.0D0/
C***********************************************************************
C READ AND ECHO MATERIAL PROPERTIES FOR VON MISES TYPE ELASTO-PLASTIC
C MATERIAL WITH NON-LINEAR (PIECEWISE LINEAR) MIXED HARDENING
C
C REFERENCE: Section 7.6.1
C***********************************************************************
 1000 FORMAT(' VON MISES elasto-plastic model with mixed hardening'/)
 1100 FORMAT(
     1' Mass density ...................................... =',G15.6/
     2' Young''s modulus ................................... =',G15.6/
     3' Poisson''s ratio ................................... =',G15.6)
 1200 FORMAT(/
     1' Number of points on hardening curves .............. =',I3//
     2'          Epbar        isotr.hard. stress ',
     3'   kin.hard. stress'/)
 1300 FORMAT(3(5X,G15.6))
C
C Model currently implemented for plane strain and axisymmetric states
C only
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('ED0156')
C Only small strain implementation is currently available
      IF(NLARGE.EQ.1)CALL ERRPRT('ED0157')
C
C Set unsymmetric tangent stiffness flag
      UNSYM=.FALSE.
C Read and echo some of the real properties
      WRITE(16,1000)
      READ(15,*,ERR=100,END=100)DENSE
      READ(15,*,ERR=100,END=100)YOUNG,POISS
      WRITE(16,1100)DENSE,YOUNG,POISS
      IF(YOUNG.LT.R0)CALL ERRPRT('ED0158')
C number of points on hardening curve
      READ(15,*,ERR=100,END=100)NHARD
      WRITE(16,1200)NHARD
      IF(NHARD.LT.2) CALL ERRPRT('ED0159')
C check dimensions of IPROPS
      IF(MIPROP.LT.3)CALL ERRPRT('ED0160')
      IPROPS(3)=NHARD
C check dimensions of RPROPS
      NRPROP=IPHARD+NHARD*4-1
      IF(NRPROP.GT.MRPROP)CALL ERRPRT('ED0161')
C
      RPROPS(1)=DENSE
      RPROPS(2)=YOUNG
      RPROPS(3)=POISS
C Read and set isotropic and kinematic hardening curves
      IPIHAR=IPHARD
      IPKHAR=IPHARD+2*NHARD
      DO 10 IHARD=1,NHARD
        READ(15,*,ERR=100,END=100)RPROPS(IPIHAR+IHARD*2-2),
     1                            RPROPS(IPIHAR+IHARD*2-1),
     2                            RPROPS(IPKHAR+IHARD*2-1)
        RPROPS(IPKHAR+IHARD*2-2)=RPROPS(IPIHAR+IHARD*2-2)
        WRITE(16,1300)RPROPS(IPHARD+IHARD*2-2),RPROPS(IPHARD+IHARD*2-1),
     1                                         RPROPS(IPKHAR+IHARD*2-1)
   10 CONTINUE
C Check dimension of RSTAVA and LALGVA
      IF(NRSTAV.GT.MRSTAV)CALL ERRPRT('ED0162')
      IF(NLALGV.GT.MLALGV)CALL ERRPRT('ED0163')
C
      GOTO 200
C Issue error message and abort program execution in case of I/O error
  100 CALL ERRPRT('ED0212')
C
  200 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RDVMMX
