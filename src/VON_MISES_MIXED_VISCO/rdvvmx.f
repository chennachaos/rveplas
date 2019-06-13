CDOC BEGIN_SUBROUTINE RDVVMX
CDOC Read data for a von Mises viscoplastic model with mixed hardening.
CDOC
CDOC This routine reads from the data file and echos to the results file
CDOC the material parameters necessary for
CDOC the von Mises elasto-viscoplastic model with piece-wise linear
CDOC mixed isotropic/kinematic hardening and Peric's power law for
CDOC viscoplasticity.
CDOC It also sets the array of real properties and some components of
CDOC the array of integer material properties.
CDOC These arrays are used by subroutines SUVVMX and CTVVMX.
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
CDOC INTEGER          NLARGE >  Large strain analysis flag.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC DOUBLE_PRECISION RPROPS <  Array of real material properties.
CDOC LOGICAL          UNSYM  <  Unsymmetric tangent stiffness flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto & F.Adziman, September 2012: Initial coding
CHST
      SUBROUTINE RDVVMX
     1(   IPROPS     ,MIPROP     ,MLALGV     ,MRPROP     ,MRSTAV     ,
     2    NLARGE     ,NTYPE      ,RPROPS     ,UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UNSYM
      PARAMETER( IPHARD=6  ,NLALGV=2  ,NRSTAV=9 )
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)
      DATA R0   /0.0D0/
C***********************************************************************
C READ AND ECHO MATERIAL PROPERTIES FOR VON MISES TYPE
C ELASTO-VISCOPLASTIC MATERIAL MODEL WITH NON-LINEAR (PIECEWISE LINEAR)
C MIXED HARDENING AND PERIC'S POWER LAW FOR VISCOPLASTIC FLOW.
C
C REFERENCE: Section 11.6 (this model has a generalisation, which
C            includes mixed hardening, of the original isotropically
C            hardening model of Section 11.6). See also Section 7.6.1
C            for the analogous rate-independent model with mixed
C            hardening.
C***********************************************************************
 1000 FORMAT(' VON MISES viscoplastic Peric power law model'
     1       ' with mixed hardening'/)
 1100 FORMAT(
     1' Mass density ...................................... =',G15.6/
     2' Young''s modulus ................................... =',G15.6/
     3' Poisson''s ratio ................................... =',G15.6/
     4' Viscosity parameter (mu) .......................... =',G15.6/
     5' Rate sensitivity parameter (epsilon) .............. =',G15.6/)
 1200 FORMAT(/
     1' Number of points on hardening curves .............. =',I3//
     2'          Epbar        isotr.hard. stress ',
     3'   kin.hard. stress'/)
 1300 FORMAT(3(5X,G15.6))
C
C Model currently implemented for plane strain and axisymmetric states
C only
      IF(NTYPE.NE.2.AND.NTYPE.NE.3)CALL ERRPRT('ED0228')
C Only small strain implementation is currently available
      IF(NLARGE.EQ.1)CALL ERRPRT('ED0229')
C
C Set unsymmetric tangent stiffness flag
      UNSYM=.FALSE.
C Read and echo some of the real properties
      WRITE(16,1000)
      READ(15,*,ERR=100,END=100)DENSE
      READ(15,*,ERR=100,END=100)YOUNG,POISS,VISCO,RTSEN
      WRITE(16,1100)DENSE,YOUNG,POISS,VISCO,RTSEN
      IF(YOUNG.LT.R0)CALL ERRPRT('ED0230')
C number of points on hardening curve
      READ(15,*,ERR=100,END=100)NHARD
      WRITE(16,1200)NHARD
      IF(NHARD.LT.2) CALL ERRPRT('ED0231')
C check dimensions of IPROPS
      IF(MIPROP.LT.3)CALL ERRPRT('ED0232')
      IPROPS(3)=NHARD
C check dimensions of RPROPS
      NRPROP=IPHARD+NHARD*4-1
      IF(NRPROP.GT.MRPROP)CALL ERRPRT('ED0233')
C
      RPROPS(1)=DENSE
      RPROPS(2)=YOUNG
      RPROPS(3)=POISS
      RPROPS(4)=VISCO
      RPROPS(5)=RTSEN
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
      IF(NRSTAV.GT.MRSTAV)CALL ERRPRT('ED0234')
      IF(NLALGV.GT.MLALGV)CALL ERRPRT('ED0235')
C
      GOTO 200
C Issue error message and abort program execution in case of I/O error
  100 CALL ERRPRT('ED0236')
C
  200 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RDVVMX
