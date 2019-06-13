CDOC BEGIN_SUBROUTINE RDOGD
CDOC Reads data for the regularised Ogden hyperelastic material model.
CDOC
CDOC This routine reads from the data file and echoes to the results
CDOC file the material parameters required by the regularised Ogden
CDOC hyperelastic material model.  It also sets the array of real
CDOC properties and some components of the array of integer material
CDOC properties.  These arrays are used by subroutines SUOGD and CTOGD.
CDOC Note that under plane stress the model implemented here is exactly
CDOC incompressible, whereas under plane strain and axisymmetric
CDOC conditions the regularised (penalty type) approach is used to
CDOC enforce incompressibility. This routine also sets the unsymmetric
CDOC tangent stiffness flag.
CDOC
CDOC BEGIN_PARAMETERS
CDOC INTEGER          IPROPS <  Array of integer material properties.
CDOC INTEGER          MIPROP >  Dimension of the global array of integer
CDOC C                          material properties.
CDOC INTEGER          MRPROP >  Dimension of the global array of real
CDOC C                          material variables.
CDOC INTEGER          MRSTAV >  Dimension of the global array of real
CDOC C                          state variables.
CDOC DOUBLE_PRECISION RPROPS <  Array of real material properties.
CDOC LOGICAL          UNSYM  <  Unsymmetric tangent stiffness flag.
CDOC C                          Always returned as .FALSE..
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, June 1996: Initial coding
CHST
CHST E.de Souza Neto, October 2008: Dimensioning checks added
CHST
CHST E.de Souza Neto, April 2011: I/O error message added
CHST
      SUBROUTINE RDOGD
     1(   IPROPS     ,MIPROP     ,MRPROP     ,MRSTAV     ,RPROPS     ,
     2    UNSYM      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL UNSYM
      PARAMETER( IPOGDC=2  ,NIPROP=3  ,NRSTAV=4 )
      DIMENSION
     1    IPROPS(*)          ,RPROPS(*)
C***********************************************************************
C READS AND ECHOES MATERIAL PROPERTIES FOR OGDEN TYPE HYPERELASTIC
C MATERIAL MODEL
C
C REFERENCE: Section 13.2.2
C***********************************************************************
 1000 FORMAT(' Ogden type hyperelastic material'/)
 1010 FORMAT(
     1' Mass density ...................................... =',G15.6)
 1020 FORMAT(/
     1' Number of terms in Ogden''s strain-energy function.. =',I3//
     2'            mu                  alpha'/)
 1030 FORMAT(2(5X,G15.6))
 1040 FORMAT(/
     1' Bulk modulus ...................................... =',G15.6/)
C
C Set unsymmetric tangent stiffness flag
      UNSYM=.FALSE.
C Read and echo some of the real properties
      WRITE(16,1000)
      READ(15,*,ERR=100,END=100)DENSE
      WRITE(16,1010)DENSE
C Ogden's constants
      READ(15,*,ERR=100,END=100)NOGTRM
      WRITE(16,1020)NOGTRM
      IF(NOGTRM.LT.1)CALL ERRPRT('ED0067')
C Check dimension of IPROPS
      IF(MIPROP.LT.NIPROP)CALL ERRPRT('ED0195')
      IPROPS(3)=NOGTRM
C Check dimension of RPROPS
      NRPROP=IPOGDC+NOGTRM*2
      IF(NRPROP.GT.MRPROP)CALL ERRPRT('ED0196')
C Store real properties in RPROPS
      RPROPS(1)=DENSE
      DO 10 I=1,NOGTRM
        READ(15,*,ERR=100,END=100)
     1                RPROPS(IPOGDC+I*2-2),RPROPS(IPOGDC+I*2-1)
        WRITE(16,1030)RPROPS(IPOGDC+I*2-2),RPROPS(IPOGDC+I*2-1)
   10 CONTINUE
C Bulk modulus
      READ(15,*,ERR=100,END=100)BULK
      RPROPS(IPOGDC+NOGTRM*2)=BULK
      WRITE(16,1040)BULK
C Check dimension of RSTAVA
      IF(NRSTAV.GT.MRSTAV)CALL ERRPRT('ED0197')
C
      GOTO 200
C Issue error message and abort program execution in case of I/O error
  100 CALL ERRPRT('ED0209')
C
  200 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RDOGD
