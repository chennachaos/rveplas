CDOC BEGIN_SUBROUTINE CTTRPN
CDOC Consistent tangent matrix for the Tresca model in plane stress.
CDOC
CDOC The tangent matrix computed in this routine is consistent with the
CDOC nested iteration algorithm for the Tresca model in plane stress
CDOC coded in subroutine SUTRPN.
CDOC It returns either the elastic tangent or the elasto-plastic
CDOC consistent tangent matrix depending on the input value of
CDOC the logical argument EPFLAG.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DMATX  <  Consistent tangent matrix.
CDOC LOGICAL          EPFLAG >  Elasto-plastic flag.
CDOC C                          If .FALSE., DMATX returns as the elastic
CDOC C                          as the elastic matrix. If .TRUE., DMATX
CDOC C                          returns as the elasto-plastic tangent
CDOC C                          consistent with the return mapping
CDOC C                          algorithm implemented in routine SUTR.
CDOC INTEGER          IPROPS >  Array of integer material properties.
CDOC C                          This array is set in routines INDATA
CDOC C                          and RDTR.
CDOC C                          The number of points on the piece-wise
CDOC C                          linear hardening curve is the only
CDOC C                          element of this array used here.
CDOC LOGICAL          LALGVA >  Array of logical algorithmic flags.
CDOC C                          See list of arguments of SUTR.
CDOC DOUBLE_PRECISION RPROPS >  Array of real material properties.
CDOC C                          Same as in the argument list of
CDOC C                          subroutine SUTR.
CDOC DOUBLE_PRECISION RSTAVA >  Array of real state variables other
CDOC C                           than the stress tensor components.
CDOC C                          Output of SUTR.
CDOC DOUBLE_PRECISION STRAT  >  Array of elastic trial (engineering)
CDOC C                          strain components. Same as in the input
CDOC C                          SUTR.
CDOC DOUBLE_PRECISION STRES  >  Array of current stress tensor
CDOC C                          components.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, Sept 1998: Initial coding
CHST
      SUBROUTINE CTTRPN
     1(   DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,RPROPS     ,
     2    RSTAVA     ,STRAT      ,STRES      )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(  MSTRE=4 )
C Arguments
      LOGICAL EPFLAG, LALGVA(4)
      DIMENSION
     1    DMATX(MSTRE,MSTRE) ,IPROPS(*)          ,RPROPS(*)          ,
     2    RSTAVA(MSTRE+1)    ,STRAT(*)           ,STRES(*)
C Local arrays and variables
      DIMENSION
     1    D12(3)             ,D21(3)
C***********************************************************************
C COMPUTATION OF CONSISTENT TANGENT MODULUS FOR TRESCA TYPE
C ELASTO-PLASTIC MATERIAL WITH PIECE-WISE LINEAR ISOTROPIC HARDENING.
C PLANE STRAIN AND AXISYMMETRIC IMPLEMENTATIONS.
C
C REFERENCE: Section 9.2
C***********************************************************************
      STRAT(4)=RSTAVA(4)
C Compute the axisymmetric tangent matrix
      NTYPAX=3
      CALL CTTR
     1(   DMATX      ,EPFLAG     ,IPROPS     ,LALGVA     ,NTYPAX     ,
     2    RPROPS     ,RSTAVA     ,STRAT      ,STRES      )
C Decompose into submatrices
      D12(1)=DMATX(1,4)
      D12(2)=DMATX(2,4)
      D12(3)=DMATX(3,4)
      D21(1)=DMATX(4,1)
      D21(2)=DMATX(4,2)
      D21(3)=DMATX(4,3)
      D22=DMATX(4,4)
C Assemble plane stress consistent tangent matrix
      DO 20 I=1,3
        DO 10 J=1,3
          DMATX(I,J)=DMATX(I,J)-D12(I)*D21(J)/D22
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE CTTRPN
