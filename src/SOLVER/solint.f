CDOC BEGIN_SUBROUTINE SOLINT
CDOC Interface for the solution of the linear system of equations
CDOC
CDOC This routine calls the corresponding function to perform the linear
CDOC system solution, according to the selected method (currently
CDOC frontal or MA41 - sparse multifrontal method).
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DTIME  >  Time increment.
CDOC INTEGER          IITER  >  Current equilibrium iteration number.
CDOC INTEGER          KRESL  >  Equation resolution index.
CDOC INTEGER          IFNEG  <  Signum (+1/-1) of the determinant of the
CDOC C                          stiffness matrix.
CDOC INTEGER          KUNLD  <> Unloading flag.
CDOC DOUBLE_PRECISION MXFRON >  Maximum frontwidth encountered in the
CDOC C                          system of linear finite element
CDOC C                          equations. Used only by FRONT.
CDOC LOGICAL          UNSYM  >  Stiffness matrix unsymmetry flag.
CDOC LOGICAL          INCCUT <  Load increment cutting flag.
CDOC INTEGER          NSOLVE >  Solver type flag.
CDOC END_PARAMETERS
CHST
CHST F.M. Andrade Pires, February 2002: Initial coding as SOLVER
CHST
CHST M.F. Adziman, D. de Bortoli, E.A. de Souza Neto, July 2013:
CHST      Modified to make it compatible with HYPLAS v3.1
CHST
      SUBROUTINE SOLINT
     1(   DTIME      ,IITER      ,KRESL      ,IFNEG      ,KUNLD      , 
     2    MXFRON     ,UNSYM      ,INCCUT     ,NSOLVE     )
      IMPLICIT NONE
C
C Arguments
C
      DOUBLE PRECISION DTIME
      INTEGER IITER, KRESL, IFNEG, KUNLD, MXFRON, NSOLVE
      LOGICAL UNSYM, INCCUT
      
C
C Call solver routine according to the type of solver
C ---------------------------------------------------
C
      IF(NSOLVE.EQ.1)THEN
C Solves the system of equations by the frontal method
        CALL FRONT ( DTIME  ,IITER ,KRESL ,IFNEG ,KUNLD  ,
     1               MXFRON ,UNSYM ,INCCUT  )
      ELSEIF(NSOLVE.EQ.2)THEN
C Multifrontal sparse Gaussian elimination: MA41 from HSL
        CALL INMA41 ( DTIME, IFNEG, IITER, INCCUT, KRESL, KUNLD,
     1                UNSYM )
      ELSE
C Incorrect solver specification
        CALL ERRPRT('ED0302')
      ENDIF
C
      RETURN
      END
CDOC END_SUBROUTINE SOLINT
