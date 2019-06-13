CDOC BEGIN_SUBROUTINE RVSTRA
CDOC Reads macroscopic strain input.
CDOC
CDOC This routine reads prescribed array of engineering strain
CDOC for small strains and prescribed array of deformation gradient
CDOC for large strains from the input data file.
CDOC
CHST
CHST M.F. Adziman    ,  July 2013: Initial coding 
CHST                               (replacing INLOAD in HYPLAS)
CHST
      SUBROUTINE RVSTRA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C Hyplas global database
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
      INCLUDE '../RVE.INC'
C
      LOGICAL FOUND, NOSTR, NODFG
      CHARACTER*80 INLINE
      DIMENSION
     1    IWBEG(40)     ,IWEND(40)               
C***********************************************************************
C READS STRAIN LOADINGS (ENGINEERING STRAIN) FROM INPUT DATA FILE 
C
C***********************************************************************
 1000 FORMAT(///
     1' Prescribed (macroscopic) infinitesimal strain tensor'/
     2' ===================================================='//
     3' exx ...........................=',G15.6/
     4' eyy ...........................=',G15.6/
     5' exy ...........................=',G15.6)
 1050 FORMAT(///
     1' Prescribed (macroscopic) infinitesimal strain tensor'/
     2' ===================================================='//
     3' exx ...........................=',G15.6/
     4' eyy ...........................=',G15.6/
     5' ezz ...........................=',G15.6/
     6' exy ...........................=',G15.6/
     7' eyz ...........................=',G15.6/
     8' exz ...........................=',G15.6)
 1100 FORMAT(///
     1' Prescribed (macroscopic) deformation gradient'/
     2' ============================================='//
     3' Fxx ...........................=',G15.6/
     4' Fxy ...........................=',G15.6/
     5' Fyx ...........................=',G15.6/
     6' Fyy ...........................=',G15.6)
 1150 FORMAT(///
     1' Prescribed (macroscopic) deformation gradient'/
     2' ============================================='//
     3' Fxx ...........................=',G15.6/
     4' Fxy ...........................=',G15.6/
     5' Fxz ...........................=',G15.6/
     6' Fyx ...........................=',G15.6/
     7' Fyy ...........................=',G15.6/
     8' Fyz ...........................=',G15.6/
     9' Fzx ...........................=',G15.6/
     O' Fzy ...........................=',G15.6/
     1' Fzz ...........................=',G15.6)
C
C Initialise strain vector and logical conditions 
C ===============================================
C
      CALL RVZERO(DRV_STRAIN,9)
      NOSTR=.FALSE.
      NODFG=.FALSE.
C
C Read prescribed array of engineering strains for small strains
C ==============================================================
C
      CALL FNDKEY
     1(   FOUND   ,IWBEG  ,IWEND  ,'PRESCRIBED_INFINITESIMAL_STRAIN',
     2    INLINE  ,15     ,NWRD   )
      IF(.NOT.FOUND)THEN 
        NOSTR=.TRUE.
        ELSE
          IF(NLARGE.NE.0)CALL ERRPRT('ED0405')
C Read the prescribed array of engineering strains
          IF(NTYPE==4)THEN
            NSTRA=6
            READ(15,*,ERR=900,END=900)DRV_STRAIN(1:NSTRA)
            WRITE(16,1050)DRV_STRAIN(1:NSTRA)
          ELSE
            NSTRA=3
            READ(15,*,ERR=900,END=900)DRV_STRAIN(1:NSTRA)
            WRITE(16,1000)DRV_STRAIN(1:NSTRA)
          ENDIF
      ENDIF
C 
C Read prescribed deformation gradient for large strains
C ======================================================
C
      CALL FNDKEY
     1(   FOUND   ,IWBEG  ,IWEND  ,'PRESCRIBED_DEFORMATION_GRADIENT',
     2    INLINE  ,15     ,NWRD   )
      IF(.NOT.FOUND)THEN 
        NODFG=.TRUE.
        ELSE
          IF(NLARGE.NE.1)CALL ERRPRT('ED0407')
C Read the prescribed deformation gradient
          IF(NTYPE==4)THEN
            NSTRA=9
            READ(15,*,ERR=910,END=910)DRV_STRAIN(1:NSTRA)
            WRITE(16,1150)DRV_STRAIN(1:NSTRA)
          ELSE
            NSTRA=4
            READ(15,*,ERR=910,END=910)DRV_STRAIN(1:NSTRA) 
c            WRITE(16,1100)DRV_STRAIN(1:NSTRA)
          ENDIF
      ENDIF
C
C Check whether macroscopic strains has been inputted  
      IF(NOSTR.AND.NODFG)CALL ERRPRT('ED0403')
C 
      GOTO 950
  900 CALL ERRPRT('ED0406')
  910 CALL ERRPRT('ED0408')
C
  950 CONTINUE 
      RETURN
      END
CDOC END_SUBROUTINE RVSTRA
