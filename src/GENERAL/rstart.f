CDOC BEGIN_SUBROUTINE RSTART
CDOC Reads (writes) data to input (output) re-start file
CDOC
CDOC Depending on the entry value of its argument MODE, this
CDOC routine either dumps the complete database and current information 
CDOC about the analysis into an output re-start file or reads this
CDOC information from a previously created re-start file.
CDOC The information dumped on the output re-start file will enable
CDOC HYPLAS to re-start from the point in the analysis where the file
CDOC was generated. To re-start, this file will be read as the input
CDOC re-start file.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION DFOLD  <> Incremental load factor of the previous
CDOC C                          load increment.
CDOC DOUBLE_PRECISION DLENG  <> Current arc-length (used for arc-length
CDOC C                          method only).
CDOC DOUBLE_PRECISION DLENGO <> Previous arc-length (used for arc-length
CDOC C                          method only).
CDOC DOUBLE_PRECISION DLENM  <> Current maximum arc-length (used for
CDOC C                          arc-length method only).
CDOC DOUBLE_PRECISION DLAMD  <> Iterative load factor obtained by the
CDOC C                          arc-length method only.
CDOC INTEGER          IFNEG  <> Signum (+1/-1) of the tangent stifness
CDOC C                          matrix determinant (used only by the
CDOC C                          arc-length method).
CDOC INTEGER          IINCS  <> Number of the current global
CDOC C                          equilibrium iteration.
CDOC INTEGER          MXFRON <> Maximum front encoutered.
CDOC INTEGER          NOUTP  >  Array of output flags for the current
CDOC C                          load increment.
CDOC DOUBLE_PRECISION TFACT  <> Current total load factor.
CDOC DOUBLE_PRECISION TFACTO <> Previous total load factor.
CDOC DOUBLE_PRECISION TTIME  <> Current time.
CDOC DOUBLE_PRECISION TTIMEO <> Time at the end of previous increment.
CDOC LOGICAL          UNSYM  <> Global unsymmetric solver flag.
CDOC CHARACTER        RSTINP >  Current input re-start file name.
CDOC CHARACTER        RSTOUT >  Current output re-start file name.
CDOC INTEGER          MODE   >  If MODE=0, reads from current
CDOC C                          input re-start file. If MODE=1,
CDOC C                          writes to the current output re-start
CDOC C                          file.
CDOC INTEGER          INCRST <  Increment number corresponding to the
CDOC C                          current output re-start file.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, September 1996: Initial coding
CHST
CHST E.de Souza Neto, June      1998: Common block ARCLEN removed
CHST                                  and global array DTANG (part of
CHST                                  ARCLEN) moved into common block
CHST                                  RESUL.
CHST
CHST E.de Souza Neto, April     2011: I/O error message added for I/O
CHST                                  error while reading.
CHST
CHST E.de Souza Neto & F.Adziman,
CHST                  September 2012: Time info added.
CHST
CHST M.F. Adziman, D. de Bortoli, July 2013: 
CHST      - Variable NSOLVE related to solver choice added to common 
CHST        block CONTRL.
CHST      - Unused array ANGFIX removed from common block MESH (now only
CHST        ANGLE is to apply prescribed displacements at an angle)
CHST
      SUBROUTINE RSTART
     1(  DFOLD      ,DLENG      ,DLENGO     ,DLENM      ,DLAMD      ,
     2   IFNEG      ,IINCS      ,MXFRON     ,NOUTP      ,TFACT      ,
     3   TFACTO     ,TTIME      ,TTIMEO     ,UNSYM      ,RSTINP     ,
     4   RSTOUT     ,MODE       ,INCRST     )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE '../MAXDIM.INC'
      INCLUDE '../MATERIAL.INC'
      INCLUDE '../ELEMENTS.INC'
      INCLUDE '../GLBDBASE.INC'
C
      DIMENSION NOUTP(5)
      LOGICAL UNSYM
      CHARACTER*256 RSTINP,RSTOUT
C***********************************************************************
C WRITE DATA TO OUTPUT RE-START FILE AND READ DATA FROM INPUT RE-START
C FILE
C
C REFERENCE: Figures 5.1-3
C***********************************************************************
 1000 FORMAT(//15X,
     1       'Writing converged results to re-start file ...'//)
 1010 FORMAT(////15X,
     1       'Reading data from input re-start file ...')
C
C
      IF(MODE.EQ.1)THEN
C
C Writing (output) mode
C =====================
C
C Check re-start output flag
C
        IF(NOUTP(5).EQ.0)GOTO 999
        IF(NOUTP(5).NE.0.AND.NALGO.LT.0)THEN
          IF(MOD(IINCS,NOUTP(5)).NE.0)GOTO 999
        ENDIF
C
        INCRST=IINCS
        WRITE(*,1000)
        WRITE(16,1000)
        WRITE(18,1000)
        OPEN(UNIT=17,FILE=RSTOUT,STATUS='UNKNOWN',FORM='UNFORMATTED')
C
C Write some global variables first
C ---------------------------------
        WRITE(17)DFOLD,DLENG,DLENGO,DLENM,DLAMD,IFNEG,IINCS,MXFRON,
     1           TFACT,TFACTO,TTIME,TTIMEO,UNSYM
C
C Then write all common blocks 
C ----------------------------
C COMMON/CONTRL/
        WRITE(17)
     1    NDOFN      ,NELEM      ,NGRUP      ,NPOIN      ,NTOTV      ,
     2    NVFIX      ,NTYPE      ,NALGO      ,NARCL      ,NDIME      ,
     3    NLARGE     ,NAXIS      ,NSOLVE
C COMMON/CORE  /
        WRITE(17)
     1    FIXED                        ,STFOR                        ,
     2    TOFOR                        ,ELOAD                        ,
     3    ELOADO                       ,RLOAD
C COMMON/MATERL/
        WRITE(17)
     1    RPROPS              ,IPROPS
C COMMON/MESH  /
        WRITE(17)
     1    ANGLE              ,COORD              ,PRESC              ,
     2    IELTID             ,IFFIX              ,IGRPID             ,
     3    LNODS              ,MASTER             ,MATTID             ,
     4    NOFIX              ,NVALEN
C COMMON/ELEMEN/
        WRITE(17)
     1    RELPRP              ,IELPRP
C COMMON/RESULT/
        WRITE(17)
     1    DITER              ,DINCR              ,DINCRO             ,
     2    DTANG              ,TDISP              ,TDISPO             ,
     3    TREAC
C COMMON/STATE /
        WRITE(17)
     1    RALGVA                       ,RSTAVA                       ,
     2    STRSG                        ,THKGP                        ,
     3    LALGVA
C
C
      ELSEIF(MODE.EQ.0)THEN
C
C Reading (input) mode
C ====================
C
        WRITE(*,1010)
        OPEN(UNIT=17,FILE=RSTINP,STATUS='OLD',FORM='UNFORMATTED')
C
C Read some global variables first
C --------------------------------
        READ(17,END=900,ERR=900)DFOLD,DLENG,DLENGO,DLENM,DLAMD,IFNEG,
     1                          IINCS,MXFRON,TFACT,TFACTO,TTIME,TTIMEO,
     2                          UNSYM
C
C Then read all common blocks 
C ---------------------------
C COMMON/CONTRL/
        READ(17,END=900,ERR=900)
     1    NDOFN      ,NELEM      ,NGRUP      ,NPOIN      ,NTOTV      ,
     2    NVFIX      ,NTYPE      ,NALGO      ,NARCL      ,NDIME      ,
     3    NLARGE     ,NAXIS      ,NSOLVE
C COMMON/CORE  /
        READ(17,END=900,ERR=900)
     1    FIXED                        ,STFOR                        ,
     2    TOFOR                        ,ELOAD                        ,
     3    ELOADO                       ,RLOAD
C COMMON/MATERL/
        READ(17,END=900,ERR=900)
     1    RPROPS              ,IPROPS
C COMMON/MESH  /
        READ(17,END=900,ERR=900)
     1    ANGLE              ,COORD              ,PRESC              ,
     2    IELTID             ,IFFIX              ,IGRPID             ,
     3    LNODS              ,MASTER             ,MATTID             ,
     4    NOFIX              ,NVALEN
C COMMON/ELEMEN/
        READ(17,END=900,ERR=900)
     1    RELPRP              ,IELPRP
C COMMON/RESULT/
        READ(17,END=900,ERR=900)
     1    DITER              ,DINCR              ,DINCRO             ,
     2    DTANG              ,TDISP              ,TDISPO             ,
     3    TREAC
C COMMON/STATE /
        READ(17,END=900,ERR=900)
     1    RALGVA                       ,RSTAVA                       ,
     2    STRSG                        ,THKGP                        ,
     3    LALGVA
C
        GOTO 910
C Issue error message and abort program execution in case of I/O error
  900   CLOSE(UNIT=17,STATUS='KEEP')
        CALL ERRPRT('ED0213')
C
  910   CONTINUE
C
      ENDIF
C
      CLOSE(UNIT=17,STATUS='KEEP')
C
  999 CONTINUE
      RETURN
      END
CDOC END_SUBROUTINE RSTART
