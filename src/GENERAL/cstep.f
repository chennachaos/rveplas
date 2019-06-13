CDOC BEGIN_SUBROUTINE CSTEP
CDOC Consistent spatial tangent modulus for finite elasto-plastic models
CDOC
CDOC This routine computes the spatial tangent modulus, a, for
CDOC hyperelastic based (logarithmic strain-based) finite
CDOC elasto-plasticity models in 3-D and 2-D: Plane strain, plane stress
CDOC and axisymmetric states.
CDOC
CDOC BEGIN_PARAMETERS
CDOC DOUBLE_PRECISION AMATX  <  Matrix of components of the spatial
CDOC C                          tangent modulus, a.
CDOC DOUBLE_PRECISION BETRL  >  Array of components of the elastic trial
CDOC C                          left Cauchy-Green strain tensor.
CDOC DOUBLE_PRECISION DMATX  >  Array of components of the infinitesimal
CDOC C                          consistent tangent modulus.
CDOC DOUBLE_PRECISION STRES  >  Array of Cauchy stress tensor
CDOC C                          components
CDOC DOUBLE_PRECISION DETF   >  Determinant of the current deformation
CDOC C                          gradient.
CDOC INTEGER          NTYPE  >  Stress state type flag.
CDOC END_PARAMETERS
CHST
CHST E.de Souza Neto, August 1996: Initial coding
CHST D. de Bortoli  , March  2015: 3-D case added 
CHST                               (renamed from CSTEP2 to CSTEP)
CHST
      SUBROUTINE CSTEP
     1(   AMATX      ,BETRL      ,DMATX      ,STRES      ,DETF       ,
     2    NTYPE       )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER
     1(   MADIM=9    ,MSTRE=6    )
      EXTERNAL  DDLGD2 ,DLGD2
      LOGICAL   OUTOFP ,IS2D
      DIMENSION
     1    AMATX(MADIM,MADIM)        ,BETRL(MSTRE)              ,
     2    DMATX(MSTRE,MSTRE)        ,STRES(MSTRE)
      DIMENSION
     1    AUXMTX(MADIM,MADIM)       ,BMTX(MADIM,MADIM)         ,
     2    DLGAUX(MSTRE,MSTRE)       ,DLGMTX(MADIM,MADIM)       ,
     3    DMATX2(MADIM,MADIM)       ,IG(5)                     ,
     4    IG3D(9)
      DATA
     1    R1   ,R2    /
     2    1.0D0,2.0D0 /
C 2-D indexing from DMATX ordering to AMATX ordering:
C              from (11,22,12,33) to (11,21,12,22,33)
      DATA
     1    IG(1),IG(2),IG(3),IG(4),IG(5)  /
     2    1    ,3    ,3    ,2    ,4      /
C 3-D indexing from DMATX ordering to AMATX ordering:
C              from (11,22,33,12,23,13) to (11,21,31,12,22,32,13,23,33)
      DATA
     1    IG3D(1),IG3D(2),IG3D(3),IG3D(4),IG3D(5),IG3D(6),IG3D(7)  /
     2    1      ,4      ,6      ,4      ,2      ,5      ,6        /
     3    IG3D(8),IG3D(9)  /
     4    5      ,3        /
C***********************************************************************
C COMPUTE THE CONSISTENT SPATIAL TANGENT MODULUS 'a' FOR LARGE STRAIN
C HYPERELASTIC-BASED ELASTOPLASTIC MATERIAL MODELS
C
C REFERENCE: Section 14.5
C***********************************************************************
      IF(NTYPE.EQ.3)THEN
        OUTOFP=.TRUE.
        NADIM=5
        NSTRE=4
        IS2D=.TRUE.
      ELSEIF(NTYPE.EQ.1.OR.NTYPE.EQ.2)THEN
        OUTOFP=.FALSE.
        NADIM=4
        NSTRE=4
        IS2D=.TRUE.
      ELSEIF(NTYPE.EQ.4)THEN
        OUTOFP=.FALSE.
        NADIM=9
        NSTRE=6
        IS2D=.FALSE.
      ELSE
        CALL ERRPRT('EI0020')
      ENDIF
C                           e trial   e trial
C Compute the derivative  dE       /dB
      IF(IS2D)THEN
        CALL DISO2
     1(   DLGAUX(1:4,1:4),DDLGD2     ,DLGD2      ,OUTOFP     ,BETRL    )
      ELSE
        CALL DISO3
     1(   DLGAUX(1:6,1:6),DDLGD2     ,DLGD2      ,BETRL      )
      ENDIF
C
      FACTOR=R1/DETF
      DO 20 I=1,NSTRE
        DO 10 J=1,NSTRE
          DLGAUX(I,J)=FACTOR*DLGAUX(I,J)
   10   CONTINUE
   20 CONTINUE
C                                      e trial   e trial
C Rearrange components of DMATX and  dE       /dB         into the
C ordering compatible with the discrete gradient G, that is:
C   - in 2-D: (11,21,12,22,33)
C   - in 3-D: (11,21,31,12,22,32,13,23,33)
      CALL RVZERO(DMATX2,MADIM*MADIM)
      CALL RVZERO(DLGMTX,MADIM*MADIM)
C
      IF(IS2D)THEN
        DO 40 INEW=1,NADIM
          DO 30 JNEW=1,NADIM
            IOLD=IG(INEW)
            JOLD=IG(JNEW)
            DLGMTX(INEW,JNEW)=DLGAUX(IOLD,JOLD)
            DMATX2(INEW,JNEW)=DMATX(IOLD,JOLD)
   30     CONTINUE
   40   CONTINUE
      ELSE
        DO 41 INEW=1,NADIM
          DO 31 JNEW=1,NADIM
            IOLD=IG3D(INEW)
            JOLD=IG3D(JNEW)
            DLGMTX(INEW,JNEW)=DLGAUX(IOLD,JOLD)
            DMATX2(INEW,JNEW)=DMATX(IOLD,JOLD)
   31     CONTINUE
   41   CONTINUE
      ENDIF
C Compute remaining needed matrix [DELTA_ik BETRL_jl+DELTA_jk BETRL_il]
C (expression 14.102)
      CALL RVZERO(BMTX,MADIM*MADIM)
      IF(IS2D)THEN
        BMTX(1,1)=R2*BETRL(1)
        BMTX(1,3)=R2*BETRL(3)
        BMTX(2,1)=BETRL(3)
        BMTX(2,2)=BETRL(1)
        BMTX(2,3)=BETRL(2)
        BMTX(2,4)=BETRL(3)
        BMTX(3,1)=BETRL(3)
        BMTX(3,2)=BETRL(1)
        BMTX(3,3)=BETRL(2)
        BMTX(3,4)=BETRL(3)
        BMTX(4,2)=R2*BETRL(3)
        BMTX(4,4)=R2*BETRL(2)
        IF(OUTOFP)BMTX(5,5)=R2*BETRL(4)
      ELSE
C                              ijkl                   B_ijlk
C                              -------------------------------------
        BMTX(1,1)=R2*BETRL(1) !1111  -> i==j==k    -> 2*Be_il=2*Be_11
C       BMTX(1,2)=R0          !1121  -> i/=j, j/=k ->   0
C       BMTX(1,3)=R0          !1131  -> i/=j, j/=k ->   0
        BMTX(1,4)=R2*BETRL(4) !1112  -> i==j==k    -> 2*Be_il=2*Be_12
C       BMTX(1,5)=R0          !1122  -> i/=j, j/=k ->   0
C       BMTX(1,6)=R0          !1132  -> i/=j, j/=k ->   0
        BMTX(1,7)=R2*BETRL(6) !1113  -> i==j==k    -> 2*Be_il=2*Be_13
C       BMTX(1,8)=R0          !1123  -> i/=j, j/=k ->   0
C       BMTX(1,9)=R0          !1133  -> i/=j, j/=k ->   0
C
        BMTX(2,1)=BETRL(4)    !2111  -> i/=k, j==k -> Be_il=Be_21
        BMTX(2,2)=BETRL(1)    !2121  -> i==k, j/=k -> Be_jl=Be_11
C       BMTX(2,3)=R0          !2131  -> i/=j, j/=k ->   0
        BMTX(2,4)=BETRL(2)    !2112  -> i/=k, j==k -> Be_il=Be_22
        BMTX(2,5)=BETRL(4)    !2122  -> i==k, j/=k -> Be_jl=Be_12
C       BMTX(2,6)=R0          !2132  -> i/=j, j/=k ->   0
        BMTX(2,7)=BETRL(5)    !2113  -> i/=k, j==k -> Be_il=Be_23
        BMTX(2,8)=BETRL(6)    !2123  -> i==k, j/=k -> Be_jl=Be_13
C       BMTX(2,9)=R0          !2133  -> i/=j, j/=k ->   0
C
        BMTX(3,1)=BETRL(6)    !3111  -> i/=k, j==k -> Be_il=Be_31
C       BMTX(3,2)=R0          !3121  -> i/=j, j/=k ->   0
        BMTX(3,3)=BETRL(1)    !3131  -> i==k, j/=k -> Be_jl=Be_11
        BMTX(3,4)=BETRL(5)    !3112  -> i/=k, j==k -> Be_il=Be_32
C       BMTX(3,5)=R0          !3122  -> i/=j, j/=k ->   0
        BMTX(3,6)=BETRL(4)    !3132  -> i==k, j/=k -> Be_jl=Be_12
        BMTX(3,7)=BETRL(3)    !3113  -> i/=k, j==k -> Be_il=Be_33
C       BMTX(3,8)=R0          !3123  -> i/=j, j/=k ->   0
        BMTX(3,9)=BETRL(6)    !3133  -> i==k, j/=k -> Be_jl=Be_13
C
        BMTX(4,1)=BETRL(4)    !1211  -> i==k, j/=k -> Be_jl=Be_21
        BMTX(4,2)=BETRL(1)    !1221  -> i/=k, j==k -> Be_il=Be_11
C       BMTX(4,3)=R0          !1231  -> i/=j, j/=k ->   0
        BMTX(4,4)=BETRL(2)    !1212  -> i==k, j/=k -> Be_jl=Be_22
        BMTX(4,5)=BETRL(4)    !1222  -> i/=k, j==k -> Be_il=Be_12
C       BMTX(4,6)=R0          !1232  -> i/=j, j/=k ->   0
        BMTX(4,7)=BETRL(5)    !1213  -> i==k, j/=k -> Be_jl=Be_23
        BMTX(4,8)=BETRL(6)    !1223  -> i/=k, j==k -> Be_il=Be_13
C       BMTX(4,9)=R0          !1233  -> i/=j, j/=k ->   0
C
C       BMTX(5,1)=R0          !2211  -> i/=j, j/=k ->   0
        BMTX(5,2)=R2*BETRL(4) !2221  -> i==j==k    -> 2*Be_il=2*Be_21
C       BMTX(5,3)=R0          !2231  -> i/=j, j/=k ->   0
C       BMTX(5,4)=R0          !2212  -> i/=j, j/=k ->   0
        BMTX(5,5)=R2*BETRL(2) !2222  -> i==j==k    -> 2*Be_il=2*Be_22
C       BMTX(5,6)=R0          !2232  -> i/=j, j/=k ->   0
C       BMTX(5,7)=R0          !2213  -> i/=j, j/=k ->   0
        BMTX(5,8)=R2*BETRL(5) !2223  -> i==j==k    -> 2*Be_il=2*Be_23
C       BMTX(5,9)=R0          !2233  -> i/=j, j/=k ->   0
C
C       BMTX(6,1)=R0          !3211  -> i/=j, j/=k ->   0
        BMTX(6,2)=BETRL(6)    !3221  -> i/=k, j==k -> Be_il=Be_31
        BMTX(6,3)=BETRL(4)    !3231  -> i==k, j/=k -> Be_jl=Be_21
C       BMTX(6,4)=R0          !3212  -> i/=j, j/=k ->   0
        BMTX(6,5)=BETRL(5)    !3222  -> i/=k, j==k -> Be_il=Be_32
        BMTX(6,6)=BETRL(2)    !3232  -> i==k, j/=k -> Be_jl=Be_22
C       BMTX(6,7)=R0          !3213  -> i/=j, j/=k ->   0
        BMTX(6,8)=BETRL(3)    !3223  -> i/=k, j==k -> Be_il=Be_33
        BMTX(6,9)=BETRL(5)    !3233  -> i==k, j/=k -> Be_jl=Be_23
C
        BMTX(7,1)=BETRL(6)    !1311  -> i==k, j/=k -> Be_jl=Be_31
C       BMTX(7,2)=R0          !1321  -> i/=j, j/=k ->   0
        BMTX(7,3)=BETRL(1)    !1331  -> i/=k, j==k -> Be_il=Be_11
        BMTX(7,4)=BETRL(5)    !1312  -> i==k, j/=k -> Be_jl=Be_32
C       BMTX(7,5)=R0          !1322  -> i/=j, j/=k ->   0
        BMTX(7,6)=BETRL(4)    !1332  -> i/=k, j==k -> Be_il=Be_12
        BMTX(7,7)=BETRL(3)    !1313  -> i==k, j/=k -> Be_jl=Be_33
C       BMTX(7,8)=R0          !1323  -> i/=j, j/=k ->   0
        BMTX(7,9)=BETRL(6)    !1333  -> i/=k, j==k -> Be_il=Be_13
C
C       BMTX(8,1)=R0          !2311  -> i/=j, j/=k ->   0
        BMTX(8,2)=BETRL(6)    !2321  -> i==k, j/=k -> Be_jl=Be_31
        BMTX(8,3)=BETRL(4)    !2331  -> i/=k, j==k -> Be_il=Be_21
C       BMTX(8,4)=R0          !2312  -> i/=j, j/=k ->   0
        BMTX(8,5)=BETRL(5)    !2322  -> i==k, j/=k -> Be_jl=Be_32
        BMTX(8,6)=BETRL(2)    !2332  -> i/=k, j==k -> Be_il=Be_22
C       BMTX(8,7)=R0          !2313  -> i/=j, j/=k ->   0
        BMTX(8,8)=BETRL(3)    !2323  -> i==k, j/=k -> Be_jl=Be_33
        BMTX(8,9)=BETRL(5)    !2333  -> i/=k, j==k -> Be_il=Be_23
C
C       BMTX(9,1)=R0          !3311  -> i/=j, j/=k ->   0
C       BMTX(9,2)=R0          !3321  -> i/=j, j/=k ->   0
        BMTX(9,3)=R2*BETRL(6) !3331  -> i==j==k    -> 2*Be_il=2*Be_31
C       BMTX(9,4)=R0          !3312  -> i/=j, j/=k ->   0
C       BMTX(9,5)=R0          !3322  -> i/=j, j/=k ->   0
        BMTX(9,6)=R2*BETRL(5) !3332  -> i==j==k    -> 2*Be_il=2*Be_32
C       BMTX(9,7)=R0          !3313  -> i/=j, j/=k ->   0
C       BMTX(9,8)=R0          !3323  -> i/=j, j/=k ->   0
        BMTX(9,9)=R2*BETRL(3) !3333  -> i==j==k    -> 2*Be_il=2*Be_33
      ENDIF
C Assemble the spatial tangent modulus a
C --------------------------------------
C compute the product  D:L:B
      CALL RVZERO(AUXMTX,MADIM*MADIM)
      DO 70 I=1,NADIM
        DO 60 J=1,NADIM
          DO 50 K=1,NADIM
            AUXMTX(I,J)=AUXMTX(I,J)+DMATX2(I,K)*DLGMTX(K,J)
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
      CALL RVZERO(AMATX,MADIM*MADIM)
      DO 100 I=1,NADIM
        DO 90 J=1,NADIM
          DO 80 K=1,NADIM
            AMATX(I,J)=AMATX(I,J)+AUXMTX(I,K)*BMTX(K,J)
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
C subtract  [SIGMA_il DELTA_jk]
      IF(IS2D)THEN
        AMATX(1,1)=AMATX(1,1)-STRES(1)
        AMATX(1,3)=AMATX(1,3)-STRES(3)
        AMATX(2,1)=AMATX(2,1)-STRES(3)
        AMATX(2,3)=AMATX(2,3)-STRES(2)
        AMATX(3,2)=AMATX(3,2)-STRES(1)
        AMATX(3,4)=AMATX(3,4)-STRES(3)
        AMATX(4,2)=AMATX(4,2)-STRES(3)
        AMATX(4,4)=AMATX(4,4)-STRES(2)
        IF(OUTOFP)AMATX(5,5)=AMATX(5,5)-STRES(4)
      ELSE
        AMATX(1,1)=AMATX(1,1)-STRES(1) !1111 -> =AMATX_1111 - STRES_11
        AMATX(1,4)=AMATX(1,4)-STRES(4) !1112 -> =AMATX_1112 - STRES_12
        AMATX(1,7)=AMATX(1,7)-STRES(6) !1113 -> =AMATX_1113 - STRES_13
        AMATX(2,1)=AMATX(2,1)-STRES(4) !2111 -> =AMATX_2111 - STRES_21
        AMATX(2,4)=AMATX(2,4)-STRES(2) !2112 -> =AMATX_2112 - STRES_22
        AMATX(2,7)=AMATX(2,7)-STRES(5) !2113 -> =AMATX_2113 - STRES_23
        AMATX(3,1)=AMATX(3,1)-STRES(6) !3111 -> =AMATX_3111 - STRES_31
        AMATX(3,4)=AMATX(3,4)-STRES(5) !3112 -> =AMATX_3112 - STRES_32
        AMATX(3,7)=AMATX(3,7)-STRES(3) !3113 -> =AMATX_3113 - STRES_33
C
        AMATX(4,2)=AMATX(4,2)-STRES(1) !1221 -> =AMATX_1221 - STRES_11
        AMATX(4,5)=AMATX(4,5)-STRES(4) !1222 -> =AMATX_1222 - STRES_12
        AMATX(4,8)=AMATX(4,8)-STRES(6) !1223 -> =AMATX_1223 - STRES_13
        AMATX(5,2)=AMATX(5,2)-STRES(4) !2221 -> =AMATX_2221 - STRES_21
        AMATX(5,5)=AMATX(5,5)-STRES(2) !2222 -> =AMATX_2222 - STRES_22
        AMATX(5,8)=AMATX(5,8)-STRES(5) !2223 -> =AMATX_2223 - STRES_23
        AMATX(6,2)=AMATX(6,2)-STRES(6) !3221 -> =AMATX_3221 - STRES_31
        AMATX(6,5)=AMATX(6,5)-STRES(5) !3222 -> =AMATX_3222 - STRES_32
        AMATX(6,8)=AMATX(6,8)-STRES(3) !3223 -> =AMATX_3223 - STRES_33
C
        AMATX(7,3)=AMATX(7,3)-STRES(1) !1331 -> =AMATX_1331 - STRES_11
        AMATX(7,6)=AMATX(7,6)-STRES(4) !1332 -> =AMATX_1332 - STRES_12
        AMATX(7,9)=AMATX(7,9)-STRES(6) !1333 -> =AMATX_1333 - STRES_13
        AMATX(8,3)=AMATX(8,3)-STRES(4) !2331 -> =AMATX_2331 - STRES_21
        AMATX(8,6)=AMATX(8,6)-STRES(2) !2332 -> =AMATX_2332 - STRES_22
        AMATX(8,9)=AMATX(8,9)-STRES(5) !2333 -> =AMATX_2333 - STRES_23
        AMATX(9,3)=AMATX(9,3)-STRES(6) !3331 -> =AMATX_3331 - STRES_31
        AMATX(9,6)=AMATX(9,6)-STRES(5) !3332 -> =AMATX_3332 - STRES_32
        AMATX(9,9)=AMATX(9,9)-STRES(3) !3333 -> =AMATX_3333 - STRES_33
      ENDIF
      RETURN
      END
CDOC END_SUBROUTINE CSTEP
