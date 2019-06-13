C 
C
C Originally GAUSEL
C WARNING: use only for small matrices
C A is M by M, but only the first N by N values are considered in the
C subroutine
      SUBROUTINE GAUSIN
     1(   A          ,AINV            ,N          ,M       ,SINGUL     )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C Arguments
      LOGICAL  SINGUL
      DIMENSION
     1    A(M,M)             ,ACOPY(M,M)          ,AINV(M,M)
C Local definitions
      DATA
     1    SMALL  , R1   , R0   /
     2    1.0D-10, 1.0D0, 0.0D0/
C***********************************************************************
C Routine for solution of simultaneous linear algebraic equations:
C
C                                A x = b
C
C by the GAUSS ELIMINATION method WITH ROW SWAPPING.
C***********************************************************************
C Initialise singular matrix flag
      SINGUL=.FALSE.
C Get norm of matrix A
      ANORM=NORM2(A(1:N,1:N))
C
C Initialise identity matrix
      AINV=R0
      DO I=1,N
        AINV(I,I)=R1
      ENDDO
C Create a copy of matrix A
      ACOPY=A
C=============================
C Start loop over pivot rows |
C=============================
C
      DO 40 IPASS=1,N
C
        IF(ABS(ACOPY(IPASS,IPASS))/ANORM.LT.SMALL)THEN
C Current pivot is zero (too small): Perform row swapping
          CALL ROWSWP2(ACOPY,ANORM,AINV,IPASS,N,M,SINGUL)
          IF(SINGUL)THEN
C Matrix is singular: exit without solving the system
            GOTO 999
          ENDIF
        ENDIF
C
C STEP 1.  Divide the entire pivot row by the pivot element to get a 1
C in the diagonal position of the pivot row
C
        PIVOT=ACOPY(IPASS,IPASS)
        DO 10 ICOL=1,N
          ACOPY(IPASS,ICOL)=ACOPY(IPASS,ICOL)/PIVOT
C... the same for the inverse
          AINV(IPASS,ICOL)=AINV(IPASS,ICOL)/PIVOT
   10   CONTINUE
C
C STEP 2.  Replace each row other than the pivot row by that row plus a
C multiple of the pivot row to get a 0 in the pivot column
C
        DO 30 IROW=1,N
          IF(IROW.NE.IPASS)THEN
            FACTOR=ACOPY(IROW,IPASS)
            DO 20 ICOL=1,N
              ACOPY(IROW,ICOL)=ACOPY(IROW,ICOL)-FACTOR*ACOPY(IPASS,ICOL)
C... the same for the right hand side vector
              AINV(IROW,ICOL)=AINV(IROW,ICOL)-FACTOR*AINV(IPASS,ICOL) 
   20       CONTINUE
          ENDIF
   30   CONTINUE
C
   40 CONTINUE
C
C==============================
C End of loop over pivot rows |
C==============================
C
C Now, "ACOPY" is the IDENTITY matrix and "AINV" is its inverse.
C
  999 CONTINUE
      RETURN
      END



C Swaps rows of A and AINV to avoid zero pivots. If a non-zero pivot is
C not available, flag the matrix as singular.
      SUBROUTINE ROWSWP2
     1(   A          ,ANORM      ,AINV         ,I          ,N          ,
     2    M          ,SINGUL     )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C Arguments
      LOGICAL  SINGUL
      DIMENSION
     1    A(M,M)             ,AINV(M,M)
C Local definitions
      DATA
     1    SMALL  /
     2    1.0D-10/
C***********************************************************************
C If row "i" has a zero (sufficiently small) entry on the diagonal
C (column "i") then swap that row with the next row below row "i" with
C a non-zero (sufficiently large) entry in column "i"
C***********************************************************************
      DO 20 J=I,N
        IF(ABS(A(J,I))/ANORM.GE.SMALL)THEN
C Non-zero (sufficiently large) element found: Swap rows "i" and "j"
          DO 10 K=1,N
            ATMP=A(I,K)
            A(I,K)=A(J,K)
            A(J,K)=ATMP
C
            ATMP=AINV(I,K)
            AINV(I,K)=AINV(J,K)
            AINV(J,K)=ATMP
   10     CONTINUE
C Row swapping complete: break loop and exit successfully
          SINGUL=.FALSE.
          GOTO 30
        ENDIF
   20 CONTINUE
C Row swapping failed (no non-zero element found): matrix is singular
      SINGUL=.TRUE.
C
   30 CONTINUE
      RETURN
      END
