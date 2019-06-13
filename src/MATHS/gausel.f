      SUBROUTINE GAUSEL
     1(   A          ,B          ,N          ,SINGUL     )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C Arguments
      LOGICAL  SINGUL
      DIMENSION
     1    A(N,N)             ,B(N)
C Local definitions
      DATA
     1    SMALL  /
     2    1.0D-10/
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
      ANORM=SQRT(SCAPRD(A,A,N*N))
C
C=============================
C Start loop over pivot rows |
C=============================
C
      DO 40 IPASS=1,N
C
        IF(ABS(A(IPASS,IPASS))/ANORM.LT.SMALL)THEN
C Current pivot is zero (too small): Perform row swapping
          CALL ROWSWP(A,ANORM,B,IPASS,N,SINGUL)
          IF(SINGUL)THEN
C Matrix is singular: exit without solving the system
            GOTO 999
          ENDIF
        ENDIF
C
C STEP 1.  Divide the entire pivot row by the pivot element to get a 1
C in the diagonal position of the pivot row
C
        PIVOT=A(IPASS,IPASS)
        DO 10 ICOL=1,N
          A(IPASS,ICOL)=A(IPASS,ICOL)/PIVOT
   10   CONTINUE
C... the same for the right hand side vector
        B(IPASS)=B(IPASS)/PIVOT
C
C STEP 2.  Replace each row other than the pivot row by that row plus a
C multiple of the pivot row to get a 0 in the pivot column
C
        DO 30 IROW=1,N
          IF(IROW.NE.IPASS)THEN
            FACTOR=A(IROW,IPASS)
            DO 20 ICOL=1,N
              A(IROW,ICOL)=A(IROW,ICOL)-FACTOR*A(IPASS,ICOL) 
   20       CONTINUE
C... the same for the right hand side vector
            B(IROW)=B(IROW)-FACTOR*B(IPASS)
          ENDIF
   30   CONTINUE
C
   40 CONTINUE
C
C==============================
C End of loop over pivot rows |
C==============================
C
C Now, "A" is the IDENTITY matrix and "B" is the solution
C vector.  Print solution vector.
C
  999 CONTINUE
      RETURN
      END




      SUBROUTINE ROWSWP
     1(   A          ,ANORM      ,B          ,I          ,N          ,
     2    SINGUL     )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C Arguments
      LOGICAL  SINGUL
      DIMENSION
     1    A(N,N)             ,B(N)
C Local definitions
      DATA
     1    SMALL  /
     2    1.0d-10/
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
   10     CONTINUE
          BTMP=B(I)
          B(I)=B(J)
          B(J)=BTMP
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
