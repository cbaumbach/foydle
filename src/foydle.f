C     ==================================================================

      SUBROUTINE RVAL(XMAT, YMAT, ZC, XCOL, YCOL, N, R, CORES)

      INTEGER XCOL, YCOL, N, CORES
      DOUBLE PRECISION XMAT(XCOL * N), YMAT(YCOL * N), ZC(N)
      DOUBLE PRECISION R(XCOL * YCOL)
      DOUBLE PRECISION X(N), Y(N), Z(N), YZ(N), NORM, PROD
      DOUBLE PRECISION YY, ZZ, XNORM, YZNORM

C$OMP PARALLEL DO
C$OMP& PRIVATE(X, Y, Z, YY, YZ, ZZ, XNORM, YZNORM),
C$OMP& NUM_THREADS(CORES)
      DO 20, J = 1, YCOL
         CALL CPYCOL(YMAT, J, N, Y)
         CALL CPYCOL(ZC, 1, N, Z)
         CALL MULT(Y, Z, YZ, N)
         CALL CENTER(YZ, N, 1)
         YY = PROD(Y, Y, N)
         CALL ORTHO(Z, Y, YY, N)
         ZZ = PROD(Z, Z, N)
         CALL ORTHO(YZ, Y, YY, N)
         CALL ORTHO(YZ, Z, ZZ, N)
         YZNORM = NORM(YZ, N)
         DO 10, I = 1, XCOL
            CALL CPYCOL(XMAT, I, N, X)
            CALL ORTHO(X, Y, YY, N)
            CALL ORTHO(X, Z, ZZ, N)
            XNORM = NORM(X, N)
            R((J - 1) * XCOL + I) = PROD(X, YZ, N) / XNORM / YZNORM
 10      CONTINUE
 20   CONTINUE

      RETURN
      END

C     ==================================================================

      SUBROUTINE ORTHO(X, Y, YY, N)

      INTEGER N
      DOUBLE PRECISION X(N), Y(N), XY, YY

      XY = 0.0
      DO 10, I = 1, N
         XY = XY + X(I) * Y(I)
 10   CONTINUE

      DO 20, I = 1, N
         X(I) = X(I) - XY / YY * Y(I)
 20   CONTINUE

      RETURN
      END

C     ==================================================================

      DOUBLE PRECISION FUNCTION NORM(X, N)

      INTEGER N
      DOUBLE PRECISION X(N), SS

      SS = 0.0
      DO 10, I = 1, N
         SS = SS + X(I)**2
 10   CONTINUE
      NORM = DSQRT(SS)

      END

C     ==================================================================

      SUBROUTINE CPYCOL(MAT, COL, ROWS, X)

      INTEGER COL, ROWS
      DOUBLE PRECISION MAT(COL * ROWS), X(ROWS)

      DO 10, I = 1, ROWS
         X(I) = MAT((COL - 1) * ROWS + I)
 10   CONTINUE

      RETURN
      END

C     ==================================================================

      SUBROUTINE MULT(X, Y, Z, N)

      INTEGER N
      DOUBLE PRECISION X(N), Y(N), Z(N)

      DO 10, I = 1, N
         Z(I) = X(I) * Y(I)
 10   CONTINUE

      RETURN
      END

C     ==================================================================

      DOUBLE PRECISION FUNCTION PROD(X, Y, N)

      INTEGER N
      DOUBLE PRECISION X(N), Y(N)

      PROD = 0.0
      DO 10, I = 1, N
         PROD = PROD + X(I) * Y(I)
 10   CONTINUE

      END

C     ==================================================================

      SUBROUTINE CENTER(MAT, NROW, NCOL)

      INTEGER NROW, NCOL
      DOUBLE PRECISION MAT(NROW * NCOL), SUM, MEAN

      DO 20, I = 1, NCOL
         SUM = 0.0
         DO 10, J = (I - 1) * NROW + 1, I * NROW
            SUM = SUM + MAT(J)
 10      CONTINUE
         MEAN = SUM / DBLE(NROW)
         DO 11, J = (I - 1) * NROW + 1, I * NROW
            MAT(J) = MAT(J) - MEAN
 11      CONTINUE
 20   CONTINUE

      RETURN
      END
