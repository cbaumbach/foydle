C     ==================================================================

      SUBROUTINE RVAL(XMAT, YMAT, ZC, XCOL, YCOL, N, R, CORES)

      INTEGER XCOL, YCOL, N, CORES
      DOUBLE PRECISION XMAT(XCOL * N), YMAT(YCOL * N), ZC(N)
      DOUBLE PRECISION R(XCOL * YCOL)
      DOUBLE PRECISION X(N), Y(N), Z(N), YZ(N), PROD
      DOUBLE PRECISION YY, ZZ

C$OMP PARALLEL DO PRIVATE(X, Y, Z, YY, YZ, ZZ), NUM_THREADS(CORES)
      DO 20, J = 1, YCOL
         CALL COPYCOL(YMAT, J, N, Y)
         CALL COPYCOL(ZC, 1, N, Z)
         CALL MULT(Y, Z, YZ, N)
         CALL CENTER(Y, N)
         CALL CENTER(Z, N)
         CALL CENTER(YZ, N)
         YY = PROD(Y, Y, N)
         CALL ORTHO(Z, Y, YY, N)
         ZZ = PROD(Z, Z, N)
         CALL ORTHO(YZ, Y, YY, N)
         CALL ORTHO(YZ, Z, ZZ, N)
         CALL NORM(YZ, N)
         DO 10, I = 1, XCOL
            CALL COPYCOL(XMAT, I, N, X)
            CALL CENTER(X, N)
            CALL ORTHO(X, Y, YY, N)
            CALL ORTHO(X, Z, ZZ, N)
            CALL NORM(X, N)
            R((J - 1) * XCOL + I) = PROD(X, YZ, N)
 10      CONTINUE
 20   CONTINUE

      RETURN
      END

C     ==================================================================

      SUBROUTINE CENTER(X, N)

      INTEGER N
      DOUBLE PRECISION X(N), MEAN, SUM

      SUM = 0.0
      DO 10, I = 1, N
         SUM = SUM + X(I)
 10   CONTINUE

      MEAN = SUM / N

      DO 20, I = 1, N
         X(I) = X(I) - MEAN
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

      SUBROUTINE NORM(X, N)

      INTEGER N
      DOUBLE PRECISION X(N), S1, S2

C     Compute sum of squares.
      S1 = 0.0
      DO 10, I = 1, N
         S1 = S1 + X(I)**2
 10   CONTINUE

      S2 = DSQRT(S1)
      DO 20, I = 1, N
         X(I) = X(I) / S2
 20   CONTINUE

      RETURN
      END

C     ==================================================================

      SUBROUTINE COPYCOL(MAT, COL, ROWS, X)

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
