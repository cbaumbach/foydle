C     ==================================================================

      SUBROUTINE RVAL(XS, YS, ZS, N, R)

      INTEGER N
      DOUBLE PRECISION XS(N), YS(N), ZS(N), R
      DOUBLE PRECISION X(N), Y(N), Z(N), YZ(N)

C     Copy variables.
      DO 10, I = 1, N
         X(I) = XS(I)
         Y(I) = YS(I)
         Z(I) = ZS(I)
 10   CONTINUE

      DO 20, I = 1, N
         YZ(I) = Y(I) * Z(I)
 20   CONTINUE

      CALL CENTER(X, N)
      CALL CENTER(Y, N)
      CALL CENTER(Z, N)
      CALL CENTER(YZ, N)

      CALL ORTHO(Z, Y, N)
      CALL ORTHO(YZ, Y, N)
      CALL ORTHO(YZ, Z, N)
      CALL ORTHO(X, Y, N)
      CALL ORTHO(X, Z, N)

      CALL NORM(YZ, N)
      CALL NORM(X, N)

      R = 0.0
      DO 30, I = 1, N
         R = R + X(I) * YZ(I)
 30   CONTINUE

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

      SUBROUTINE ORTHO(X, Y, N)

      INTEGER N
      DOUBLE PRECISION X(N), Y(N), XY, YY

      XY = 0.0
      DO 10, I = 1, N
         XY = XY + X(I) * Y(I)
 10   CONTINUE

      YY = 0.0
      DO 20, I = 1, N
         YY = YY + Y(I) * Y(I)
 20   CONTINUE

      DO 30, I = 1, N
         X(I) = X(I) - XY / YY * Y(I)
 30   CONTINUE

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
