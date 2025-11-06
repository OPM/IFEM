      DOUBLE PRECISION FUNCTION DNORM3 (V)
      DOUBLE PRECISION V(3),DV
C
C     DNORM3 : Calculate  |v|  and  1/|v| * {v}
C
      DOUBLE PRECISION   ZERO        , EPS          , ONE
      PARAMETER        ( ZERO = 0.0D0, EPS = 1.0D-15, ONE = 1.0D0 )
C
      DV = SQRT(V(1)*V(1) + V(2)*V(2) + V(3)*V(3))
      IF (DV .LT. EPS) THEN
         DNORM3 = ZERO
      ELSE
         DNORM3 = DV
         V(1)   = V(1) / DV
         V(2)   = V(2) / DV
         V(3)   = V(3) / DV
         IF (ABS(V(1)) .GT. ONE-EPS) THEN
            V(1) = SIGN(ONE,V(1))
            V(2) = ZERO
            V(3) = ZERO
         ELSE IF (ABS(V(2)) .GT. ONE-EPS) THEN
            V(1) = ZERO
            V(2) = SIGN(ONE,V(2))
            V(3) = ZERO
         ELSE IF (ABS(V(3)) .GT. ONE-EPS) THEN
            V(1) = ZERO
            V(2) = ZERO
            V(3) = SIGN(ONE,V(3))
         ENDIF
      ENDIF
C
      RETURN
      END
