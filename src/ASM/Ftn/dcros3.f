      SUBROUTINE DCROS3 (X,Y,Z)
C
C   Calculate   {Z} = {X} x {Y}
C
      DOUBLE PRECISION X(3),Y(3),Z(3)
C
      Z(1) = X(2)*Y(3) - Y(2)*X(3)
      Z(2) = X(3)*Y(1) - Y(3)*X(1)
      Z(3) = X(1)*Y(2) - Y(1)*X(2)
C
      RETURN
      END
