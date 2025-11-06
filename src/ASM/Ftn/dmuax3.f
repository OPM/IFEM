      SUBROUTINE DMUAX3 (A,X,Y)
C
C   Calculate   {Y} = [A] * {X}
C
      DOUBLE PRECISION A(3,3),X(3),Y(3)
C
      Y(1) = A(1,1)*X(1) + A(1,2)*X(2) + A(1,3)*X(3)
      Y(2) = A(2,1)*X(1) + A(2,2)*X(2) + A(2,3)*X(3)
      Y(3) = A(3,1)*X(1) + A(3,2)*X(2) + A(3,3)*X(3)
C
      RETURN
      END
