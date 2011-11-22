      SUBROUTINE DSLBLN (IPSW,IWR,EPS,A,B,NSOL,X,Y)
C
C $Id$
C***************************** USERS' PART *****************************
C
C A F E M  ROUTINE : DSLBLN
C
C PURPOSE:
C     DSLBLN : Solve a bi-linear set of equations in x and y.
C
C METHOD:
C     The following set of equations are solved for x and y:
C
C        A1 * x*y  +  A2 * x  +  A3 * y  =  A4
C        B1 * x*y  +  B2 * x  +  B3 * y  =  B4
C
C INPUT ARGUMENTS:
C     IPSW   - Print switch
C     IWR    - Print unit number
C     EPS    - Relative tolerance for equality check with zero
C     A, B   - Coefficients of the bi-linear equations
C
C INPUT/OUTPUT ARGUMENTS:
C     None
C
C OUTPUT ARGUMENTS:
C     NSOL   - Number of solutions found
C     X, Y   - The solution(s)
C
C RESTRICTIONS:
C     None
C
C************************** MAINTENANCE PART ***************************
C
C INTERNAL VARIABLES:
C     TOL    - Absolute tolerance for equality check with zero
C     Z      - Solution(s) of second-order equation with one unknown
C
C PRINT SWITCH:
C     .LT. 0 - No print
C     .EQ. 0 - Print error messages only
C     .GE. 5 - Print ENTER and LEAVE
C     .GE. 6 - Print in addition INPUT parameters
C     .GE. 7 - Print in addition OUTPUT parameters
C
C SUBROUTINES CALLED:
C     DSOLV3
C
C CODING:
C     18-SEP-96,  Knut Morten Okstad,   NTNU
C
C REVISION:
C
C********************** END SUBROUTINE DESCRIPTION *********************
C
      IMPLICIT NONE
C
C                                      * Global variables
      INTEGER            IPSW,IWR,NSOL
      DOUBLE PRECISION   EPS,A(4),B(4),X(4),Y(4)
C
C                                      * Internal variables
      INTEGER            I,J,NX,NY
      DOUBLE PRECISION   DET,Q0,Q1,Q2,TOL,Z(2)
C
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D0 )
C
C**************************** END VARIABLES ****************************
C
C        Entry section
C
      IF (IPSW .GE. 5)                           WRITE(IWR,9900)
      IF (IPSW .GE. 6)                           WRITE(IWR,9920) EPS,A,B
C
C
C        Program logic section
C
C                                      * Initiate geometric tolerance
      NSOL = 0
      TOL  = ZERO
      DO 1000 I = 1,4
         TOL = MAX(TOL,EPS*ABS(A(I)),EPS*ABS(B(I)))
 1000 CONTINUE
C
      DET = A(2)*B(3) - B(2)*A(3)
      IF (ABS(A(1)) .LT. TOL .AND. ABS(B(1)) .LT. TOL) THEN
C
C                                      * Linear set of equations
C        | A2 A3 | * | X | = | A4 |
C        | B2 B3 |   | Y | = | B4 |
C
         IF (ABS(DET) .LT. TOL*TOL)              GOTO 8000
C
         NSOL = 1
         X(1) = ( B(3)*A(4) - A(3)*B(4)) / DET
         Y(1) = (-B(2)*A(4) + A(2)*B(4)) / DET
C
      ELSE
C
C                                      * Solve the second order equation
C                                           Q2*x^2 + Q1*x + Q0 = 0
         Q0 = B(3)*A(4) - A(3)*B(4)
         Q1 = B(1)*A(4) - A(1)*B(4) - DET
         Q2 = A(1)*B(2) - B(1)*A(2)
         CALL DSOLV3 (IPSW,IWR,EPS,ZERO,Q2,Q1,Q0,NX,Z)
C
C                                      * Find y = (A4-A2*x) / (A1*x+A3)
         DO 2000 I = 1,NX
            Q0 = A(1)*Z(I) + A(3)
            IF (ABS(Q0) .LE. TOL)                GOTO 2000
            NSOL    = NSOL + 1
            X(NSOL) = Z(I)
            Y(NSOL) = (A(4)-A(2)*Z(I)) / Q0
 2000    CONTINUE
C
C                                      * Solve the second order equation
C                                           Q2*y^2 + Q1*y + Q0 = 0
         Q0 = B(2)*A(4) - A(2)*B(4)
         Q1 = B(1)*A(4) - A(1)*B(4) + DET
         Q2 = A(1)*B(3) - B(1)*A(3)
         CALL DSOLV3 (IPSW,IWR,EPS,ZERO,Q2,Q1,Q0,NY,Z)
C
C                                      * Find x = (A4-A3*y) / (A1*y+A2)
         NX = NSOL
         DO 3100 I = 1,NY
            Q0 = A(1)*Z(I) + A(2)
            IF (ABS(Q0) .LE. TOL)                GOTO 3100
C
C                                      * Avoid solutions equal to those
C                                        found by solving for x first
            DO 3000 J = 1,NX
               Q1 = MAX(ABS(Y(J)),ABS(Z(I)))
               IF (ABS(Y(J)-Z(I)) .LE. EPS*Q1)   GOTO 3100
 3000       CONTINUE
C
            NSOL    = NSOL + 1
            Y(NSOL) = Z(I)
            X(NSOL) = (A(4)-A(3)*Z(I)) / Q0
 3100    CONTINUE
C
      ENDIF
C
C
C        Closing section
C
 8000 CONTINUE
      IF (IPSW .GE. 5)                           WRITE(IWR,9910)
      IF (NSOL .LT. 1 .AND. IPSW .GE. 0)         THEN
         WRITE(IWR,9970)
         WRITE(IWR,9920) EPS,A,B
      ELSE IF (IPSW .GE. 7)                      THEN
         WRITE(IWR,9930) NSOL,(I,X(I),Y(I),I=1,NSOL)
      ENDIF
C
      RETURN
C
 9900 FORMAT( / ' ENTERING SUBROUTINE DSLBLN' )
 9910 FORMAT( / ' LEAVING SUBROUTINE DSLBLN' )
 9920 FORMAT(   '     WITH INPUT PARAMETERS:'
     +        / '     EPS    = ', 1P,E13.5,
     +        / '     A      = ',   4E13.5,
     +        / '     B      = ',   4E13.5 )
 9930 FORMAT(   '     WITH OUTPUT PARAMETERS:'
     +        / '     NSOL   = ', I6,
     +        /('     X(',I1,')   = ', 1P,2E13.5) )
 9970 FORMAT( / ' *** DSLBLN FOUND NO SOLUTIONS' )
C
      END
