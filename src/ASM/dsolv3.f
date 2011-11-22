      SUBROUTINE DSOLV3 (IPSW,IWR,EPS,A,B,C,D,NSOL,X)
C
C $Id$
C***************************** USERS' PART *****************************
C
C A F E M  ROUTINE : DSOLV3
C
C PURPOSE:
C     DSOLV3 : Solve the cubic equation  A*x^3 + B*x^2 + C*x + D = 0.
C
C METHOD:
C     See K. Rottmann, "Matematische Formelsammlung" (1960) pp. 13-16.
C
C INPUT ARGUMENTS:
C     IPSW   - Print switch
C     IWR    - Print unit
C     EPS    - Tolerance for equality check with zero
C     A - D  - Coefficients of the cubic equation
C
C INPUT/OUTPUT ARGUMENTS:
C     None
C
C OUTPUT ARGUMENTS:
C     NSOL   - Number of solutions found
C     X      - Array containing the solutions
C
C RESTRICTIONS:
C     None
C
C************************** MAINTENANCE PART ***************************
C
C INTERNAL VARIABLES:
C     None
C
C PRINT SWITCH:
C     .LT. 0 - No print
C     .EQ. 0 - Print error messages only
C     .GE. 5 - Print ENTER and LEAVE
C     .GE. 6 - Print in addition INPUT parameters
C     .GE. 7 - Print in addition OUTPUT parameters
C
C SUBROUTINES CALLED:
C     None
C
C CODING:
C     18-MAY-96,  Knut Morten Okstad,   NTNU
C
C REVISION:
C
C********************* END SUBROUTINE DESCRIPTION **********************
C
      IMPLICIT NONE
C
C                                      * Global variables
      INTEGER            IPSW,IWR,NSOL
      DOUBLE PRECISION   EPS,A,B,C,D,X(3)
C
C                                      * Internal variables
      INTEGER            I
      DOUBLE PRECISION   ATOL,FI,PI,P,Q,W,Y1,Y2,Y3
C
C**************************** END VARIABLES ****************************
C
C        Entry section
C
      IF (IPSW .GE. 5)                           WRITE(IWR,9900)
      IF (IPSW .GE. 6)                           THEN
         WRITE(IWR,9920) EPS,A,B,C,D
      ENDIF
C
C
C        Program logic section
C
C                                      * Initiate
      NSOL = 0
      ATOL = EPS * MAX(ABS(A),ABS(B),ABS(C),ABS(D))
      IF (ATOL .LT. EPS*EPS)                     GOTO 8000
C
      IF (ABS(A) .GT. ATOL)                      THEN
C
C                                      * Cubic equation
         NSOL = 3
         P =  (C - B*B/(3.0D0*A)) / (3.0D0*A)
         Q = ((2.0D0*B*B/(27.0D0*A) - C/3.0D0) * B/A + D) / (A+A)
         W = Q*Q + P*P*P
C
         IF (W .LE. -EPS)                        THEN
            PI = ACOS(-1.0D0)
            FI = ACOS(-Q/SQRT(-P*P*P))
            Y1 =  2.0D0*SQRT(-P) * COS(FI/3.0D0)
            Y2 = -2.0D0*SQRT(-P) * COS((FI+PI)/3.0D0)
            Y3 = -2.0D0*SQRT(-P) * COS((FI-PI)/3.0D0)
         ELSE IF (W .LT. EPS)                    THEN
            Y1 = 2.0D0 * (-Q)**(1.0D0/3.0D0)
            Y2 = -Y1 / 2.0D0
            Y3 =  Y2
         ELSE
            NSOL = -3
         ENDIF
C
         IF (NSOL .EQ. 3)                        THEN
            W    = B/(3.0D0*A)
            X(1) = Y1 - W
            X(2) = Y2 - W
            X(3) = Y3 - W
         ENDIF
C
      ELSE IF (ABS(B) .GT. ATOL)                 THEN
C
C                                      * Quadratic equation
         NSOL = 2
         P = C*C - 4.0D0*B*D
         IF (P .GE. ATOL*ATOL)                   THEN
            Q    = SQRT(P)
            X(1) = (-C+Q)/(B+B)
            X(2) = (-C-Q)/(B+B)
         ELSE IF (P .GT. -ATOL*ATOL)             THEN
            X(1) = -C/(B+B)
            X(2) = X(1)
         ELSE
            NSOL = -2
         ENDIF
C
      ELSE IF (ABS(C) .GT. ATOL)                 THEN
C
C                                      * Linear equation
         NSOL = 1
         X(1) = -D/C
C
      ENDIF
C
C
C        Closing section
C
 8000 CONTINUE
      IF (IPSW .GE. 5)                           WRITE(IWR,9910)
      IF (NSOL .EQ. 0 .AND. IPSW .GE. 0)         THEN
         WRITE(IWR,9970)
         WRITE(IWR,9920) EPS,A,B,C,D
      ELSE IF (NSOL .LT. 0 .AND. IPSW .GE. 0)    THEN
         WRITE(IWR,9971) ABS(NSOL)
         WRITE(IWR,9920) EPS,A,B,C,D
      ELSE IF (IPSW .GE. 7)                      THEN
         WRITE(IWR,9930) NSOL,(I,X(I),I=1,NSOL)
      ENDIF
C
      RETURN
C
 9900 FORMAT( / ' ENTERING SUBROUTINE DSOLV3' )
 9910 FORMAT( / ' LEAVING SUBROUTINE DSOLV3' )
 9920 FORMAT(   '     WITH INPUT PARAMETERS:'
     +        / '     EPS    = ', 1P,E13.5
     +        / '     ABCD   = ',   4E13.5 )
 9930 FORMAT(   '     WITH OUTPUT PARAMETERS:'
     +        / '     NSOL   = ', I6,
     +        /('     X(',I1,')   = ', 1P,E13.5) )
 9970 FORMAT( / ' *** DSOLV3 FOUND NO SOLUTIONS' )
 9971 FORMAT( / ' *** DSOLV3 DETECTED',I2,' COMPLEX ROOTS' )
C
      END
