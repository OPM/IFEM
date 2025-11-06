      SUBROUTINE Q4_INV (NDIM,EPSO,EPSZ,X,TENC,XI,IPSW,IWR,IERR)
C
C $Id$
C***************************** USERS' PART *****************************
C
C A F E M  ROUTINE : Q4_INV
C
C PURPOSE:
C     Q4_INV : Find natural coordinates for a point within a Q4-element.
C
C METHOD:
C     The natural coordinates {XI} corresponding to the
C     point {X} are calculated. The element nodes are assumed to have
C     natural coordinates (-1,-1), (1,-1), (1,1), (-1,1), respectively.
C
C INPUT ARGUMENTS:
C     NDIM   - Number of spatial dimensions
C     EPSO   - Tolerence for checking if point is outside the element
C     EPSZ   - Tolerance for equality check with zero
C     X      - Cartesian coordinates of a point within the element
C     TENC   - Table of element nodal coordinates
C     IPSW   - Print switch
C     IWR    - Print unit number
C
C INPUT/OUTPUT ARGUMENTS:
C     None
C
C OUTPUT ARGUMENTS:
C     XI     - Natural coordinates of the point {X}
C     IERR   - Error flag
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
C     .GE. 7 - Print ENTER and LEAVE
C     .GE. 8 - Print in addition INPUT parameters
C     .GE. 9 - Print in addition OUTPUT parameters
C
C SUBROUTINES CALLED:
C     DDOT   (double precision function)
C     DMUAX3
C     DSLBLN
C     SHAPE2
C     SHLTRI
C
C CODING:
C     19-SEP-96,  Knut Morten Okstad,   NTNU
C
C REVISION:
C
C********************** END SUBROUTINE DESCRIPTION *********************
C
      IMPLICIT NONE
C
C                                      * Global variables
C
      INTEGER            NDIM,IPSW,IWR,IERR
      DOUBLE PRECISION   EPSO,EPSZ,X(NDIM),TENC(NDIM,4),XI(2)
C
C                                      * Internal variables
      LOGICAL            LXY
      INTEGER            i,IERRL,j,NSOL
      DOUBLE PRECISION   A(4,2),D1,D2,TOL,TGL(9),TLX(2,5),
     +                   X1(4),X2(4),XL(3)
C
      DOUBLE PRECISION   ZERO        , C         , ONE
      PARAMETER        ( ZERO = 0.0D0, C = 0.25D0, ONE = 1.0D0 )
C
      DOUBLE PRECISION   DDOT
C
C**************************** END VARIABLES ****************************
C
C        Entry section
C
      IF (IPSW .GE. 7)                           WRITE(IWR,9900)
      IF (IPSW .GE. 8)                           THEN
         WRITE(IWR,9920) EPSO,EPSZ,X
         CALL DPRINT (TENC,NDIM,4,'TENC  ',IWR)
      ENDIF
C
C
C        Program logic section
C
C                                      * Check input arguments
      IERR = -1
      IF (EPSZ .LE. ZERO .OR. EPSO .LT. EPSZ)    GOTO 7200
C
      IF (NDIM .EQ. 2)                           THEN
         LXY = .TRUE.
      ELSE IF (NDIM .EQ. 3)                      THEN
         TOL = MAX(TENC(3,1),TENC(3,2),TENC(3,3),TENC(3,4))
     +       - MIN(TENC(3,1),TENC(3,2),TENC(3,3),TENC(3,4))
         LXY = TOL .LT. EPSZ
      ELSE
                                                 GOTO 7200
      ENDIF
C
      IF (LXY)                                   THEN
C
C                                      * The element is parallel
C                                        to the global XY-plane
         DO 1000 j = 1,2
            A(1,j) = ( TENC(j,1) - TENC(j,2) + TENC(j,3) - TENC(j,4))*C
            A(2,j) = (-TENC(j,1) + TENC(j,2) + TENC(j,3) - TENC(j,4))*C
            A(3,j) = (-TENC(j,1) - TENC(j,2) + TENC(j,3) + TENC(j,4))*C
            A(4,j) = ( TENC(j,1) + TENC(j,2) + TENC(j,3) + TENC(j,4))*C
 1000    CONTINUE
C
         A(4,1) = X(1) - A(4,1)
         A(4,2) = X(2) - A(4,2)
C
      ELSE
C
C                                      * The element is not parallel
C                                        to the global XY-plane, need
C                                        to calculate local coordinates
         DO 2000 j = 1,NDIM
            X1(j)  = TENC(j,3) - TENC(j,1)
            X2(j)  = TENC(j,4) - TENC(j,2)
 2000    CONTINUE
         CALL SHLTRI (2,X1,X2,TGL,IPSW,IWR,IERRL)
         IF (IERRL .LT. 0)                       GOTO 7100
C
         DO 2100 j = 1,4
            CALL DMUAX3 (TGL,TENC(1,j),TLX(1,j))
 2100    CONTINUE
         CALL DMUAX3 (TGL,X,XL)
C
C                                      * Initiate the aux. variables
         DO 2200 j = 1,2
            A(1,j) = ( TLX(j,1) - TLX(j,2) + TLX(j,3) - TLX(j,4)) * C
            A(2,j) = (-TLX(j,1) + TLX(j,2) + TLX(j,3) - TLX(j,4)) * C
            A(3,j) = (-TLX(j,1) - TLX(j,2) + TLX(j,3) + TLX(j,4)) * C
            A(4,j) = ( TLX(j,1) + TLX(j,2) + TLX(j,3) + TLX(j,4)) * C
 2200    CONTINUE
C
         A(4,1) = XL(1) - A(4,1)
         A(4,2) = XL(2) - A(4,2)
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C     We have to solve the following set of equations (j=1,2):
C
C        A1j*XI*ETA + A2j*XI + A3j*ETA = A4j
C
C     The way that we may solve this nonlinear (in XI and ETA)
C     set of equations depends on the coefficients Aij.
C     The solution is unique for proper input. See DSLBLN for details.
C
C-----------------------------------------------------------------------
C
      CALL DSLBLN (IPSW,IWR,EPSZ,A(1,1),A(1,2),NSOL,X1,X2)
      IF (NSOL .LE. 0)                           GOTO 7100
C
C                                      * Check that the solutions
C                                        are "inside" the element
      TOL = ONE + EPSO
      DO 3000 j = 1,NSOL
         IF (ABS(X1(j)).LT.TOL .AND. ABS(X2(j)).LT.TOL) THEN
            IERR = IERR + 1
            IF (IERR .GT. 0)                     THEN
C
C                                      * We have multiple solutions.
C                                        Check if they are almost equal.
C
               IF (ABS(XI(1)-X1(j)) .LT. EPSZ*(TOL+TOL) .AND.
     +             ABS(XI(2)-X2(j)) .LT. EPSZ*(TOL+TOL)) THEN
                  XI(1) = 0.5D0*(XI(1)+X1(j))
                  XI(2) = 0.5D0*(XI(2)+X2(j))
                  IERR  = IERR - 1
                                                 GOTO 3000
               ENDIF
C
C                                      * Choose the one closest to X
C
               CALL SHAPE2 (0,4,XI,TGL,IPSW,IWR,IERRL)
               IF (IERRL .LT. 0)                 GOTO 7100
C
               D1 = ZERO
               DO 2500 i = 1,NDIM
                  D1 = D1 + (X(i)-DDOT(4,TGL,1,TENC(i,1),NDIM))**2.0D0
 2500          CONTINUE
C
               XL(1) = X1(j)
               XL(2) = X2(j)
               CALL SHAPE2 (0,4,XL,TGL,IPSW,IWR,IERRL)
               IF (IERRL .LT. 0)                 GOTO 7100
C
               D2 = ZERO
               DO 2600 i = 1,NDIM
                  D2 = D2 + (X(i)-DDOT(4,TGL,1,TENC(i,1),NDIM))**2.0D0
 2600          CONTINUE
               IF (D1 .LT. D2)                   GOTO 3000
C
            ENDIF
            XI(1) = X1(j)
            XI(2) = X2(j)
         ENDIF
 3000 CONTINUE

      IF (IERR .GT. 0)                           GOTO 6200
      IF (IERR .EQ. 0)                           GOTO 8000
C
C
C        Warning section
C
 6100 CONTINUE
      IERR = 1
      IF (IPSW .GE. 0)                           THEN
         WRITE(IWR,9961) (j,X1(j),X2(j),j=1,NSOL)
         WRITE(IWR,9920) EPSO,EPSZ,X
         CALL DPRINT (TENC,NDIM,4,'TENC  ',IWR)
      ENDIF
                                                 GOTO 8000
 6200 CONTINUE
      IERR = 2
      IF (IPSW .GE. 0)                           THEN
         WRITE(IWR,9962) (j,X1(j),X2(j),j=1,NSOL)
         WRITE(IWR,9963) XI,MIN(D1,D2)
         WRITE(IWR,9920) EPSO,EPSZ,X
         CALL DPRINT (TENC,NDIM,4,'TENC  ',IWR)
      ENDIF
                                                 GOTO 8000
C
C        Error section
C
 7100 CONTINUE
      IERR = -1
                                                 GOTO 7900
 7200 CONTINUE
      IERR = -2
C
 7900 CONTINUE
      IF (IPSW .GE. 0)                           THEN
         CALL ERRMSG ('Q4_INV',IWR,IERR)
      ENDIF
C
C
C        Closing section
C
 8000 CONTINUE
      IF (IPSW .GE. 7)                           WRITE(IWR,9910)
      IF (IPSW .GE. 9)                           WRITE(IWR,9930) IERR,XI
C
      RETURN
C
 9900 FORMAT( / ' ENTERING SUBROUTINE Q4_INV' )
 9910 FORMAT( / ' LEAVING SUBROUTINE Q4_INV' )
 9920 FORMAT(   '     WITH INPUT PARAMETERS:'
     +        / '     EPS    = ', 1P,2E13.5,
     +        / '     X      = ',    3E13.5 )
 9930 FORMAT(   '     WITH OUTPUT PARAMETERS:'
     +        / '     IERR   = ', I6,
     +        / '     XI,ETA = ', 1P,2E13.5 )
 9961 FORMAT( / '  ** WARNING FROM SUBROUTINE Q4_INV'
     +        / '     THE POINT(S) ARE OUTSIDE THE ELEMENT',
     +        /('     XI(',I1,')  = ', 1P,2E13.5) )
 9962 FORMAT( / '  ** WARNING FROM SUBROUTINE Q4_INV'
     +        / '     THE SOLUTION IS NOT UNIQUE'
     +        /('     XI(',I1,')  = ', 1P,2E13.5) )
 9963 FORMAT(   '     CHOSE THE SOLUTION CLOSEST TO {X}'
     +        / '     XI,ETA = ', 1P,2E13.5, '     DELTA  = ', E12.5 )
C
      END
