      SUBROUTINE T3_INV (NDIM,EPSO,EPSZ,X,TENC,XI,IPSW,IWR,IERR)
C
C $Id$
C***************************** USERS' PART *****************************
C
C A F E M  ROUTINE : T3_INV
C
C PURPOSE:
C     T3_INV : Find natural coordinates for a point within a triangle.
C
C METHOD:
C     The natural (area) coordinates {XI}={A1,A2} corresponding to the
C     point {X} are calculated. The element nodes are assumed to have
C     the natural coordinates (1,0), (0,1) and (0,0), respectively.
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
C     AREA   - Element area
C     VN     - Element normal vector
C     V12    - Vector from node 1 to node 2
C     V13    - Vector from node 1 to node 3
C     V1X    - Vector from node 1 to the point X
C
C PRINT SWITCH:
C     .LT. 0 - No print
C     .EQ. 0 - Print error messages only
C     .GE. 7 - Print ENTER and LEAVE
C     .GE. 8 - Print in addition INPUT parameters
C     .GE. 9 - Print in addition OUTPUT parameters
C
C SUBROUTINES CALLED:
C     DCROS3
C     DNORM3
C
C CODING:
C      4-NOV-94,  Knut Morten Okstad,   NTH
C
C REVISION:
C      2-OCT-96,  Knut Morten Okstad,   NTNU
C        The routine is included in the AFEM package.
C
C********************** END SUBROUTINE DESCRIPTION *********************
C
      IMPLICIT NONE
C
C                                      * Global variables
C
      INTEGER            NDIM,IPSW,IWR,IERR
      DOUBLE PRECISION   EPSO,EPSZ,X(NDIM),TENC(NDIM,3),XI(2)
C
C                                      * Internal variables
C
      DOUBLE PRECISION   AUX(3),A2,A3,AREA,TOL,
     +                   V12(3),V13(3),V1X(3),VN(3)
C
      DOUBLE PRECISION   ZERO        , ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
C
      DOUBLE PRECISION   DNORM3
C
C**************************** END VARIABLES ****************************
C
C        Entry section
C
      IF (IPSW .GE. 7)                           WRITE(IWR,9900)
      IF (IPSW .GE. 8)                           THEN
         WRITE(IWR,9920) NDIM,EPSO,EPSZ,X
         CALL DPRINT (TENC,NDIM,3,'TENC  ',IWR)
      ENDIF
C
C
C        Program logic section
C
C                                      * Check input arguments
      IERR = 0
      IF (EPSZ .LE. ZERO .OR. EPSO .LT. EPSZ)    GOTO 7200
C
      IF (NDIM .EQ. 2)                           THEN
C
C                                      * Calculate the element area
C
         AREA = ( (TENC(1,1)-TENC(1,3))*(TENC(2,2)-TENC(2,3)) -
     +            (TENC(1,2)-TENC(1,3))*(TENC(2,1)-TENC(2,3)) )
         IF (AREA .LE. EPSZ)                     GOTO 7200
C
C                                      * Calculate the area coordinates
C
         XI(1) = ( (   X(1)  -TENC(1,3))*(TENC(2,2)-TENC(2,3)) -
     +             (TENC(1,2)-TENC(1,3))*(   X(2)  -TENC(2,3)) ) / AREA
         XI(2) = ( (TENC(1,1)-TENC(1,3))*(   X(2)  -TENC(2,3)) -
     +             (   X(1)  -TENC(1,3))*(TENC(2,1)-TENC(2,3)) ) / AREA
C
      ELSE IF (NDIM .EQ. 3)                      THEN
C
C                                      * Calculate the element normal
C                                        vector and the element area
         V12(1) = TENC(1,2) - TENC(1,1)
         V12(2) = TENC(2,2) - TENC(2,1)
         V12(3) = TENC(3,2) - TENC(3,1)
         V13(1) = TENC(1,3) - TENC(1,1)
         V13(2) = TENC(2,3) - TENC(2,1)
         V13(3) = TENC(3,3) - TENC(3,1)
         CALL DCROS3 (V12,V13,VN)
         AREA   = DNORM3(VN)
         IF (AREA .LE. EPSZ)                     GOTO 7200
C
C                                      * Check the distance from point X
C                                        to the element plane
         V1X(1) = X(1) - TENC(1,1)
         V1X(2) = X(2) - TENC(2,1)
         V1X(3) = X(3) - TENC(3,1)
         A3     = VN(1)*V1X(1) + VN(2)*V1X(2) + VN(3)*V1X(3)
         IF (ABS(A3) .GT. SQRT(AREA)*EPSO)       GOTO 7300
C
C                                      * Calculate the area coordinates
         V1X(1) = V1X(1) - A3*VN(1)
         V1X(2) = V1X(2) - A3*VN(2)
         V1X(3) = V1X(3) - A3*VN(3)
C
         CALL DCROS3 (V12,V1X,AUX)
         A3 = SIGN(SQRT(AUX(1)*AUX(1)+AUX(2)*AUX(2)+AUX(3)*AUX(3)),
     +                  AUX(1)*VN(1) +AUX(2)*VN(2) +AUX(3)*VN(3))
         CALL DCROS3 (V1X,V13,AUX)
         A2 = SIGN(SQRT(AUX(1)*AUX(1)+AUX(2)*AUX(2)+AUX(3)*AUX(3)),
     +                  AUX(1)*VN(1) +AUX(2)*VN(2) +AUX(3)*VN(3))
C
         XI(1) = (AREA-A2-A3) / AREA
         XI(2) = A2 / AREA
C
      ELSE
                                                 GOTO 7200
      ENDIF
C
C
C                                      * Check that the coordinates are
C                                        "inside" the element
      TOL = ONE + EPSO
      IF (XI(1)       .LT. -EPSO .OR. XI(1)       .GT. TOL) GOTO 6100
      IF (      XI(2) .LT. -EPSO .OR.       XI(2) .GT. TOL) GOTO 6100
      IF (XI(1)+XI(2) .LT. -EPSO .OR. XI(1)+XI(2) .GT. TOL) GOTO 6100
C
                                                 GOTO 8000
C
C        Warning section
C
 6100 CONTINUE
      IERR = 1
      IF (IPSW .GE. 0)                           THEN
         WRITE(IWR,9961) XI,X
         CALL DPRINT (TENC,NDIM,3,'TENC  ',IWR)
      ENDIF
                                                 GOTO 8000
C
C        Error section
C
 7200 CONTINUE
      IERR = -2
                                                 GOTO 7900
 7300 CONTINUE
      IERR = -3
      IF (IPSW .GE. 0)                           THEN
         WRITE(IWR,9973) AREA
         CALL DPRINT (TENC,NDIM,3,'TENC  ',IWR)
      ENDIF
C
 7900 CONTINUE
      IF (IPSW .GE. 0)                           THEN
         CALL ERRMSG ('T3_INV',IWR,IERR)
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
 9900 FORMAT( / ' ENTERING SUBROUTINE T3_INV' )
 9910 FORMAT( / ' LEAVING SUBROUTINE T3_INV' )
 9920 FORMAT(   '     WITH INPUT PARAMETERS:'
     +        / '     NDIM   = ', I6,
     +        / '     EPS    = ', 1P,2E13.5,
     +        / '     X      = ',    3E13.5 )
 9930 FORMAT(   '     WITH OUTPUT PARAMETERS:'
     +        / '     IERR   = ', I6,
     +        / '     XI,ETA = ', 1P,2E13.5 )
 9961 FORMAT( / '  ** WARNING FROM SUBROUTINE T3_INV'
     +        / '     THE POINT IS OUTSIDE THE ELEMENT',
     +        / '     XI,ETA = ', 1P,2E13.5,
     +        / '     X      = ',    3E13.5 )
 9973 FORMAT( / ' *** THE ELEMENT IS DEGENERATED',
     +        / '     AREA   = ', 1P,E13.5 )
C
      END
