      SUBROUTINE SHLTRI (IOP,VX,V3,TGL,IPSW,IWR,IERR)
C
C $Id$
C***************************** USERS' PART *****************************
C
C A F E M  ROUTINE : SHLTRI
C
C PURPOSE:
C     SHLTRI : Calculate a local triad on a shell surface.
C
C METHOD:
C     The local Z-axis is defined by the given vector V3. For IOP=0, the
C     local X-axis is defined by the projection of the global X-axis
C     onto the shell surface. If this projection is zero, the local
C     Y-axis is instead defined by the projection of the global Y-axis.
C     For IOP=1, the local X-axis is defined by the projection of the
C     given vector VX onto the shell surface. For IOP=2, the local
C     X-axis is defined by the given vector VX and the local Z-axis is
C     defined by the cross product between VX and V3.
C
C INPUT ARGUMENTS:
C     IOP    - Option telling how to construct the triad (see above)
C     VX     - Vector defining the local X-axis (for IOP>0 only)
C     IPSW   - Print switch
C     IWR    - Print unit number
C
C INPUT/OUTPUT ARGUMENTS:
C     V3     - Vector defining the local Z-axis
C
C OUTPUT ARGUMENTS:
C     TGL    - Transformation matrix
C     IERR   - Error flag
C
C RESTRICTIONS:
C     None
C
C************************** MAINTENANCE PART ***************************
C
C INTERNAL VARIABLES:
C     V1, V2 - In-plane base vectors
C
C PRINT SWITCH:
C     .LT. 0 - No print
C     .EQ. 0 - Print error messages only
C     .GE. 6 - Print ENTER and LEAVE
C     .GE. 7 - Print in addition INPUT parameters
C     .GE. 8 - Print in addition OUTPUT parameters
C
C SUBROUTINES CALLED:
C     DCROS3
C     DNORM3
C
C CODING:
C      3-JUN-96,  Knut Morten Okstad,   NTNU
C
C REVISION:
C
C********************** END SUBROUTINE DESCRIPTION *********************
C
      IMPLICIT NONE
C
C                                      * Global variables
C
      INTEGER            IOP,IPSW,IWR,IERR
      DOUBLE PRECISION   VX(3),V3(3),TGL(3,3)
C
C                                      * Internal variables
      INTEGER            I,J
C
      DOUBLE PRECISION   ZERO        , EPS
      PARAMETER        ( ZERO = 0.0D0, EPS = 1.0D-8 )
C
      DOUBLE PRECISION   V1(3),V2(3)
C
      DOUBLE PRECISION   DNORM3
C
C**************************** END VARIABLES ****************************
C
C        Entry section
C
      IF (IPSW .GE. 6)                           WRITE(IWR,9900)
      IF (IPSW .GE. 7)                           THEN
         WRITE(IWR,9920) IOP,V3
         IF (IOP .EQ. 1)                         WRITE(IWR,9921) VX
      ENDIF
C
C
C        Program logic section
C
      i = 1
      IERR = 0
      IF (IOP .NE. 2)                            THEN
         IF (DNORM3(V3) .LE. ZERO)               GOTO 7200
      ENDIF
C
C                                      * Set up a local Cartesian basis
C                                        (V1,V2,V3) where V3 coincide
C                                        with the plane normal
C
      IF (IOP .EQ. 1)                            THEN
C
C                                      * Define V1 by the projection of VX
C                                        onto the plane, then V2 = V3 x V1
         CALL DCROS3 (V3,VX,V2)
         CALL DCROS3 (V2,V3,V1)
C
      ELSE IF (IOP .EQ. 2)                       THEN
C
C                                      * Define V2 as V1 x V3,
C                                        then V3 = V1 x V2
         CALL DCROS3 (V1,V3,V2)
         V3(1) = V2(1)
         V3(2) = V2(2)
         V3(3) = V2(3)
         CALL DCROS3 (V1,V2,V3)
         IF (DNORM3(V3) .LE. ZERO)               GOTO 7200
C
      ELSE IF (ABS(V3(2)) .GT. EPS .OR. ABS(V3(3)) .GT. EPS) THEN
C
C                                      * Define V1 by projecting the
C                                        global X-axis onto the plane,
C                                        then V2 = V3 x V1
C
         V1(1) =  V3(2)*V3(2) + V3(3)*V3(3)
         V1(2) = -V3(1)*V3(2)
         V1(3) = -V3(1)*V3(3)
         CALL DCROS3 (V3,V1,V2)
C
      ELSE
C
C                                      * Define V2 by projecting the
C                                        global Y-axis onto the plane,
C                                        then V1 = V2 x V3
         V2(1) = -V3(2)*V3(1)
         V2(2) =  V3(1)*V3(1) + V3(3)*V3(3)
         V2(3) = -V3(2)*V3(3)
         CALL DCROS3 (V2,V3,V1)
C
      ENDIF
C
C                                      * Normalize the base vectors
      i = 2
      IF (DNORM3(V1) .LE. ZERO)                  GOTO 7200
      IF (DNORM3(V2) .LE. ZERO)                  GOTO 7200
C
C                                      * Establish the global-to-local
C                                        transformation matrix, Tgl
      DO 1000 I = 1,3
         TGL(1,I) = V1(I)
         TGL(2,I) = V2(I)
         TGL(3,I) = V3(I)
 1000 CONTINUE
C
                                                 GOTO 8000
C
C         Error section
C
 7200 CONTINUE
      IERR = -2
      IF (IPSW .GE. 0)                           THEN
         CALL ERRMSG ('SHLTRI',IWR,IERR)
         WRITE(IWR,9920) IOP,V3
         IF (IOP .EQ. 1 .OR. IOP .EQ. 2)         WRITE(IWR,9921) VX
         IF (i .EQ. 2)                           WRITE(IWR,9922) V1,V2
      ENDIF
C
C
C        Closing section
C
 8000 CONTINUE
      IF (IPSW .GE. 6)                           WRITE(IWR,9910)
      IF (IPSW .GE. 8)                           THEN
         WRITE(IWR,9930) IERR,((TGL(i,j),j=1,3),i=1,3)
      ENDIF
C
      RETURN
C
 9900 FORMAT( / ' ENTERING SUBROUTINE SHLTRI' )
 9910 FORMAT( / ' LEAVING SUBROUTINE SHLTRI' )
 9920 FORMAT(   '     WITH INPUT PARAMETERS:'
     +        / '     IOP    = ', I6,
     +        / '     V3     = ', 1P,3E13.5 )
 9921 FORMAT(   '     VX     = ', 1P,3E13.5 )
 9922 FORMAT(   '     V1     = ', 1P,3E13.5
     +        / '     V2     = ',    3E13.5 )
 9930 FORMAT(   '     WITH OUTPUT PARAMETERS:'
     +        / '     IERR   = ', I6,
     +        / '     TGL    = ', 1P,3E13.5, 2(/14X,3E13.5) )
C
      END
