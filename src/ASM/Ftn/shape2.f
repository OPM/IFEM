      SUBROUTINE SHAPE2 (IOP,NENOD,XI,S,IPSW,IWR,IERR)
      IMPLICIT NONE
      INTEGER            IOP,NENOD,IPSW,IWR,IERR
      DOUBLE PRECISION   XI(2),S(NENOD)
C
      DOUBLE PRECISION   UMX1,UMX2,UPX1,UPX2
C
      DOUBLE PRECISION   QUART         , ONE
      PARAMETER        ( QUART = 0.25D0, ONE = 1.0D0 )
C
      IF (IOP .EQ. 0 .AND. NENOD .EQ. 3) THEN
         IERR = 0
         S(1) = XI(1)
         S(2) = XI(2)
         S(3) = ONE - XI(1) - XI(2)
      ELSE IF (IOP .EQ. 0 .AND. NENOD .EQ. 4) THEN
         IERR = 0
         UMX1 = ONE - XI(1)
         UMX2 = ONE - XI(2)
         UPX1 = ONE + XI(1)
         UPX2 = ONE + XI(2)
         S(1) = QUART * UMX1*UMX2
         S(2) = QUART * UPX1*UMX2
         S(3) = QUART * UPX1*UPX2
         S(4) = QUART * UMX1*UPX2
      ELSE
         IERR = -2
         IF (IPSW .GE. 0) THEN
            write(IWR,*) '*** SHAPE2:',IOP,NENOD,'not supported'
         ENDIF
      ENDIF
C
      RETURN
      END
