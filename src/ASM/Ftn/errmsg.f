      SUBROUTINE ERRMSG (PRGNAM,IWR,IERR)
      CHARACTER          PRGNAM*(*)
      INTEGER                   IWR,IERR
C
C     ERRMSG : Print an error/warning message.
C
C     PRGNAM - Name of subroutine where the error occurred
C     IERR   - Error flag
C              >  1 : Warning
C              = -1 : Traceback
C              < -1 : Error
C
      IF (IERR .LT. -1)                          THEN
         WRITE(IWR,9901) PRGNAM
      ELSE IF (IERR .GT. 1)                      THEN
         WRITE(IWR,9902) PRGNAM
      ELSE IF (IERR .EQ. -1)                     THEN
         WRITE(IWR,9903) PRGNAM
      ENDIF
C
      IF (ABS(IERR) .LE. 1)                      GOTO 8000
C
      IF (IERR .EQ. -2)                          THEN
         WRITE(IWR,9904) 'INVALID INPUT ARGUMENTS'
      ELSE IF (IERR .EQ. -98)                    THEN
         WRITE(IWR,9904) 'ARRAY IS NOT ALLOCATED'
      ELSE IF (IERR .EQ. -99)                    THEN
         WRITE(IWR,9904) 'UNABLE TO ALLOCATE MEMORY'
      ELSE
         WRITE(IWR,9905) IERR
      ENDIF
C
 8000 CONTINUE
      RETURN
C
 9901 FORMAT( / ' *** ERROR RETURN FROM SUBROUTINE ', A )
 9902 FORMAT( / '  ** WARNING FROM SUBROUTINE ', A )
 9903 FORMAT(   '   * TRACEBACK: ', A )
 9904 FORMAT( 5X, A )
 9905 FORMAT( 5X, 'ERROR FLAG = ', I6 )
C
      END
