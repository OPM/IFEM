      SUBROUTINE RPRIN0 (A,MA,NA,CHHEAD,IW)
      INTEGER              MA,NA,       IW
      DOUBLE PRECISION   A(MA,* )
      CHARACTER                  CHHEAD*(*)
C
C @(#)rprin0.f 1.3 90/11/19  Copyright (c) 1989 by A.S Veritas Research.
C-----------------------------------------------------------------------
C
C F E N R I S   ROUTINE :     RPRIN0                      NTH/SINTEF/DNV
C
C PURPOSE:
C     RPRIN0: Print the real matrix A(MA,NA).
C
C METHOD:
C     THE MATRIX IS PRINTED AS A TWO-DIMENSIONAL MATRIX.ROWS AND COLUMNS
C     ARE NUMBERED,AND EACH ELEMENT IS PRINTED WITH THE FORMAT 1PE14.6 ,
C     I.E. 7 SIGNIFICANTS DIGITS.
C     A HEADING ( E.G. THE SYMBOL NAME "A" OF THE MATRIX ) CONSISTING OF
C     UP TO 6 CHARACTERS IS PRINTED.
C     EXAMPLE OF A CALL TO RPRIN0:  CALL RPRIN0(A,MA,NA,'A     ',IW)
C
C ARGUMENTS INPUT:
C     A      - MATRIX TO BE PRINTED
C     MA     - NUMBER OF ROWS
C     NA     - NUMBER OF COLUMNS
C     CHHEAD - HEADING, E.G. "NAME" OF MATRIX
C     IW     - LINE PRINTER UNIT NUMBER
C
C ARGUMENTS OUTPUT:
C     NONE
C
C COMMON BLOCKS:
C     NONE
C
C PERIPHERAL UNITS:
C     NONE
C
C LIMITATIONS:
C     NONE
C
C CODING:
C     STANDARD ADHERED TO: COMPUTAS REPORT NO 78-949
C     CODED: DECEMBER 5, 1974     K.BELL   SINTEF - DIV. OF STRUCT. ENG.
C
C REVISIONS:
C     DATE: AUGUST 27,1980    ALF ENGSETH  SINTEF - DIV. OF STRUCT. ENG.
C     DATE: DECEMBER 20,1982    MAGNE NYGAARD  VERITAS
C     12-MAR-87  M.K.Nygaard  VERITEC      86.1-0
C        Changed format 9100
C      2-DEC-87  K.M.Mathisen NTH          5.2-04
C        Included test on the value of NA and changed dimensioning of A.
C     16-NOV-88  Ole Stroem / Robert O. Bjaerum, Veritas Research
C        Any dimension of CHHEAD OK. Please keep it under 120.
C
C-----------------------------------------------------------------------
C
      INTEGER     NCOL,J,I,K,L
      PARAMETER ( NCOL = 9 )
C
C EXECUTABLE SECTION
C
C         PROGRAM LOGIC SECTION
C
      IF (NA .LE. 0 .OR. MA .LE. 0)               GO TO 8000
C
C         PRINT HEADING
C
      WRITE(IW,9050) CHHEAD
C
C         PRINT MATRIX
C
      DO 4000 J=1,NA,NCOL
      K = J+(NCOL-1)
      IF ( K .GT. NA ) THEN
        K = NA
      END IF
C
      WRITE(IW,9000) (L,L=J,K)
C
      DO 3000 I=1,MA
      WRITE(IW,9100) I,(A(I,L),L=J,K)
 3000 CONTINUE
 4000 CONTINUE
C
 8000 CONTINUE
C
      RETURN
C
 9000 FORMAT ( / 3X , 9( I9 , 5X ) )
 9050 FORMAT ( / 5X , A , ' = ' )
 9100 FORMAT ( 1X , I3 , 1P , 9E14.6 )
C
      END
