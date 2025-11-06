      SUBROUTINE DPRINT (A,M,N,CHEAD,IWR)
      INTEGER              M,N,      IWR
      DOUBLE PRECISION   A(M,*)
      CHARACTER                CHEAD*(*)
C
C PURPOSE:
C     DPRINT : Print the double precision array A(M,N).
C
C METHOD:
C     The array is printed as a two-dimensional matrix.
C     Rows and columns are numbered and each element is printed with
C     the format 1PE14.6, i.e. 7 significant digits.
C     A heading (e.g. the array name) is also printed.
C
C INPUT ARGUMENTS:
C     A      - Array to be printed
C     M      - Number of rows
C     N      - Number of columns
C     CHEAD  - Heading
C     IWR    - Print unit number
C
      INTEGER     NCOL,J,I,K,L
      PARAMETER ( NCOL = 9 )
C
      IF (N .LE. 0 .OR. M .LE. 0)                GOTO 8000
C
      WRITE(IWR,9000) CHEAD
C
      DO 2000 J = 1,N,NCOL
         K = MIN(N,J+NCOL-1)
C
         WRITE(IWR,9100) (L,L=J,K)
         WRITE(IWR,9200) ('--------------',L=J,K)
C
         DO 1000 I = 1,M
            WRITE(IWR,9300) I,(A(I,L),L=J,K)
 1000    CONTINUE
C
 2000 CONTINUE
C
 8000 CONTINUE
      RETURN
C
 9000 FORMAT( / 5X, A, ' = ' )
 9100 FORMAT( / 4X, 9(I9,5X) )
 9200 FORMAT(   4X, 9A14 )
 9300 FORMAT( I3, '|', 1P,9E14.6 )
C
      END
