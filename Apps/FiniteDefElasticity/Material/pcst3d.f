      subroutine pcst3d (aa,dd,ipswb3,iwr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   PCST3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     PCST3D: Project 6x6 constitutive tensor (AA) onto a 7x7 equivalent
C             constitutive tensor for mixed model (DD)
C
C METHOD:
C                     | D_11   D_12 |
C                DD = |             |
C                     | D_21   D_22 |
C              where:
C                D_11  = 6 x 6 Deviatoric part of matrix
C                D_12  = 6 x 1 Coupling   part of matrix
C                D_21  = 1 x 6 Coupling   part of matrix
C                D_22  = 1 x 1 Volumetric part of matrix
C
C
C ARGUMENTS INPUT:
C     aa(6,6) - Material tangent matrix (based on F)
C
C ARGUMENTS OUTPUT:
C     dd(7,7) - Mixed material tangent for stiffness computations
C
C COMMON BLOCKS:
C     None
C
C PERIPHERAL UNITS:
C     None
C
C INTERNAL VARIABLES:
C     None
C
C PRINT SWITCH:
C     ipswb3 = 0  Gives no print
C     ipswb3 = 2  Gives enter and leave
C     ipswb3 = 3  Gives in addition parameters on input
C     ipswb3 = 5  Gives in addition parameters on output
C
C LIMITS:
C
C CODING:
C     STANDART ADHERED TO: Computas report no 78-949
C     AUTHOR: Kjell Magne Mathisen
C             NTNU, Department of structural engineering
C             September 14, 2010
C
C INSPECTED:
C
C REVISIONS:
C     xx-OCT-2010  K.M.Mathisen NTNU
C        Changed XG1-XG8 to XG and UG1-UG8 to UG.
C
C ---------------------------------------------------------------------
C
      implicit  none
C
      integer   ipswb3, iwr, j
      real*8    aa(6,6), dd(7,7)
C
      include 'const.h'
C
C
C         Entry section
C
      if (ipswb3 .gt. 0)                          then
         write(iwr,9010) 'ENTERING SUBROUTINE PCST3D'
         if (ipswb3 .gt. 2)                       then
            write(iwr,9010) 'WITH INPUT ARGUMENTS'
            call rprin0(aa, 6, 6, 'aa    ', iwr)
         endif
      endif
C
C
C         Load moduli from constitution
C
      dd(1:6,1:6) = aa
C
C     Compute left and right multiples with trace
C
      dd(1:6,7) = (aa(:,1) + aa(:,2) + aa(:,3))*one3
      dd(7,1:6) = (aa(1,:) + aa(2,:) + aa(3,:))*one3
C
C     Convert upper 6 x 6 to a deviatoric D_11
C
      do j = 1, 3
         dd(1:6,j) = dd(1:6,j) - dd(1:6,7)
         dd(j,1:6) = dd(j,1:6) - dd(7,1:6)
      end do
C
C     Form last term, D_22
C
      dd(7,7) = (dd(1,7) + dd(2,7) + dd(3,7))*one3
C
C     Final update to form D_12 and D_21
C
      dd(1:3,  7) = dd(1:3,  7) - dd(7,7)
      dd(  7,1:3) = dd(  7,1:3) - dd(7,7)
      dd(1:3,1:3) = dd(1:3,1:3) + dd(7,7)
C
C
C         Closing section
C
 8000 continue
      if (ipswb3 .gt. 0)                          then
         write(iwr,9010) 'LEAVING SUBROUTINE PCST3D'
         if (ipswb3 .gt. 3)                       then
            write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
            call rprin0(dd, 7, 7, 'dd    ', iwr)
         endif
         call flush(iwr)
      endif
C
      return
C
 9010 format(3X,A)
C
      end
