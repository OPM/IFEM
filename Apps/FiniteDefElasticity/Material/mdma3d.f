      subroutine mdma3d (P_bar, P_mix, Sig, DD, ipsw, iwr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   MDMA3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     MDMA3D: Compute finite deformation mixed D arrays
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
C ARGUMENTS INPUT:
C     P_bar   - Constitutive pressure
C     P_mix   - Mixed pressure
C     Sig(6)  - Constitutive Cauchy stresses
C     DD(7,7) - Modified constitutive array
C
C ARGUMENTS OUTPUT:
C     DD(7,7) - Mixed material tangent terms
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
C     ipsw = 0  Gives no print
C     ipsw = 2  Gives enter and leave
C     ipsw = 3  Gives in addition parameters on input
C     ipsw = 5  Gives in addition parameters on output
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
      integer   ipsw, iwr, j
      real*8    DD(7,7), P_bar, P_mix, Sig(6)
      real*8    fac1, Sig_D(6)
C
      include 'include/feninc/const.h'
C
C         Entry section
C
      if (ipsw .gt. 0)                            then
         write(iwr,9010) 'ENTERING SUBROUTINE MDMA3D'
         if (ipsw .gt. 2)                         then
            write(iwr,9010) 'WITH INPUT ARGUMENTS'
            write(iwr,9030) 'P_bar  =', P_bar
            write(iwr,9030) 'P_mix  =', P_mix
            call rprin0(Sig, 1, 6, 'Sig   ', iwr)
            call rprin0(DD , 7, 7, 'DD    ', iwr)
         endif
      endif
C
C
C         Compute deviatoric stress
C
      Sig_D(1:3) = two3 *(Sig(1:3) - P_bar)
      Sig_D(4:6) = two3 * Sig(4:6)
C
C         D_11: B_matrix part
C
      fac1 = P_mix - two3 * P_bar
C
      DD(1:3,1:3) = DD(1:3,1:3) + fac1
C
      do j = 1,6
         DD(1:3,j) = DD(1:3,j) - Sig_D(j)
         DD(j,1:3) = DD(j,1:3) - Sig_D(j)
      end do
C
      fac1 = P_bar - P_mix
C
      do j = 1,3
         DD(j  ,j  ) = DD(j  ,j  ) + fac1 * two
         DD(j+3,j+3) = DD(j+3,j+3) + fac1
      end do ! j
C
C         D_12: Coupling matrix with
C
      DD(7,1:6) = DD(7,1:6) + Sig_D
      DD(1:6,7) = DD(1:6,7) + Sig_D
C
C         D_22: Volumetric part
C
      DD(7,7) = DD(7,7) - one3 * P_bar
C
C
C         Closing section
C
 8000 continue
      if (ipsw .gt. 0)                            then
         write(iwr,9010) 'LEAVING SUBROUTINE MDMA3D'
         if (ipsw .gt. 3)                         then
            write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
            call rprin0(DD , 7, 7, 'DD    ', iwr)
         endif
         call flush(iwr)
      endif
C
      return
C
 9010 format(3X,A)
 9030 format(3X,A,E15.6)
C
      end
