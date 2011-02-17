      subroutine liel3d (detF, F, Bmod, Smod, Engy, Sig, Cst,
     &                   ipswb3, iwr, ierr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   LIEL3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     LIEL3D: Compute Cauchy stresses and spatial constitutive tensor
C             for linear elastic isotropic 3D materials
C
C METHOD:
C     Material and spatial stress components are stored as:
C
C                   | Stress(1)  Stress(4)  Stress(6) |
C     Stress(3,3) = | Stress(4)  Stress(2)  Stress(5) |
C                   | Stress(6)  Stress(5)  Stress(3) |
C     
C
C ARGUMENTS INPUT:  
C     detF    - Determinant of the deformation gradient
C     F       - Deformation gradients 
C     Bmod    - Bulk modulus
C     Smod    - Shear modulus
C
C ARGUMENTS OUTPUT:
C     Engy     - Strain energy density
C     Sig(6)   - Cauchy stresses (spatial stresses)
C     Cst(6,6) - Spatial constitutive tensor
C     ierr     - Error flag
C               = 1 : Illegal number of element nodes
C
C COMMON BLOCKS:
C     const
C
C PERIPHERAL UNITS:
C     None
C
C INTERNAL VARIABLES:
C     GLe(6)   - Green-Lagrangse strains
C     PK2(6)   - 2nd Piola-Kirchhoff stresses (material stresses)
C     Cmt(6,6) - Material constitutive tensor
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
      integer   ipswb3, iwr, ierr
      integer   i, j, i3
C
      real*8    detF, Engy, F(9), Bmod, Smod, Sig(6), Cst(6,6)
      real*8    c1, c2
      real*8    GLe(6), PK2(6), Cmt(6,6)
C
      include 'const.h'
C
C         Entry section
C
      ierr = 0
C
      if (ipswb3 .gt. 0)                          then 
          write(iwr,9010) 'ENTERING SUBROUTINE LIEL3D'
          if (ipswb3 .gt. 2)                      then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9030) 'detF   =', detF
              write(iwr,9030) 'Bmod   =', Bmod
              write(iwr,9030) 'Smod   =', Smod     
              call rprin0(F   , 3,    3, 'F     ', iwr)
          endif
      endif
C
      Cmt = zero
C
C         Express the material constitutive tensor 
C         in terms of the Shear and Bulk modulus
C
      c1 = Bmod + four3 * Smod
      c2 = Bmod - two3 * Smod      
C
      do i  = 1, 3
         j  = mod(i,3) + 1
         i3 = i + 3
         Cmt(i ,i ) = c1
         Cmt(i ,j ) = c2
         Cmt(j ,i ) = c2
         Cmt(i3,i3) = Smod
      end do ! i
C
C         Compute Green-Lagrange strains from the deformation gradient
C
C         E = (F^T * F - 1) / 2
C
      GLe(1) = one2 * (F(1)*F(1) + F(2)*F(2) + F(3)*F(3) - one)
      GLe(2) = one2 * (F(4)*F(4) + F(5)*F(5) + F(6)*F(6) - one)
      GLe(3) = one2 * (F(7)*F(7) + F(8)*F(8) + F(9)*F(9) - one)
      GLe(4) =        (F(1)*F(4) + F(2)*F(5) + F(3)*F(6))
      GLe(5) =        (F(4)*F(7) + F(5)*F(8) + F(6)*F(9))
      GLe(6) =        (F(1)*F(7) + F(2)*F(8) + F(3)*F(9))
C
C         Compute 2nd Piola Kirchhoff stresses from product of
C         material constitutive tensor and Green-Lagrange strains
C
C         PK2 = C_mt * GLe
C
      PK2 = matmul(Cmt,GLe)
C
C         Compute strain energy density
C
      Engy = one2 * dot_product(PK2,GLe)
C
C         Push stresses and constitutive tensor forward from
C         material to spatial coordinates in current configuration
C
      call push3D (detF, F, PK2, Cmt, Sig, Cst, ipswb3, iwr)
C
                                                  go to 8000
C
C         Error section
C
 7000 continue
      ierr =-1
                                                  go to 7900
 7100 continue
      ierr =-2
                                                  go to 7900
C     ----------
 7900 continue
      write(iwr,9020) '*** ERROR IN SUBROUTINE LIEL3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      write(iwr,9030) 'detF   =', detF
      write(iwr,9030) 'Bmod   =', Bmod
      write(iwr,9030) 'Smod   =', Smod     
      call rprin0(F   , 3,    3, 'F     ', iwr)
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
      write(iwr,9030) 'Engy   =', Engy
      call rprin0(Sig, 1, 6,'Sig   ', iwr)
      call rprin0(Cst, 6, 6,'Cst   ', iwr)
                                                  go to 8010
C
C         Closing section
C
 8000 continue
C
      if (ipswb3 .gt. 0)                          then
          write(iwr,9010) 'LEAVING SUBROUTINE LIEL3D'
          if (ipswb3 .gt. 3)                      then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              write(iwr,9030) 'Engy   =', Engy
              call rprin0(Sig, 1, 6,'Sig   ', iwr)
              call rprin0(Cst, 6, 6,'Cst   ', iwr)
          endif
      endif
C
 8010 continue
C
      return
C
 9010 format(3X,A)
 9020 format(3X,A,I6)
 9030 format(3X,A,E15.6)
C
      end
