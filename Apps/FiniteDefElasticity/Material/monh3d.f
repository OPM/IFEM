      subroutine monh3d (ipsw, iwr, iVF, detF, F, Bmod, Smod, Engy,
     &                   Sig, Cst, ierr)
C 
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   MONH3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     MONH3D: Compute Cauchy stresses and spatial constitutive tensor
C             for a modified NeoHookean isotropic 3D material
C
C METHOD:
C     Energy function: Compressible Neohookean model
C                      with J_2/3 regularization
C                   _ __      _            _       __
C                 W(J,be) = U(J) + 0.5*mu*(J^(-2/3)be:1 - 3)
C
C     Volumetric function type: U_1 = K*[ 0.25*( J^2 - 1 )-0.5*ln(J)]
C                               U_2 = K*[ 0.5*( J - 1 )^2 ]
C                               U_3 = K*[ 0.5*( ln(J) )^2 ]
C                               U_4 = K*[ 2.0*(J - 1 - ln(J)) ]
C     Free energy function:
C
C     Psi   = U(J) -3*alpha*J*U'(J)*[T - T_ref]
C
C     up    = partial_J ( Psi )
C     upp   = [ partial_J ( J partial_J Psi ) - up ]/J
C           =   partial^2_J ( Psi )/JC
C     
C
C ARGUMENTS INPUT:  
C     ipsw    - Print switch
C     iwr     - Write unit number
C     iVF     - Volumetric function type
C     detF    - Determinant of the deformation gradient
C     F       - Deformation gradient
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
C     nSTR     - Number of stress/strain components 
C     b(6)     - Left Cauchy-Green deformation tensor
C     bm(6)    - Modified Left Cauchy-Green deformation tensor
C     
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
      integer   ipsw, iwr, iVF, ierr
      integer   i, j, ierrl
C
      real*8    detF, Engy, Bmod, Smod, F(3,3), Sig(6), Cst(6,6)
      real*8    c1, c2, c3, Cstv, b(6), bm(6), detFi, Press, trbm3
C
      include 'include/feninc/const.h'
C
C         Entry section
C
      ierr = 0
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'ENTERING SUBROUTINE MONH3D'
          if (ipsw .gt. 2)                        then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9020) 'iVF    =', iVF
              write(iwr,9030) 'detF   =', detF
              write(iwr,9030) 'Bmod   =', Bmod
              write(iwr,9030) 'Smod   =', Smod
              call rprin0(F   , 3,    3, 'F     ', iwr)
          endif
      endif
C
C         Initialize spatial constitutive tensor
C
      Cst = zero
C
C         Compute Left Cauchy-Green deformation tensor
C   
C         b = F * F^T
C
      b(1) = F(1,1)*F(1,1) + F(1,2)*F(1,2) + F(1,3)*F(1,3)
      b(2) = F(2,1)*F(2,1) + F(2,2)*F(2,2) + F(2,3)*F(2,3)
      b(3) = F(3,1)*F(3,1) + F(3,2)*F(3,2) + F(3,3)*F(3,3)
      b(4) = F(1,1)*F(2,1) + F(1,2)*F(2,2) + F(1,3)*F(2,3)
      b(5) = F(2,1)*F(3,1) + F(2,2)*F(3,2) + F(2,3)*F(3,3)
      b(6) = F(1,1)*F(3,1) + F(1,2)*F(3,2) + F(1,3)*F(3,3)
C
C         Compute modified Left Cauchy-Green deformation tensor
C               _   
C         bm =  b = J^(-2/3) * b
C
      detFi = one / detF
      bm    = b * (detFi ** two3)
C
C         Compute one third of the trace of the modified 
C         Left Cauchy-Green deformation tensor
C                     _   
C         trbm3 = tr( b ) / 3
C     
      trbm3 = (bm(1) + bm(2) + bm(3)) * one3 
      bm(1) = bm(1) - trbm3
      bm(2) = bm(2) - trbm3
      bm(3) = bm(3) - trbm3
C
C         Compute deviatoric part of the constitutive tensor
C
C         Part 1:                     _             _     _
C         Rank one update: -2/3 G * ( b x g +  g x  b ) / J
C
      c1 = two3 * Smod
C
      do i = 1, 6
         c2 = c1 * bm(i)
         do j = 1,3
            Cst(i,j) =  Cst(i,j) - c2
            Cst(j,i) =  Cst(j,i) - c2
         end do ! j
      end do ! i
C
C         Part 2:           __                     _
C         Deviatoric term 2 mu [ I - 1/3 g x g ] / J
C
      c1 = Smod * trbm3
      c2 = c1 + c1
      c3 = c2 * one3
C
      do i = 1,3
         Cst(i  ,i  )  = Cst(i  ,i  ) + c2
         Cst(i+3,i+3)  = Cst(i+3,i+3) + c1
         do j = 1,3
            Cst(i ,j ) = Cst(i ,j )   - c3
         end do ! j
      end do ! i
C
C         Compute deviatoric Kirchhoff stress tensor 
C                   _        _
C         tau_dev = b -  tr( b ) / 3 
C
      Sig = Smod * bm 
C
C         Compute spatial deviatoric (isochoric) 
C         stresses and material moduli
C     
      Sig = Sig * detFi
      Cst = Cst * detFi
C
      if (ipsw .gt. 3)                            then
         write(iwr,9030) 'Press =', Press
         call rprin0(Sig, 1, 6,'Sig   ', iwr)
         call rprin0(Cst, 6, 6,'Cst   ', iwr)
      endif
C
C         Compute volumetric stresses (pressure) and material moduli
C         according to the volumetric function type
C
      call vfnh3d (ipsw, iwr, iVF, detF, detFi, Bmod, Engy, Press, Cstv,
     &             ierrl)
      if (ierrl .lt. 0)                           go to 7000
C
      Sig(1:3) = Sig(1:3) + Press
C
      Cstv     = Cstv * detF
C
      Cst(1,1) = Cst(1,1) - Press + Cstv 
      Cst(1,2) = Cst(1,2) + Press + Cstv
      Cst(1,3) = Cst(1,3) + Press + Cstv

      Cst(2,1) = Cst(2,1) + Press + Cstv
      Cst(2,2) = Cst(2,2) - Press + Cstv
      Cst(2,3) = Cst(2,3) + Press + Cstv

      Cst(3,1) = Cst(3,1) + Press + Cstv
      Cst(3,2) = Cst(3,2) + Press + Cstv
      Cst(3,3) = Cst(3,3) - Press + Cstv

      Cst(4,4) = Cst(4,4) - Press
      Cst(5,5) = Cst(5,5) - Press
      Cst(6,6) = Cst(6,6) - Press
C
C         Compute strain energy density
C
      Engy = Engy + onep5 * Smod * (trbm3 - one)
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
      write(iwr,9020) '*** ERROR IN SUBROUTINE MONH3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      write(iwr,9020) 'iVF    =', iVF
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
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'LEAVING SUBROUTINE MONH3D'
          if (ipsw .gt. 3)                        then
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
