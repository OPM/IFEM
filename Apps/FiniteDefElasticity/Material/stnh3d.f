      subroutine stnh3d (iVF, detF, F, lambda, mu, Engy, Sig, Cst,
     &                   ipswb3, iwr, ierr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   STNH3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     MONH3D: Compute Cauchy stresses and spatial constitutive tensor
C             for a standard NeoHookean isotropic 3D material
C
C METHOD:
C     Energy function: Standard Neohookean model
C
C                 W(J,b) = lambda*U(J) + 0.5*mu*(I_c - 3 - 2*ln(J))
C
C                 Sig = [mu*(b - 1) + lambda*Press]/detF 
C                     =  mu*(b - 1)/detF + lambda*partial_U/partial_J
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
C     iVF     - Volumetric function type
C     detF    - Determinant of the deformation gradient
C     F       - Deformation gradients 
C     lambda  - Lame's constant (=K-2*G/3)
C     mu      - Shear modulus
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
C     b(6)    - Left Cauchy-Green deformation tensor
C     
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
      integer   iVF, ipswb3, iwr, ierr
      integer   i, i3, ierrl
C
      real*8    detF, Engy, F(3,3), lambda, mu, Sig(6), Cst(6,6)
      real*8    c1, c2, b(6), Cstv, detFi, Press
C
      include 'const.h'
C
C         Entry section
C
      ierr = 0
C
      if (ipswb3 .gt. 0)                          then 
          write(iwr,9010) 'ENTERING SUBROUTINE STNH3D'
          if (ipswb3 .gt. 2)                      then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9020) 'iVF    =', iVF
              write(iwr,9030) 'detF   =', detF
              write(iwr,9030) 'lambda =', lambda
              write(iwr,9030) 'mu     =', mu    
              call rprin0(F   , 3,    3, 'F     ', iwr)
          endif
      endif
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
C         Compute volumetric stresses (pressure) and material moduli
C         according to the volumetric function type
C
      detFi = one / detF
C
      call vfnh3d (iVF, detF, detFi, lambda, Engy, Press, Cstv,
     &             ipswb3, iwr, ierrl)
      if (ierrl .lt. 0)                     go to 7000
C
C         Initialize spatial constitutive tensor
C
      Cst = zero
C     
C         Accumulate deviatoric and volumetric contribution to
C         spatial stresses and material moduli
C
      c1    = mu * detFi
C      c2    = two * mu
      c2    = two * c1
      Cstv  = Cstv * detF
C
C        Diagonal terms
C
      do i          = 1, 3
         i3         = i + 3
         Sig(i)     = c1 * b(i) - c1 + Press
         Sig(i3)    = c1 * b(i3)
         Cst(i ,i ) = c2 - Press + Cstv
         Cst(i3,i3) = c1 - Press
      end do ! i
C
C       Off-diagonal terms
C
      c1 = Press + Cstv
C
      Cst(1,2) = c1
      Cst(1,3) = c1
      Cst(2,1) = c1     
      Cst(2,3) = c1
      Cst(3,1) = c1
      Cst(3,2) = c1
C
C       Compute strain energy desnity
C
      Engy = Engy + mu*(one2*(b(1)+b(2)+b(3)) - onep5 - log(abs(detF)))
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
      write(iwr,9020) '*** ERROR IN SUBROUTINE STNH3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      write(iwr,9020) 'iVF    =', iVF
      write(iwr,9030) 'detF   =', detF
      write(iwr,9030) 'lambda =', lambda
      write(iwr,9030) 'mu     =', mu   
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
          write(iwr,9010) 'LEAVING SUBROUTINE STNH3D'
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
