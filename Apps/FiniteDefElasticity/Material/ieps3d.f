      subroutine ieps3d (detF, F, lambda, mu, Engy, Sig, Cst,
     &                   ipswb3, iwr, ierr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   IEPS3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     IEPS3D: Compute Cauchy stresses and spatial constitutive tensor
C             for a isotropic elastic 3D material formulated in
C             principal stretches
C
C METHOD:
C     
C     
C
C ARGUMENTS INPUT:  
C     detF     - Determinant of the deformation gradient
C     F(3,3)   - Deformation gradient 
C     lambda   - Lame's constant (=K-2*G/3)
C     mu       - Shear modulus
C
C ARGUMENTS OUTPUT:
C     Engy     - Strain energy density
C     Sig(6)   - Cauchy stresses (spatial stresses)
C     Cst(6,6) - Spatial constitutive tensor
C     ierr     - Error flag
C               = 1 : 
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
C             November 14, 2010
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
      integer   i, j, k, l, iVF, i3, ierrl, nrotd
C
      real*8    detF, Engy, F(3,3), lambda, mu, Sig(6), Cst(6,6)
      real*8    aap(6,6), b(6), bpr(6), Cstv, detFi, EngU, EngW,
     &          epsd(3), Press, q(6,6), qen(3,3), taup(3)
C
      include 'include/feninc/const.h'
C
C         Entry section
C
      ierr = 0
C
      if (ipswb3 .gt. 0)                          then 
          write(iwr,9010) 'ENTERING SUBROUTINE IEPS3D'
          if (ipswb3 .gt. 2)                      then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9030) 'detF   =', detF
              write(iwr,9030) 'lambda =', lambda
              write(iwr,9030) 'mu     =', mu    
              call rprin0(F , 3, 3, 'F     ', iwr)
          endif
      endif
C
C         Initialize Cauchy stresses and spatial constitutive tensor
C
      Sig = zero
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
C         Compute principal stretches and directions
C
C         Initialization: Upon entry to EIGS3D array qen(3,3) contains 
C                         left Cauchy-Green tensor b; upon exit, 
C                         qen(3,3) contains eigenvectors
C
      do i  = 1,3
C
         i3 = i + 3
         j  = 1 + mod(i,3)
C
         qen(i,i) = b(i)
         qen(i,j) = b(i3)
         qen(j,i) = b(i3)
C
      end do ! i
C
      call eigs3d(qen, bpr, nrotd, ipswb3, iwr, ierr)
      if (ierrl .lt. 0)                           go to 7000
C
C         Compute transformation matrix from
C         Principal to Standard Basis
C
      do i = 1,3
C
         k = 1 + mod(i,3)
C
         do j = 1,3
C
            l = 1 + mod(j,3)
C
            q(i,j)     = qen(i,j)*qen(i,j)
C
            q(i+3,j)   = qen(i,j)*qen(k,j)
C
            q(j,i+3)   = qen(j,i)*qen(j,k)*two
C
            q(j+3,i+3) = qen(j,i)*qen(l,k) + qen(j,k)*qen(l,i)
C
         end do ! j
C
      end do ! i
C
C         Compute Deviatoric Strains, deviatoric Kirchhoff 
C         Stresses and Material Moduli in Principal Basis
C 
      call wlog3d (mu, bpr, epsd, taup, aap, EngW, ipswb3, iwr, ierrl)
      if (ierrl .lt. 0)                           go to 7000
C
C         Transform deviatoric Kirchhoff Stresses and Material Moduli
C         to Cauchy stresses and spatial material moduli
C 
      call amat3d (taup, aap, q, Sig, Cst, ipswb3, iwr, ierrl)
      if (ierrl .lt. 0)                           go to 7000
C
C         Transform deviatoric Kirchhoff Stresses and Material Moduli
C         to Cauchy stresses and spatial material moduli
C     
      detFi = one / detF
      Sig   = Sig * detFi
      Cst   = Cst * detFi
C
C         Compute volumetric stresses (pressure) and material moduli
C
      iVF   = 1
C
      call vfnh3d (iVF, detF, detFi, lambda, EngU, Press, Cstv,
     &             ipswb3, iwr, ierrl)
      if (ierrl .lt. 0)                           go to 7000
C     
C         Accumulate deviatoric and volumetric contribution to
C         spatial stresses and material moduli
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
C       Compute strain energy desnity
C
      Engy = EngU + EngW
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
      write(iwr,9020) '*** ERROR IN SUBROUTINE IEPS3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      write(iwr,9030) 'detF   =', detF
      write(iwr,9030) 'lambda =', lambda
      write(iwr,9030) 'mu     =', mu    
      call rprin0(F , 3, 3, 'F     ', iwr)
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
          write(iwr,9010) 'LEAVING SUBROUTINE IEPS3D'
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
