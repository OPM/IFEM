      subroutine cons3d (mTYP, mVER, detF, F, pmat, Engy, Sig, Cst,
     &                   ipswb3, iwr, ierr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   CONS3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     CONS3D: Driver for finite deformation constitutive models for
C             3D solids
C
C METHOD:
C     Compute Cauchy stresses and spatial constitutive tensor
C     for finite deformation constitutive models for 3D solids.
C
C     The following material models are available:
C     mTYP = 10  Linear elastic material, 3-D stress
C     mTYP = 20  Hyperelastic compressible neoHookean model, 3-D stress
C     mTYP = 40  Elastic-plastic strain hardening material,
C                combined isotropic and kinematic hardening,
C                3-D stress
C
C ARGUMENTS INPUT:
C     mTYP    - Material type number
C     mVER    - Material version number
C     detF    - Determinant of the deformation gradient
C     F       - Deformation gradient
C     pmat(2) - Material property array
C
C ARGUMENTS OUTPUT:
C     Engy     - Strain energy density
C     Sig(6)   - Cauchy stresses (spatial stresses)
C     Cst(6,6) - Spatial constitutive tensor
C     ierr     - Error flag
C               = -1 : Illegal number of element nodes
C
C COMMON BLOCKS:
C     const
C
C PERIPHERAL UNITS:
C     None
C
C INTERNAL VARIABLES:
C     Bmod     - Material property (Bulk modulus)
C     Emod     - Material property (Youngs modulus)
C     Pratio   - Material property (Poissons ratio)
C     Smod     - Material property (Shear modulus)
C     incomp   - logical variable
C                = .true. if incompressible material ): Pratio = 0
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
      integer   mTYP, mVER, ipswb3, iwr, ierr
      real*8    detF, Engy, F(9), pmat(2), Sig(6), Cst(6,6)
      real*8    Bmod, Emod, Pratio, Smod
      logical   incomp
C
      include  'const.h'
C
C         Entry section
C
      ierr = 0
C
      if (ipswb3 .gt. 0)                          then
          write(iwr,9010) 'ENTERING SUBROUTINE CONS3D'
          if (ipswb3 .gt. 2)                      then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9020) 'mTYP   =', mTYP
              write(iwr,9020) 'mVER   =', mVER
              write(iwr,9030) 'detF   =', detF
              write(iwr,9030) 'Pmat   =', pmat
              call rprin0(F   , 3,    3, 'F     ', iwr)
          endif
      endif
C
C
C         Compute Shear and Bulk modulus
C
      if (mTYP .lt. 0) then
         Bmod = pmat(1)
         Smod = pmat(2)
      else
         Emod   = pmat(1)
         Pratio = pmat(2)
         Smod   = one2 * Emod / (one + Pratio)
         incomp = Pratio .eq. one2
         if (incomp) then
C
C         Incompressible material
C
            Bmod = zero
         else
            Bmod = one3 * Emod / (one - two * Pratio)
         end if
      end if
C
C         Compute Cauchy stresses and spatial constitutive tensor
C
      if (abs(mTYP) .eq. 10) then
C
         if (mVER .eq. 0) then
C
C         Linear elastic isotropic material
C
            call liel3d(detF, F, Bmod, Smod, Engy, Sig, Cst,
     &                  ipswb3, iwr, ierr)
            if (ierr .lt. 0)                      go to 7000
C
         elseif (mVER .gt. 0 .and. mVER .lt. 5)   then
C
C         Standard hyperelastic neoHookean model
C
C         Compute Lame's parameter for compressibel material
C
            if (.not. incomp)  Bmod = Bmod - two3 * Smod
C
            call stnh3d(mVER, detF, F, Bmod, Smod, Engy, Sig, Cst,
     &                  ipswb3, iwr, ierr)
            if (ierr .lt. 0)                      go to 7000
C
         elseif (mVER .gt. 10 .and. mVER .lt. 15) then
C
C         Modified hyperelastic neoHookean model
C
            call monh3d(mVER-10, detF, F, Bmod, Smod, Engy, Sig, Cst,
     &                  ipswb3, iwr, ierr)
            if (ierr .lt. 0)                      go to 7000
C
         elseif (mVER .eq. 20)                    then
C
C         Isotropic elasticity formulated in terms of
C         principal logarithmic stretches
C
            call ieps3d(detF, F, Bmod, Smod, Engy, Sig, Cst,
     &                  ipswb3, iwr, ierr)
            if (ierr .lt. 0)                      go to 7000
C
         else
C
C         Illegal material version
C
                                                  go to 7100
         endif
C
      else
C
C         Illegal material type
C
                                                  go to 7200
      endif
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
 7200 continue
      ierr =-3
                                                  go to 7900
C     ----------
 7900 continue
      write(iwr,9020) '*** ERROR IN SUBROUTINE CONS3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      write(iwr,9020) 'mTYP   =', mTYP
      write(iwr,9020) 'mVER   =', mVER
      write(iwr,9030) 'detF   =', detF
      write(iwr,9030) 'Pmat   =', pmat
      call rprin0(F   , 3,    3, 'F     ', iwr)
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
      write(iwr,9030) 'Engy   =', Engy
      call rprin0(Sig , 1,    6, 'Sig   ', iwr)
      call rprin0(Cst , 6,    6, 'Cst   ', iwr)
      call flush(iwr)
                                                  go to 8010
C
C         Closing section
C
 8000 continue
C
      if (ipswb3 .gt. 0)                          then
          write(iwr,9010) 'LEAVING SUBROUTINE CONS3D'
          if (ipswb3 .gt. 3)                      then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              write(iwr,9030) 'Engy   =', Engy
              call rprin0(Sig , 1,    6, 'Sig   ', iwr)
              call rprin0(Cst , 6,    6, 'Cst   ', iwr)
          endif
          call flush(iwr)
      endif
C
 8010 continue
C
      return
C
 9010 format(3X,A)
 9020 format(3X,A,I6)
 9030 format(3X,A,1P,2E15.6)
C
      end
