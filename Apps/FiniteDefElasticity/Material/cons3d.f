      subroutine cons3d (ipsw, iter, iwr, lfirst, mTYP, mVER, nHV, nTM,
     &                   detF, F, pMAT, HV, Engy, Sig, Cst, ierr)
C
C ----------------------------------------------------------------------
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
C     ipsw     - Print switch
C     iter     - Iteration number
C     iwr      - Write unit number
C     lFIRST   - Control variable
C              = 0 : The number of computed equilibrium configurations
C                    are greater or equal to one
C              = 1 : No equilibrium configuration has been computed
C     mTYP     - Material type number
C     mVER     - Material version number
C     nHV      - Number of history parameters/stress points
C     nTM      - Number og stress/strain components per point
C     detF     - Determinant of the deformation gradient
C     F        - Deformation gradient
C     pMAT     - Material property array
C     HV       - History parameters at t_n
C
C ARGUMENTS OUTPUT:
C     HV       - History parameters at t_n+1
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
C                = .true. if incompressible material ): Pratio = 0.5
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
      integer   mTYP, mVER, ipsw, iter, iwr, lfirst, nHV, nTM, ierr
      real*8    detF, Engy, F(3,*), pMAT(*), Sig(*), Cst(*), HV(*)
C
      real*8    Bmod, Emod, Pratio, Smod
      logical   incomp
C
      include 'include/feninc/const.h'
C
C         Entry section
C
      ierr = 0
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'ENTERING SUBROUTINE CONS3D'
          if (ipsw .gt. 2)                        then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9020) 'iter   =', iter
              write(iwr,9020) 'lfirst =', lfirst
              write(iwr,9020) 'mTYP   =', mTYP
              write(iwr,9020) 'mVER   =', mVER
              write(iwr,9020) 'nHV    =', nHV
              write(iwr,9020) 'nTM    =', nTM
              write(iwr,9030) 'detF   =', detF
              write(iwr,9030) 'pMAT   =', pMAT(1:2)
              call rprin0(F   , 3,    3, 'F     ', iwr)
              call rprin0(HV  , 1,  nHV, 'HV    ', iwr)
          endif
      endif
C
C
C         Compute Shear and Bulk modulus
C
      if (mTYP .lt. 0) then
         Bmod   = pmat(1)
         Smod   = pmat(2)
         incomp = Bmod .eq. zero
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
      if (abs(mTYP) .eq. 10)                      then
C
         if (mVER .eq. 0)                         then
C
C         Linear elastic isotropic material
C
            call liel3d(ipsw, iwr, detF, F, Bmod, Smod, Engy, Sig, Cst,
     &                  ierr)
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
            call stnh3d(ipsw, iwr, mVER, detF, F, Bmod, Smod, Engy, Sig,
     &                  Cst, ierr)
            if (ierr .lt. 0)                      go to 7000
C
         elseif (mVER .gt. 10 .and. mVER .lt. 15) then
C
C         Modified hyperelastic neoHookean model
C
            call monh3d(ipsw, iwr, mVER-10, detF, F, Bmod, Smod, Engy,
     &                  Sig, Cst, ierr)
            if (ierr .lt. 0)                      go to 7000
C
         elseif (mVER .eq. 20)                    then
C
C         Isotropic elasticity formulated in terms of
C         principal logarithmic stretches
C
            call ieps3d(ipsw, iwr, detF, F, Bmod, Smod, Engy, Sig, Cst,
     &                  ierr)
            if (ierr .lt. 0)                      go to 7000
C
         else
C
C         Illegal material version
C
                                                  go to 7100
         endif
C
      elseif (abs(mTYP) .eq. 40)                  then
C
C         Finite stretch plasticity model
C
         call plas3d(ipsw, iwr, iter, lfirst, ntm, pMAT, detF, F(1,1),
     &               F(1,4), HV(1), HV(7), HV(8), Sig, Cst, ierr)
         if (ierr .lt. 0)                         go to 7000
C
         Engy = zero
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
      write(iwr,9020) 'iter   =', iter
      write(iwr,9020) 'lfirst =', lfirst
      write(iwr,9020) 'mTYP   =', mTYP
      write(iwr,9020) 'mVER   =', mVER
      write(iwr,9020) 'nHV    =', nHV
      write(iwr,9020) 'nTM    =', nTM
      write(iwr,9030) 'detF   =', detF
      write(iwr,9030) 'pMAT   =', pMAT(1:2)
      call rprin0(F   , 3,    3, 'F     ', iwr)
      call rprin0(HV  , 1,  nHV, 'HV    ', iwr)
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
      write(iwr,9030) 'Engy   =', Engy
      call rprin0(HV  , 1,  nHV, 'HV    ', iwr)
      call rprin0(Sig , 1,    6, 'Sig   ', iwr)
      call rprin0(Cst , 6,    6, 'Cst   ', iwr)
      call flush(iwr)
                                                  go to 8010
C
C         Closing section
C
 8000 continue
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'LEAVING SUBROUTINE CONS3D'
          if (ipsw .gt. 3)                        then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              write(iwr,9030) 'Engy   =', Engy
              call rprin0(HV  , 1,  nHV, 'HV    ', iwr)
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
