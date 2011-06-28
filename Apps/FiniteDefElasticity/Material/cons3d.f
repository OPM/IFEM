      subroutine cons3d (ipsw, iter, iwr, lfirst, mTYP, mVER, nHV, nSig,
     &                   nTM, detF, F, pMAT, HV, Engy, Sig, Cst, ierr)
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
C     nSig     - Number of stress components per point
C     nTM      - Number og strain components per point
C     detF     - Determinant of the deformation gradient
C     F        - Deformation gradient
C     pMAT(13) - Material property array:
C                pMAT(1) = Emod : Elastic Bulk  modulus
C                pMAT(2) = Poissons ratio
C                pMAT(5) = Bmod : Bulk modulus
C                pMAT(6) = Smod : Shear modulus
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
      integer   mTYP, mVER, ipsw, iter, iVOL, iwr, lfirst, nHV, nSig,
     &          nTM, ierr
      real*8    detF, Engy, F(3,*), pMAT(*), Sig(*), Cst(*), HV(*)
C
      real*8    Bmod, Emod, Pratio, Smod
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
              write(iwr,9020) 'nSig   =', nSig
              write(iwr,9020) 'nTM    =', nTM
              write(iwr,9030) 'detF   =', detF
              write(iwr,9030) 'pMAT   =', pMAT(1:2)
              call rprin0(F   , 3,    3, 'F     ', iwr)
              call rprin0(HV  , 1,  nHV, 'HV    ', iwr)
          endif
      endif
C
C         Compute Shear and Bulk modulus
C
      if (mTYP .lt. 0) then
         Bmod   = pMAT(1)
         Smod   = pMAT(2)
      else
         Bmod   = pMAT(5)
         Smod   = pMAT(6)
      end if
C
      iVOL = nint(pMAT(7))
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
         elseif (mVER .eq. 1)                     then
C
C         Standard hyperelastic neoHookean model
C
C         Compute Lame's parameter for compressibel material
C
C            if (.not. incomp)  Bmod = Bmod - two3 * Smod
C
            call stnh3d(ipsw, iwr, iVOL, detF, F, Bmod, Smod, Engy, Sig,
     &                  Cst, ierr)
            if (ierr .lt. 0)                      go to 7000
C
         elseif (mVER .eq. 2)                     then
C
C         Modified hyperelastic neoHookean model
C
            call monh3d(ipsw, iwr, iVOL, detF, F, Bmod, Smod, Engy,
     &                  Sig, Cst, ierr)
            if (iERR .lt. 0)                      go to 7000
C
         elseif (mVER .eq. 3)                     then
C
C         Modified Ogden hyperelastic material formulated in terms of
C         principal logarithmic stretches
C
            call ieps3d(ipsw, iwr, 1, iVOL, detF, F, Bmod, Smod, 
     &                  pmat(8), Engy, Sig, Cst, ierr)
            if (iERR .lt. 0)                      go to 7000
C
         elseif (mVER .eq. 4)                     then
C
C         Isotropic elasticity formulated in terms of
C         principal logarithmic stretches
C
            call ieps3d(ipsw, iwr, 2, iVOL, detF, F, Bmod, Smod, 
     &                  pmat(8), Engy, Sig, Cst, ierr)
            if (iERR .lt. 0)                      go to 7000
C
         else
C
C         Illegal material version
C
                                                  go to 7100
         endif
C
      elseif (abs(mTYP) .eq. 40 .and. nHV .ge. 10) then
C
C         Finite stretch plasticity model
C
         call plas3d(ipsw, iwr, iter, lfirst, nSig, nTM, pMAT, detF,
     &               F(1,1), F(1,4), HV(1), HV(7), HV(8), Sig, Cst,
     &               iERR)
         if (iERR .lt. 0)                         go to 7000
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
      write(iwr,9020) 'nSig   =', nSig
      write(iwr,9020) 'nTM    =', nTM
      write(iwr,9030) 'detF   =', detF
      write(iwr,9030) 'pMAT   =', pMAT(1:2)
      call rprin0(F   , 3,    3, 'F     ', iwr)
      call rprin0(HV  , 1,  nHV, 'HV    ', iwr)
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
      write(iwr,9030) 'Engy   =', Engy
      call rprin0(HV  , 1,  nHV, 'HV    ', iwr)
      call rprin0(Sig , 1, nSig, 'Sig   ', iwr)
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
              call rprin0(Sig , 1, nSig, 'Sig   ', iwr)
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
