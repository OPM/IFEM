      subroutine plas3d (ipsw, iwr, iter, lfirst, nSig, nTM, pMAT, detF,
     &                   Fn1, Fn, be, Epp, Epl, Sig, Cst, ierr)
C
C ----------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   PLAS3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     PLAS3D: Compute Cauchy stresses and spatial constitutive tensor
C             for Finite Deformation Isotropic I1-J2-J3 Plasticity
C             Models in Principal Logarithmic Stretches
C             Incremental  Lagrangian  Formulation
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
C     ipsw     - Print switch
C     iwr      - Write unit
C     iter     - Iteration number
C     lFIRST   - Control variable
C                = 0 : The number of computed equilibrium configurations
C                      are greater or equal to one
C                = 1 : No equilibrium configuration has been computed
C     nSig     - Number of stress components pr point
C     nTM      - Number of strain components pr point
C     pMAT(13) - Material properties:
C                pMAT(1) = Emod : Elastic Bulk  modulus
C                pMAT(2) = Poissons ratio
C                pMAT(5) = Bmod : Bulk modulus
C                pMAT(6) = Smod : Shear modulus
C                pMAT(7) = H_iso : Isotropic Hardening Modulus [times sqrt(2/3)]
C                pMAT(8) = H_kin : Kinematic Hardening Modulus [times (2/3)]
C                pMAT(9) = iYIELD (Yield Function number)
C                          iYIELD   = 1     : von Mises
C                          pMAT(10) = Y0    : Yield stress (initial)
C                          pMAT(11) = Y_inf : Yield stress (infinity)
C                          pMAT(12) = Delta : Hardening exponent
C                pMAT(13) = istrt : start state
C                           = 0 : Elastic
C                           = 1 : Last solution
C     detF     - Determinant of the deformation gradient at t_n+1
C     Fn1(3,3) - Deformation gradient                    at t_n+1
C     Fn(3,3)  - Deformation gradient                    at t_n
C History variables:
C     be(*)    - Left Cauchy-Green deformation tensor    at t_n
C     Epp      - Cumulative plastic strain               at t_n
C     Epl(*)   - Plastic Strain for Hardening            at t_n
C
C ARGUMENTS OUTPUT:
C History variables:
C     be(*)    - Left Cauchy-Green deformation tensor    at t_n+1
C     Epp      - Cumulative plastic strain               at t_n+1
C     Epl(*)   - Plastic Strain for Hardening            at t_n+1
C     Sig(*)   - Cauchy stresses (spatial stresses)
C     Cst(6,6) - Spatial constitutive tensor
C     ierr     - Error flag
C
C COMMON BLOCKS:
C     const
C
C PERIPHERAL UNITS:
C     None
C
C INTERNAL VARIABLES:
C     vol       - volumetric strain   -   vol = ( eps : 1 ) / 3
C     nn_t(3,3) - principal directions (by columns) tensor form
C     ll2(3)    - squares of principal stretches =  lambda^2
C     tau(3)    - principal values   total    Kirchhoff stress
C     tt(3)     - principal values deviatoric Kirchhoff stress
C     pp        - pressure         volumetric Kirchhoff stress
C     dtde(3,3) - Kirchhoff stress derivative
C                 dtde(a,b) = [ d tau_a / d (lambda_b^2) ] * lambda_b^2
C                           = [ d tau_a / d ( eps_b ) ] / 2.0d0
C     Cstm(6,6) - spatial elastic moduli in principal basis
C     yield     - yield function value
C     yfunct    - yield function
C                 flg = .false.  only yield function
C                 flg = .true.   yield function and derivatives
C     state    -
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
      integer   ipsw, iwr, istrt, iter, lfirst, nSig, nTM, ierr
      real*8    pMAT(*), detF, Fn1(3,3), Fn(3,3),
     &          Epp, be(*), Epl(*), Sig(*), Cst(6,6)
C
      integer   i, j, a, b, rot, it, itmax, p1(6), p2(6), ierrl, iop
      real*8    Bmod, Smod, F(3,3), Finv(3,3), detFr, detr
      real*8    pp, TwoG, K3inv, G2inv, Ypr, YY, dyld
      real*8    ll2_tr(3), ll2(3), tau(3), tt(3)
      real*8    eps_tr(3), th_tr, vol_tr, ee_tr(3), alp_n(3)
      real*8    eps_e(3),         vol_e,  ee_e(3) , alp(3), ta(3)
      real*8    I1, J2, J3, f1, f2, f3, f11, f22, f33, f12, f13, f23
      real*8    gam, Hk, Hkr, aatr, epse, Y0, Yinf, beta
      real*8    err, dsol(7), res(7), tres(7,7), fss(3,3)
      real*8    dtde(3,3), ff(6,6), Cstm(6,6), Eppn
      real*8    nn_t(3,3), nn(3), aa(3), bb(3)
      real*8    tolb, tolc, d_el, xx, f112, f123, yield, yfunct
      real*8    temp
      logical   conv, state
C
      parameter ( itmax = 50, tolb = 1.0d-8, tolc = 1.0d-9 )
C
      include 'include/feninc/const.h'
C
      data     p1    /1,2,3,1,2,3/
      data     p2    /1,2,3,2,3,1/
C
C         Entry section
C
      ierr = 0
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'ENTERING SUBROUTINE PLAS3D'
          if (ipsw .gt. 2)                        then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9020) 'iter   =', iter
              write(iwr,9020) 'lfirst =', lfirst
              write(iwr,9020) 'nSig   =', nSig
              write(iwr,9020) 'nTM    =', nTM
              write(iwr,9030) 'detF   =', detF
              write(iwr,9030) 'Epp    =', Epp
              call rprin0(pMAT, 1,  13, 'pMAT  ', iwr)
              call rprin0(Fn1 , 3,   3, 'Fn1   ', iwr)
              call rprin0(Fn  , 3,   3, 'Fn    ', iwr)
              call rprin0(be  , 1, nTM, 'be    ', iwr)
              call rprin0(Epl , 1,   3, 'Epl   ', iwr)
          endif
      endif
C
C         Initialize Left Cauchy-Green deformation tensor
C
      if (lfirst .eq. 1)                          then
         be(1:3) = one
      end if
C
      Cstm = zero
C
C     Fetch material parameters  &  CONSTANT SETUP
C
      Bmod  = pMAT(5)
      Smod  = pMAT(6)
      Y0    = sqt23 * pMAT(10)  ! N.B. Radius of yield = sqrt(2/3) Sig_y
      Hk    = two3 *  pMAT(8)   ! N.B. Kinematic hardening = 2/3 * H_kin
      istrt = nint(pMAT(13))
C
      if (Hk .gt. zero)                          then
         Hkr = one / Hk
      else
         Hkr = zero
      endif
C
      TwoG  = two * Smod
C
C     Check state for iterations
C
      state = .true.
      dyld  = zero
      if (iter .eq. 0)  then           ! First iteration in step
         if (istrt .eq. 0) then        ! Elastic state requested
            state = .false.
         else                          ! Last state requested
            dyld  = 1.0d-08*Y0
         endif
      endif
C
C     KINEMATIC COMPUTATIONS  (Incremental Lagrangian formulation)
C
C     Inverse of Fn(3,3) * J_n
C
      Finv(1,1) = Fn(2,2)*Fn(3,3) - Fn(2,3)*Fn(3,2)
      Finv(2,1) = Fn(2,3)*Fn(3,1) - Fn(2,1)*Fn(3,3)
      Finv(3,1) = Fn(2,1)*Fn(3,2) - Fn(2,2)*Fn(3,1)
C
      Finv(1,2) = Fn(3,2)*Fn(1,3) - Fn(3,3)*Fn(1,2)
      Finv(2,2) = Fn(3,3)*Fn(1,1) - Fn(3,1)*Fn(1,3)
      Finv(3,2) = Fn(3,1)*Fn(1,2) - Fn(3,2)*Fn(1,1)
C
      Finv(1,3) = Fn(1,2)*Fn(2,3) - Fn(1,3)*Fn(2,2)
      Finv(2,3) = Fn(1,3)*Fn(2,1) - Fn(1,1)*Fn(2,3)
      Finv(3,3) = Fn(1,1)*Fn(2,2) - Fn(1,2)*Fn(2,1)
C
      detFr = one /(Fn(1,1)*Finv(1,1)+
     &              Fn(1,2)*Finv(2,1)+
     &              Fn(1,3)*Finv(3,1))
C
      do j = 1,3
         do i = 1,3
            F(i,j) = (Fn1(i,1)*Finv(1,j)
     &             +  Fn1(i,2)*Finv(2,j)
     &             +  Fn1(i,3)*Finv(3,j))*detFr
         end do ! i
      end do ! j
C
C     Compute Elastic left Cauchy-Green tensor:    be = F * Cp * F^T
C     in current configuration
C
      iop = 1
      call push3d (ipsw, iop, iwr, one, F, be, Cst, be, Cst)
C      call pushr2(F,be,be,one)      ! Update configuration
C
C     Compute principal stretches and directions
C
      nn_t(1,1) = be(1)
      nn_t(2,2) = be(2)
      nn_t(3,3) = be(3)
C
      nn_t(1,2) = be(4)
      nn_t(2,1) = be(4)
C
      nn_t(2,3) = be(5)
      nn_t(3,2) = be(5)
C
      nn_t(1,3) = be(6)
      nn_t(3,1) = be(6)
C
      call eigs3d (ipsw, iwr, nn_t, ll2_tr, rot, ierrl)
      if (ierrl .lt. 0)                           go to 7000
C      call eig3(nn_t, ll2_tr, rot)
C
C     COMPUTE TRIAL KIRCHHOFF STRESS  ( pressure and deviator )
C
      do a = 1, 3
         temp = sqrt(ll2_tr(a))
         ee_tr(a) = log(temp) ! log ( lambda(a)^TR )
         if (ipsw .ge. 5)                         then
            write(iwr,9030) 'll2_tr(a) =', ll2_tr(a)
            write(iwr,9030) 'temp      =', temp
            write(iwr,9030) 'ee_tr(a)  =', ee_tr(a)
         end if
      end do ! a
C
      th_tr = ee_tr(1) + ee_tr(2) + ee_tr(3)
C
      if (ipsw .ge. 5)                            then
         write(iwr,9030) 'th_tr  =', th_tr
         call rprin0(ee_tr, 1, 3, 'ee_tr ',iwr)
      end if
C
      ee_tr = ee_tr - one3*th_tr           ! Trial deviators
C
      pp   = Bmod * th_tr                  ! Pressure: K*th_tr
      Eppn = Epp
C
      tt    = TwoG * ee_tr ! Trial deviatoric stress
      tau   = tt + pp
      alp   = Hk * Epl(1:3)
      alp_n = alp
C
C     Deviatoric: ta =  tt - alp_dev
C
      aatr = (alp(1) + alp(2) + alp(3))*one3
C
      ta = tt - alp + aatr    ! Trial SigMA = ss - alp
C
C     CHECK ELASTIC / PLASTIC STEP
C
C     Compute stress invariant and yield function
C
      I1 =  three * (pp - aatr)
      J2 = ( ta(1)*ta(1) + ta(2)*ta(2) + ta(3)*ta(3) ) * one2
      J3 = ( ta(1)*ta(1)*ta(1) + ta(2)*ta(2)*ta(2) +
     &                           ta(3)*ta(3)*ta(3)   ) * one3
C
      gam   = zero
      yield = yfunct(ipsw,iwr,pMAT,
     &               Epp,I1,J2,J3, f1,f2,f3,
     &               f11,f22,f33,f12,f13,f23,Ypr,YY,gam,.false.)
C
C     Check yield
C
      if ( (yield+dyld .gt. zero) .and. state )   then
C
C     PLASTIC STEP  -->  RETURN MAP
C
         conv   = .false.
         it     =  1
C
         K3inv  = one3 / Bmod
         G2inv  = one2 / Smod
C
         vol_tr = th_tr * one3
         vol_e  = K3inv * pp
C
         eps_tr = ee_tr + vol_tr
         ee_e   = G2inv * tt
         eps_e  = ee_e  + vol_e
C
         epse = 1.d0/(abs(eps_tr(1)) + abs(eps_tr(2)) + abs(eps_tr(3)))
C
         d_el = ( K3inv - G2inv ) * one3
C
         do while ((.not.conv).and.(it.le.itmax))

            yield = yfunct(ipsw,iwr,pMAT,
     &                     Epp,I1,J2,J3, f1,f2,f3,
     &                     f11,f22,f33,f12,f13,f23,Ypr,YY,gam,.true.)
C
            xx  = f1 - f3 * two3 * J2
C
            do a = 1, 3
               nn(a) = xx + ( f2 + f3 * ta(a) ) * ta(a)
            end do ! a
C
            res(1:3) = eps_e - eps_tr + gam * nn
            res(4:6) = (alp_n - alp)*Hkr + gam * nn
            res(7) = yield
C
            err  = (abs(res(1)) + abs(res(2)) + abs(res(3)))*epse
     &           +  abs(res(7))/Y0
            conv = err .lt. tolc
C
C         Construct local tangent matrix
C
            do a     = 1, 3
               aa(a) = f23   * ta(a) + f13
               bb(a) = ta(a) * ta(a) - two3 * J2
            end do ! a
C
            f112 = f11 - one3 * f2
            f123 = f12 - two3 * f3
C
            do a = 1, 3
C
               do b = a, 3
                  fss(a,b) = ( f112 + f22 * ta(a) * ta(b)
     &                     +   f33 * bb(a) * bb(b)
     &                     +  f123 * ( ta(b) + ta(a) )
     &                     + aa(a) * bb(b) + bb(a) * aa(b) ) * gam
                  fss(b,a) = fss(a,b)
               end do ! b
C
               fss(a,a) = fss(a,a) + gam * ( f2 + two * f3 * ta(a) )
C
            end do ! a
C
            do a = 1, 3
C
               tres(a,1:3) = fss(a,:) + d_el
C
               tres(a,a) = tres(a,a) + G2inv
C
            end do ! a
C
C         Kinematic and Isotropic Hardening
C
            if (Hk .gt. zero)                     then
C
               do a = 1,3
C
                  do b = 1,3
                     tres(a,b+3)   = -fss(a,b)
                     tres(a+3,b)   = -fss(a,b)
                     tres(a+3,b+3) =  fss(a,b)
                  end do ! b
C
                  tres(a+3,a+3) = tres(a+3,a+3) + one/Hk
                  tres(a  ,7)   =  nn(a)
                  tres(a+3,7)   = -nn(a)
                  tres(7,a  )   =  nn(a)
                  tres(7,a+3)   = -nn(a)
C
               end do ! a
C
               tres(7,7) = -Ypr
C
               call invert (7, 7, tres, ierrl)
               if (ierrl .lt. 0)                  go to 7000
C
               do i = 1, 7
                 dsol(i) = - tres(i,1)*res(1)
     &                     - tres(i,2)*res(2)
     &                     - tres(i,3)*res(3)
     &                     - tres(i,4)*res(4)
     &                     - tres(i,5)*res(5)
     &                     - tres(i,6)*res(6)
     &                     - tres(i,7)*res(7)
               end do ! i
C
C         Isotropic hardening only
C
            else
C
               do a = 1,3
                  tres(a,4) = nn(a)
                  tres(4,a) = nn(a)
               end do ! a
C
               tres(4,4) = -Ypr
C
               call invert (7, 4, tres, ierrl)
               if (ierrl .lt. 0)                   go to 7000
C
               do i = 1, 4
                  dsol(i) = - tres(i,1)*res(1)
     &                      - tres(i,2)*res(2)
     &                      - tres(i,3)*res(3)
     &                      - tres(i,4)*res(7)
               end do ! i
C
               dsol(7) = dsol(4)
               dsol(4:6) = zero
C
            endif
C
C         Update Kirchhoff stress and plastic flow
C
            if (ipsw .ge. 5)                      then
               call rprin0(tau , 1, 3, 'tau   ',iwr)
               call rprin0(dsol, 1, 3, 'dsol  ',iwr)
            end if
C
            tau = tau + dsol(1:3)
            gam = gam + dsol(7)
C
            if (ipsw .ge. 5)                      then
               call rprin0(tau , 1, 3, 'tau   ',iwr)
            end if
C
C         Update accumulated plastic strain
C
            Epp = Eppn + sqt23*gam
C
C         Update Back Stress
C
            alp = alp + dsol(4:6)
C
C         Update vol.-dev. Kirchhoff stress and stress invariants
C
            pp = ( tau(1) + tau(2) + tau(3) ) * one3
            tt = tau - pp
C
            aatr = (alp(1) + alp(2) + alp(3))*one3
C
            ta = tt - alp + aatr
C
            I1 =  three * (pp - aatr)
            J2 = ( ta(1)*ta(1)       + ta(2)*ta(2)       +
     &                                 ta(3)*ta(3)         ) * one2
            J3 = ( ta(1)*ta(1)*ta(1) + ta(2)*ta(2)*ta(2) +
     &                                 ta(3)*ta(3)*ta(3)   ) * one3
C
C         Update vol.-dev. logarithmic strain
C
            vol_e  = K3inv * pp
C
            do a = 1, 3
               ee_e(a)  = G2inv * tt(a)
               eps_e(a) = ee_e(a) + vol_e
            end do ! a
C
            it = it + 1
C
         end do ! while
C
C       Warning: check convergence
C
         if (.not.conv .and. iter.gt.0)          then
            write(  *,*) ' *WARNING* No convergence in PLASFD',err,tolc
            write(iwr,*) ' *WARNING* No convergence in PLASFD',err,tolc
C          call plstop()
         endif
C
C       Update elastic left Cauchy-Green tensor and plastic acc. strain
C
         do a = 1, 3
            ll2(a) = exp( two * eps_e(a) )
         end do ! a
C
         do i = 1, 6
            be(i) = ll2(1) * nn_t(p1(i),1) * nn_t(p2(i),1)
     &            + ll2(2) * nn_t(p1(i),2) * nn_t(p2(i),2)
     &            + ll2(3) * nn_t(p1(i),3) * nn_t(p2(i),3)
         end do ! i
C
C       Update plastic strains
C
         Epl(1:3) = Epl(1:3) + gam * nn
C
C       Compute elasto-plastic tangent
C
         dtde = one2 * tres(1:3,1:3)
C
      else
C
C     ELASTIC STEP  ( only tangent computation )
C
         dtde(1,1) = one2 * Bmod  + two3 * Smod
         dtde(1,2) = one2 * Bmod  - one3 * Smod
         dtde(1,3) = dtde(1,2)
         dtde(2,1) = dtde(1,2)
         dtde(2,2) = dtde(1,1)
         dtde(2,3) = dtde(1,2)
         dtde(3,1) = dtde(1,2)
         dtde(3,2) = dtde(1,2)
         dtde(3,3) = dtde(1,1)
C
      end if
C
C     COMPUTE CAUCHY STRESS
C
      detr = one / detF
C
      if (ipsw .ge. 5)                            then
         call rprin0(tau , 1, 3, 'tau   ',iwr)
         call rprin0(nn_t, 3, 3, 'nn_t  ',iwr)
      end if
C
      do i = 1, min(6,nSig)
         Sig(i) = (tau(1) * nn_t(p1(i),1)*nn_t(p2(i),1)
     &          +  tau(2) * nn_t(p1(i),2)*nn_t(p2(i),2)
     &          +  tau(3) * nn_t(p1(i),3)*nn_t(p2(i),3)) * detr
      end do ! i
C
      if (nSig .ge. 10) then
         Sig(7:8) = zero
         Sig(  9) = Epp                      ! Plastic strain for output
         Sig( 10) = sqrt(onep5)*(YY + yield) ! Yield stress value
      end if
C
C
C     TANGENT TRANSFORMATION
C
C     Material tangent (computation in the principal basis)
C
      do a = 1, 3
C
C       Upper 3x3 block of Cstm()
C
         Cstm(1:3,a) = two * dtde(1:3,a)
C
         Cstm(a,a) = Cstm(a,a) - two * tau(a)
C
C       Lower 3x3 block of Cstm() [ diagonal block ]
C
         b = mod(a,3) + 1
         if (abs(ll2_tr(a)-ll2_tr(b)).gt.tolb)    then
            Cstm(a+3,a+3) = ( ll2_tr(a)*tau(b) - ll2_tr(b)*tau(a) ) /
     &                      ( ll2_tr(b) - ll2_tr(a))
         else
            Cstm(a+3,a+3) =  dtde(a,a) - dtde(b,a) - tau(a)
         end if
C
      end do ! a
C
C     Transform matrix to standard basis
C
      iop = 2
      call push3d (ipsw, iop, iwr, detF, nn_t, Sig, Cstm, Sig, Cst)
C      call tranr4(nn_t,nn_t,ff)
C      call pushr4(ff,ff,dd_l,dd,detf)
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
      write(iwr,9020) '*** ERROR IN SUBROUTINE PLAS3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      write(iwr,9020) 'iter   =', iter
      write(iwr,9020) 'lfirst =', lfirst
      write(iwr,9020) 'nSig   =', nSig
      write(iwr,9020) 'nTM    =', nTM
      write(iwr,9030) 'detF   =', detF
      call rprin0(pMAT, 1, 13, 'pMAT  ', iwr)
      call rprin0(Fn1 , 3,  3, 'Fn1   ', iwr)
      call rprin0(Fn  , 3,  3, 'Fn    ', iwr)
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
      write(iwr,9030) 'Epp    =', Epp
      call rprin0(be  , 1, nTM, 'be    ', iwr)
      call rprin0(Epl , 1,   3, 'Epl   ', iwr)
      call rprin0(Sig , 1,nSig, 'Sig   ', iwr)
      call rprin0(Cst , 6,   6, 'Cst   ', iwr)
                                                  go to 8010
C
C         Closing section
C
 8000 continue
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'LEAVING SUBROUTINE PLAS3D'
          if (ipsw .gt. 3)                        then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              write(iwr,9030) 'Epp    =', Epp
              call rprin0(be  , 1, nTM, 'be    ', iwr)
              call rprin0(Epl , 1,   3, 'Epl   ', iwr)
              call rprin0(Sig , 1,nSig, 'Sig   ', iwr)
              call rprin0(Cst , 6,   6, 'Cst   ', iwr)
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
 9030 format(3X,A,E15.6)
C
      end

      function yfunct(ipsw, iwr, pMAT, Epp, I1, J2, J3, f1, f2, f3,
     &                f11, f22, f33, f12, f13, f23, Ypr, YY, gam, flg)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   YFUNCT                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     YFUNCT: Compute yield function value
C
C METHOD:
C
C
C ARGUMENTS INPUT:
C     ipsw     - Print switch
C     iwr      - Write unit
C     pMAT     - Material properties:
C                pMAT(1) = Emod : Elastic Bulk  modulus
C                pMAT(2) = Poissons ratio
C                pMAT(5) = Bmod : Bulk modulus
C                pMAT(6) = Smod : Shear modulus
C                pMAT(7) = H_iso : Isotropic Hardening Modulus [times sqrt(2/3)]
C                pMAT(8) = H_kin : Kinematic Hardening Modulus [times (2/3)]
C                pMAT(9) = iYIELD (Yield Function number)
C                          iYIELD   = 1     : von Mises
C                          pMAT(10) = Y0    : Yield stress (initial)
C                          pMAT(11) = Y_inf : Yield stress (infinity)
C                          pMAT(12) = Delta : Hardening exponent
C                pMAT(13) = istrt : start state
C     Epp      - Cumulative plastic strain
C
C ARGUMENTS OUTPUT:
C     Yield    - Yield function value
C
C COMMON BLOCKS:
C     const
C
C PERIPHERAL UNITS:
C     None
C
C INTERNAL VARIABLES:
C     yield     - yield function value
C     yfunct    - yield function
C                 flg = .false.  only yield function
C                 flg = .true.   yield function and derivatives
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
      include 'include/feninc/const.h'
C
      integer ipsw, iwr, iYIELD
      real*8  yfunct, yield, gam
C
      real*8  pMAT(*), Epp
      real*8  I1, J2, J3
      real*8  f1, f2, f3, f11, f22, f33, f12, f13, f23, Ypr
      real*8  J2_12, J2_1, J2_2
      real*8  YY, aa, bb, ss, c1, c2
      real*8  beta, Hi, Y0, Yinf
      logical flg
C
      f1  = zero
      f2  = zero
      f3  = zero
      f11 = zero
      f22 = zero
      f33 = zero
      f12 = zero
      f13 = zero
      f23 = zero
C
      iYIELD = nint(pMAT(9))     ! Type of yield function
      Y0     = sqt23 * pMAT(10)  ! Radius of yield = sqrt(2/3) Sig_y
      Yinf   = sqt23 * pMAT(11)  ! Radius of yield = sqrt(2/3) Sig_inf
      beta   = pMAT(12)
      Hi     = sqt23 * pMAT(7)   ! Isotropic hardening = sqrt(2/3)*H_iso
      YY     = Y0 + Hi*Epp       ! Sig_y(t_n) yield stress
      Ypr    = two3 * pMAT(7)    ! Isotropic yield = 2/3 H_iso
      yield  =-YY                ! Default return = no yield
C
C     VON MISES YIELD FUNCTION
C
      if (iYIELD .eq. 1)                          then
C
         aa    = (Y0 - Yinf)*exp(-beta*Epp)
         YY    = Yinf + aa + Hi*Epp
         Ypr   = Ypr - sqt23*beta*aa
         ss    = sqrt( two * J2 )
         yield = ss - YY
C
C       Compute yield function derivatives
C
         if (flg)                                 then
C
            if (ss .ne. zero)                     then
               f2  =   one / ss
               f22 = - f2**3
            end if
C
         end if
C
C     DRUCKER PRAGER   YIELD FUNCTION
C
      else if (iYIELD .eq. 2)                     then
C
         aa    = Yinf
         ss    = sqrt( two * J2 )
         yield = ss + one3 * aa * I1 - YY
C
C       Compute yield function derivatives
C
         if (flg)                                 then
C
            f1 = one3 * aa
C
            if (ss .ne. zero)                     then
               f2  =   one / ss
               f22 = - f2**3
            end if
C
         end if
C
C     PRAGER-LODE      YIELD FUNCTION
C
      else if (iYIELD .eq. 3)                     then
C
         if (J2 .gt. zero)                        then
C
            if (flg)                              then
C
C           Read in material parameters
C
               bb = sqt23*pMAT(11)! Radius of yield = sqrt(2/3) Sig_inf
C
C           Compute yield function
C
               c1 = one / sqrt(two)
               c2 = sqrt(13.5d0) * bb
C
C           yield   = J2 * ( one + bb * J3 * J2 **(-onep5)) - YY
               yield = ( sqrt(two*J2) + c2 * J3 / J2 ) - YY
C
C           Compute yield function derivatives
C
               J2_1  =   one / J2
               J2_2  =   J2_1**2
               J2_12 =   sqrt(J2_1)
C
               f2    =   c1 * J2_12            - c2 * J3 * J2_2
               f3    =                           c2      * J2_1
               f22   = - c1 * J2_12**3 * one2 + c2 * J3 * J2_1**3 * two
               f23   =                         - c2      * J2_2
C
            end if
C
         else
C
            yield = -YY
C
         end if
C
      else
C
C     NO YIELD FUNCTION SPECIFIED
C
         yield = zero
         write (iwr,9000)
C
      end if
C
      yfunct = yield
C
      return
C
C     Format
C
 9000  format(' *** NO YIELD FUNCTION SPECIFIED ***')
C
      end
