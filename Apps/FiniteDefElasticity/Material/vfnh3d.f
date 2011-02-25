      subroutine vfnh3d (ipsw, iwr, iVF, detF, detFi, lambda, U, Up,
     &                   Upp, ierr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   VFNH3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     VFNH3D: Compute volumetric stresses (pressure) and material moduli
C             according to the volumetric function type
C
C METHOD:
C     Volumetric function type: U_1 = K*[ 0.25*( J^2 - 1 )-0.5*ln(J)]
C                               U_2 = K*[ 0.5*( J - 1 )^2 ]
C                               U_3 = K*[ 0.5*( ln(J) )^2 ]
C                               U_4 = K*[ 2.0*(J - 1 - ln(J)) ]
C     Free energy function:
C     U     = U(J)
C     Up    = partial_J ( Psi )
C     Upp   = [ partial_J ( J partial_J Psi ) - up ]/J
C           =   partial^2_J ( Psi )/JC
C     
C
C ARGUMENTS INPUT:  
C     ipsw    - Prints switch
C     iwr     - Write unit number
C     iVF     - Volumetric function type
C     detF    - Determinant of the deformation gradient
C     detFi   - Inverse of the determinant of the Jacobian gradient
C     lambda  - Lame's constant
C
C ARGUMENTS OUTPUT:
C     U       - Strain energy density
C     Up      - Mean stress (pressure)
C     Upp     - Volumetric constitutive tensor
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
C     U     = U(J)     
C     Up    = partial_J ( Psi )
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
C
      real*8    detF, detFi, lambda,
     &          detF2, detFi2, detFm1, U, Up, Upp
C
      include 'include/feninc/const.h'
C
C         Entry section
C
      ierr = 0
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'ENTERING SUBROUTINE VFNH3D'
          if (ipsw .gt. 2)                        then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9020) 'iVF    =', iVF
              write(iwr,9030) 'detF   =', detF
              write(iwr,9030) 'detFi  =', detFi
              write(iwr,9030) 'lambda =', lambda
          endif
      endif
C
C         Compute volumetric stresses (pressure) and material moduli
C         according to the volumetric function type
C
      if (iVF . eq. 1)                            then
C
C         Volumetric function type 1: 
C         U(J) = lambda*0.25*(J^2 - 1 - 2*(log J))
C
        detF2  = detF**2
        detFi2 = detFi**2
        U      = one4 * lambda * ( detF2 - one  - two*log(abs(detF)) )
        Up     = one2 * lambda * ( detF - detFi )
        Upp    = one2 * lambda * ( one + detFi2 )
C
      else if (iVF .eq. 2)                        then
C
C         Volumetric function type 2: 
C         U(J) = lambda*0.5*(J-1)^2
C
        detFm1 = detF - one
        U      = one2 * lambda * (detFm1)**2
        Up     =        lambda *  detFm1
        Upp    =        lambda
C
      else if (IVF .eq. 3)                        then
C
C         Volumetric function type 3: 
C         U(J) = lambda*0.5*(log J)^2
C
        U    =  one2 * lambda * log(abs(detF))**2
        Up   =         lambda * log(abs(detF)) * detFi
        Upp  =       ( lambda * detFi - Up ) * detFi
C
      else if (iVF .eq. 4)                        then
C
C         Volumetric function type 4: 
C         U(J) = lambda*2.0*(J - 1 - log J)
C
        detFi2 = detFi**2 
        U      = two * lambda * ( detF - one - log(abs(detF)) )
        Up     =       lambda * ( two - detFi )
        Upp    =       lambda * detFi2
C
C         Illegal volumetric function
C
      else
                                                  go to 7000
      end if
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
      write(iwr,9020) '*** ERROR IN SUBROUTINE VFNH3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      write(iwr,9030) 'detF   =', detF
      write(iwr,9030) 'detFi  =', detFi
      write(iwr,9030) 'lambda =', lambda
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
      write(iwr,9030) 'U      =', U  
      write(iwr,9030) 'Up     =', Up
      write(iwr,9030) 'Upp    =', Upp
                                                  go to 8010
C
C         Closing section
C
 8000 continue
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'LEAVING SUBROUTINE VFNH3D'
          if (ipsw .gt. 3)                        then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              write(iwr,9030) 'U      =', U  
              write(iwr,9030) 'Up     =', Up
              write(iwr,9030) 'Upp    =', Upp
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
