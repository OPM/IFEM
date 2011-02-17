      subroutine wlog3d (mu, bpr, epsd, taup, aap, w,
     &                   ipswb3, iwr, ierr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   WLOG3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     WLOG3D: Compute stress and tangent in principal stretch form
C
C METHOD:
C     Stored Energy Function in Principal Stretches: 
C          W       = w(lam_1) + w(lam_2) + w(lam_3)
C          taup_i  = [w(lam_i)]' * lambda_i
C          gamm_ij = lambda_j * [taup_i],j
C
C     Logarithmic strain energy function
C          w    = mu * [log(lambda)]**2      
C
C ARGUMENTS INPUT:  
C     mu       - Shear modulus
C     bpr(3)   - Principal stretch (squared)
C
C ARGUMENTS OUTPUT:
C     epsd(3)  - Principal deviatoric strains
C     taup(3)  - Principal deviatoric Kirchhoff stresses
C     aap(6,6) - Deviatoric tangent moduli in principal basis
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
C          
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
      integer   i, j, ipswb3, iwr, ierr
C
      real*8    jthird, mu, twomu, fourmu, w, tol,
     &          bpr(3), lamt(3), epsd(3), taup(3), aap(6,6)
C
      include 'const.h'
C
      data      tol / 1.0d-08 /
C
C         Entry section
C
      ierr = 0
C
      if (ipswb3 .gt. 0)                          then 
          write(iwr,9010) 'ENTERING SUBROUTINE WLOG3D'
          if (ipswb3 .gt. 2)                      then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9030) 'mu     =', mu
              call rprin0(bpr, 1, 3, 'bpr   ', iwr)
          endif
      endif
C
C         Compute Je**(-1/3) and volume preserving stretches
C
      jthird    = (bpr(1)*bpr(2)*bpr(3))**(-one6)
      do      i = 1,3
        lamt(i) = jthird*sqrt(bpr(i))
      end do ! i
C
C         Compute predictor strain energy function derivatives
C
      twomu     = two * mu
      w         = zero
      do      i = 1,3
        epsd(i) = log(lamt(i))
        taup(i) = twomu * epsd(i)
        w       = w +  mu * epsd(i)**2
      end do ! i
C
C         Load material part of moduli in principal stretches
C
      aap    = zero
C
      twomu  = twomu * one3
      fourmu = twomu + twomu
C
      do       i = 1,3
        aap(i,i) = fourmu - two * taup(i)
      end do ! i
C
      aap(1,2) = - twomu
      aap(2,1) = - twomu
C
      aap(2,3) = - twomu
      aap(3,2) = - twomu
C
      aap(3,1) = - twomu
      aap(1,3) = - twomu
C
C         Compute deviatoric tangent moduli in principal basis
C
      do i = 1,3
         j = mod(i,3) + 1
C
         if (abs(bpr(i) - bpr(j)) .gt. tol) then
          aap(i+3,i+3) = ( bpr(i)*taup(j) - bpr(j)*taup(i) )
     &                 / ( bpr(j) - bpr(i) )
        else
          aap(i+3,i+3) = (aap(i,i) - aap(j,i))*one2
        endif
C
      end do ! i
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
      write(iwr,9020) '*** ERROR IN SUBROUTINE WLOG3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      write(iwr,9030) 'mu     =', mu
      call rprin0(bpr, 1, 3, 'bpr   ', iwr)
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
      write(iwr,9030) 'w      =', w  
      call rprin0(epsd, 1, 3, 'epsd  ', iwr)
      call rprin0(taup, 1, 3, 'taup  ', iwr)
      call rprin0(aap , 6, 6, 'aap   ', iwr)
                                                  go to 8010
C
C         Closing section
C
 8000 continue
C
      if (ipswb3 .gt. 0)                          then
          write(iwr,9010) 'LEAVING SUBROUTINE WLOG3D'
          if (ipswb3 .gt. 3)                      then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              write(iwr,9030) 'w      =', w  
              call rprin0(epsd, 1, 3, 'epsd  ', iwr)
              call rprin0(taup, 1, 3, 'taup  ', iwr)
              call rprin0(aap , 6, 6, 'aap   ', iwr)
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
