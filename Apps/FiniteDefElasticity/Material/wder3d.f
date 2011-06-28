      subroutine wder3d (ipsw, iwr, nn, pMAT, bpr, taup, aap, wengy, 
     &                   ierr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   WDER3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     WDER3D: Compute stress and tangent in principal stretch form
C
C METHOD:
C     Stored Energy Function in Principal Stretches: 
C          W       = w(xlam_1) + w(xlam_2) + w(xlam_3)
C          xlam _i = (Je)^(-1/3)*lambda_i (deviatoric stretches)
C
C     Ogden strain energy functions in terms of Seth strains
C          w_i  = c-alpha/(m-alpha)*(xlam_i**m-alpha - 1.) -> wengy     
C
C ARGUMENTS INPUT:
C     ipsw     - Print switch
C     iwr      - Write unit number 
C     nn       - Number of terms in series expansion 
C                of W(xlam) = 3, max
C     pMAT     - Ogden material parameters
C                pMAT(1,i) = ci
C                pMAT(2,i) = ni
C     bpr(3)   - Principal values of left Cauchy-Green tensor
C
C ARGUMENTS OUTPUT:
C     taup(3)  - Principal deviatoric Kirchhoff stresses
C     aap(6,6) - Deviatoric tangent moduli in principal basis
C     wengy    - Stored deviatoric strain energy
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
      integer   i, j, k, n, nn, ipsw, iwr, ierr
C
      real*8    jthrd, c1, c2, c3, c4, wengy, tol, pMAT(2,*),
     &          bpr(3), xlamt(3), xlam2(3), taup(3), aap(6,6),
     &          dwtil(4,3), dtdl(3,3)
C
      include 'include/feninc/const.h'
C
      data      tol / 1.0d-06 /
C
C         Entry section
C
      ierr = 0
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'ENTERING SUBROUTINE WDER3D'
          if (ipsw .gt. 2)                        then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9020) 'nn     =', nn
              call rprin0(pMAT, 2, nn, 'pMAT  ', iwr)
              call rprin0(bpr , 1, 3   , 'bpr   ', iwr)
          endif
      endif
C
      if (nn .lt. 1 .or. nn .gt. 3)               go to 7100
C
C         Initialize arrays
C
      dwtil = zero
      dtdl  = zero
      aap   = zero
C
C        Compute J**(-1/3)
C
      jthrd = (bpr(1)*bpr(2)*bpr(3))**(-one6)
C
C        Compute modified (Flory) stretches
C
      do i = 1,3
         xlamt(i) = jthrd*sqrt(bpr(i))
         xlam2(i) = xlamt(i)*xlamt(i)
      end do ! i
C
C        Compute strain energy function derivatives
C
      wengy = zero
C
      do n = 1, nn
C
         c1   = pMAT(1,n)
         c2   = pMAT(2,n)
C
         do i = 1, 3
C
            wengy      = wengy      + (c1/c2)*(xlamt(i)**c2 - one)
            dwtil(1,i) = dwtil(1,i) + c1*xlamt(i)**(c2-one)
            dwtil(2,i) = dwtil(2,i) + (c2-one)*c1*xlamt(i)**(c2-two)
            dwtil(3,i) = dwtil(3,i)
     &                  + (c2-two)*(c2-one)*c1*xlamt(i)**(c2-three)
            dwtil(4,i) = dwtil(4,i)
     &        + (c2-three)*(c2-two)*(c2-one)*c1*xlamt(i)**(c2-four)
        end do ! i
C
      end do ! n
C
C         Compute deviatoric principal Kirchoff stresses and
C         their derivatives wrt xlamt_j

      do i = 1, 3
C
         taup(i) = xlamt(i)*dwtil(1,i) - ( xlamt(1)*dwtil(1,1)
     &                                   + xlamt(2)*dwtil(1,2)
     &                                   + xlamt(3)*dwtil(1,3) )*one3
         dtdl(i,i) = dwtil(2,i)*xlamt(i) + dwtil(1,i)
C
         do j = 1,3
            dtdl(i,j) = dtdl(i,j)
     &                - ( dwtil(2,j)*xlamt(j) + dwtil(1,j) )*one3
         end do ! j
C
      end do ! i
C
C         Compute deviatoric tangent matrix in principal basis
C
      do i = 1, 3
C
         k = 1 + mod(i,3)
C
         if (abs(xlamt(i)-xlamt(k)) .gt. tol)     then
C
            aap(i+3,i+3) =
     &         ( xlam2(i)*taup(k) - xlam2(k)*taup(i) )/
     &         ( xlam2(k) - xlam2(i) )
         else
            c3 = one/(xlamt(i)+xlamt(k))
            c4 = xlamt(k)-xlamt(i)
            aap(i+3,i+3) = c3*xlamt(i)*xlamt(k) * ( - dwtil(1,i)
     &                   + dwtil(2,i)*xlamt(i)
     &                   + dwtil(3,i)*xlamt(i)*c4*one2
     &                   + dwtil(4,i)*xlamt(i)*c4**2*one6 )
     &                 + ( dwtil(1,1)*xlamt(1)
     &                   + dwtil(1,2)*xlamt(2)
     &                   + dwtil(1,3)*xlamt(3) )*one3
         endif
C
         aap(i,i) = -two*taup(i)
C
         do j = 1,3
C
            aap(i,j) = aap(i,j) + xlamt(j)*dtdl(i,j)
     &               - ( xlamt(1)*dtdl(i,1)
     &                 + xlamt(2)*dtdl(i,2)
     &                 + xlamt(3)*dtdl(i,3) )*one3
         end do ! j
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
      write(iwr,9020) '*** ERROR IN SUBROUTINE WDER3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      write(iwr,9020) 'nn     =', nn
      call rprin0(pMAT, 1, 2*nn, 'pMAT  ', iwr)
      call rprin0(bpr , 1, 3   , 'bpr   ', iwr)
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
      write(iwr,9030) 'wengy  =', wengy  
      call rprin0(taup, 1, 3, 'taup  ', iwr)
      call rprin0(aap , 6, 6, 'aap   ', iwr)
                                                  go to 8010
C
C         Closing section
C
 8000 continue
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'LEAVING SUBROUTINE WDER3D'
          if (ipsw .gt. 3)                        then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              write(iwr,9030) 'wengy  =', wengy  
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
