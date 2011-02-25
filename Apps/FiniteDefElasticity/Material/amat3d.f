      subroutine amat3d (ipsw, iwr, taup, aap, q, tau, aa, ierr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   AMAT3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     AMAT3D: Transform Kirchhoff stresses and material moduli from 
C             principal to standard basis
C
C METHOD:
C     
C
C ARGUMENTS INPUT:
C     ipsw     - Print switch
C     iwr      - Write unit number
C     taup(3)  - Principal deviatoric Kirchhoff stresses
C     aap(6,6) - Deviatoric tangent moduli in principal basis
C     q(6,6)   - Transformation matrix: [principal]-->[std] basis
C
C ARGUMENTS OUTPUT:
C     tau(3)   - Deviatoric Kirchhoff stresses (ref to standard basis)
C     aa(6,6)  - Tangent moduli in standard basis
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
      integer   i, j, k, ipsw, iwr, ierr
C
      real*8    taup(3), aap(6,6), q(6,6), tau(6), aa(6,6), vi(6)
C
      include 'include/feninc/const.h'
C
C         Entry section
C
      ierr = 0
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'ENTERING SUBROUTINE AMAT3D'
          if (ipsw .gt. 2)                        then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              call rprin0(taup, 1, 3, 'taup  ', iwr)
              call rprin0(aap , 6, 6, 'aap   ', iwr)
              call rprin0(q   , 6, 6, 'q     ', iwr)
          endif
      endif
C
C         Compute stress tensor
C
      do j = 1,6
        tau(j) = q(j,1)*taup(1) + q(j,2)*taup(2) + q(j,3)*taup(3)
      end do ! j
C
C         Compute tangent matrix in standard basis
C
      do i = 1,6
C
C         Left multiplication: V = Q * AL
C
        do      k = 1,3
          vi(k)   = zero
          do    j = 1,3
            vi(k) = vi(k) + q(i,j)*aap(j,k)
          end do ! j
        end do ! k
C
        do k = 4,6
            vi(k) = q(i,k)*aap(k,k)
        end do ! k
C
C         Right multiplication: AA = V * Q_t
C
        do        j = 1,i
          aa(i,j)   = zero
          do      k = 1,6
            aa(i,j) = aa(i,j) + vi(k)*q(j,k)
          end do ! k
C
C         Fill in symmetric part
C
          aa(j,i) = aa(i,j)
C
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
      write(iwr,9020) '*** ERROR IN SUBROUTINE AMAT3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      call rprin0(taup, 1, 3, 'taup  ', iwr)
      call rprin0(aap , 6, 6, 'aap   ', iwr)
      call rprin0(q   , 6, 6, 'q     ', iwr)
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS' 
      call rprin0(tau, 1, 3, 'tau   ', iwr)
      call rprin0(aa , 6, 6, 'aa    ', iwr)
                                                  go to 8010
C
C         Closing section
C
 8000 continue
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'LEAVING SUBROUTINE AMAT3D'
          if (ipsw .gt. 3)                        then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              call rprin0(tau, 1, 3, 'tau   ', iwr)
              call rprin0(aa , 6, 6, 'aa    ', iwr)
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
