      subroutine eigs3d (ipsw, iwr, v, d, nrotd, ierr)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   EIGS3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     EIGS3D: Compute eigenvalues/vectors for 3x3 symmetric matrix
C
C METHOD:
C     Storage done as follows:
C
C       | v(1,1) v(1,2) v(1,3) |     |  d(1)  a(1)  a(3)  |
C       | v(2,1) v(2,2) v(2,3) |  =  |  a(1)  d(2)  a(2)  |
C       | v(3,1) v(3,2) v(3,3) |     |  a(3)  a(2)  d(3)  |
C
C     Transformations performed on d(i) and a(i) and v(i,j) become
C     eigenvectors.  Thus, original array is destroyed.
C
C ARGUMENTS INPUT:  
C     v(3,3) - Matrix for which eigenvalues/vectors should be 
C              computed for
C
C ARGUMENTS OUTPUT:
C     ipsw   - Print switch
C     iwr    - Write unit number
C     v(3,3) - Matrix of eigenvectors (by column)
C     d(3)   - Eigenvalues associated with the columns of v
C     nrotd  - number of rotations to diagonalize
C     ierr   - Error flag
C              = -1 : 
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
      integer   nrotd, i, its, j, k, ipsw, iwr, ierr
C
      real*8    g, h, aij, sm, thresh, t, c, s, tau,
     *          v(3,3), d(3), a(3), b(3), z(3)
C
      include 'include/feninc/const.h'
C
C         Entry section
C
      ierr = 0
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'ENTERING SUBROUTINE EIGS3D'
          if (ipsw .gt. 2)                        then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              call rprin0(v, 3, 3, 'v     ', iwr)
          endif
      endif
C
C         Move array into one-d arrays
C
      a(1) = v(1,2)
      a(2) = v(2,3)
      a(3) = v(1,3)
C
      do i = 1,3
C
         d(i) = v(i,i)
         b(i) = d(i)
         z(i) = zero
C
         do j = 1,3
C
          v(i,j) = zero
C
         end do ! j
C
         v(i,i) = one
C
      end do ! i
C
C         Check for diagonal case
C
      sm = abs(a(1)) + abs(a(2)) + abs(a(3))
      g  = abs(d(1)) + abs(d(2)) + abs(d(3))
C
      if (sm .lt. 1.d-13*g)                       go to 8000
C
      nrotd = 0
C
      do its = 1, 50
C
C         Set convergence test and threshold
C
        sm = abs(a(1)) + abs(a(2)) + abs(a(3))
        if (sm .eq. zero)                         go to 8000
C
        if (its .lt. 4)                           then
C
          thresh = 0.011d0*sm
C
        else
C
          thresh = zero
C
        end if
C
C         Perform sweeps for rotations
C
        do i = 1,3
C
           j = mod(i,3) + 1
           k = mod(j,3) + 1
C
           aij  = a(i)
           g    = 100.d0*abs(aij)
C
           if (abs(d(i))+g .ne. abs(d(i)) .or.
     &         abs(d(j))+g .ne. abs(d(j)))        then
C
              if (abs(aij) .gt. thresh)           then
C
                 a(i) = zero
                 h    = d(j) - d(i)
C
                 if (abs(h)+g .eq. abs(h))        then
                    t = aij/h
                 else
                    t = sign(two,h/aij)/(abs(h/aij)
     &                     + sqrt(four+(h/aij)**2))
                 end if
C
C         Set rotation parameters
C
                 c    = one/sqrt(one+t*t)
                 s    = t*c
                 tau  = s/(one+c)
C
C         Rotate diagonal terms
C
                 h    = t*aij
                 z(i) = z(i) - h
                 z(j) = z(j) + h
                 d(i) = d(i) - h
                 d(j) = d(j) + h
C
C         Rotate off-diagonal terms
C
                 h    = a(j)
                 g    = a(k)
                 a(j) = h + s*(g - h*tau)
                 a(k) = g - s*(h + g*tau)
C
C         Rotate eigenvectors
C
                 do k      = 1,3
C
                    g      = v(k,i)
                    h      = v(k,j)
                    v(k,i) = g - s*(h + g*tau)
                    v(k,j) = h + s*(g - h*tau)
C
                 end do ! k
C
                 nrotd = nrotd + 1
C
              end if
C
           else
C
              a(i) = zero
C
           end if
C
        end do ! i
C
C         Update diagonal terms
C
        do i   = 1,3
C
          b(i) = b(i) + z(i)
          d(i) = b(i)
          z(i) = zero
C
        end do ! i
C
      end do ! its
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
      write(iwr,9020) '*** ERROR IN SUBROUTINE EIGS3D *** IERR = ', ierr
      write(iwr,9010) 'WITH INPUT ARGUMENTS'
      call rprin0(v, 3, 3, 'v     ', iwr)
C
      write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
      write(iwr,9020) 'nrotd  =', nrotd 
      call rprin0(v, 3, 3, 'v     ', iwr)
      call rprin0(d, 1, 3, 'd     ', iwr)
                                                  go to 8010
C
C         Closing section
C
 8000 continue
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'LEAVING SUBROUTINE EIGS3D'
          if (ipsw .gt. 3)                        then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              write(iwr,9020) 'nrotd  =', nrotd 
              call rprin0(v, 3, 3, 'v     ', iwr)
              call rprin0(d, 1, 3, 'd     ', iwr)
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
