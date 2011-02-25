      subroutine push3d (ipsw, iop, iwr, detF, F, PK2, Cmt, Sig, Cst)
C
C ---------------------------------------------------------------------
C
C F E N R I S   ROUTINE:   PUSH3D                        NTNU/SINTEF/DNV
C
C PURPOSE:
C     PUSH3D: Push stresses and constitutive tensor forward from
C             material to spatial coordinates in current configuration
C
C METHOD:
C     Material and spatial stress components are stored as:
C
C                   | Stress(1)  Stress(4)  Stress(6) |
C     Stress(3,3) = | Stress(4)  Stress(2)  Stress(5) |
C                   | Stress(6)  Stress(5)  Stress(3) |
C
C ARGUMENTS INPUT:
C     ipsw     - Prints switch
C     iop      - Optional parameter
C                = 0 : Push stresses and constitutive tensor
C                = 1 : Push stresses only
C                = 2 : Push constitutive tensor only
C     iwr      - Write unit number
C     detF     - Determinant of deformation gradient
C     F(3,3)   - Deformation gradient
C     PK2(6)   - 2nd Piola-Kirchhoff stresses (material stresses)
C     Cmt(6,6) - Material constitutive tensor
C
C ARGUMENTS OUTPUT:
C     Sig(6)   - Cauchy stresses (spatial stresses)
C     Cst(6,6) - Spatial constitutive tensor
C
C COMMON BLOCKS:
C     const
C
C PERIPHERAL UNITS:
C     None
C
C INTERNAL VARIABLES:
C     Tmat(6,6) - Transformation matrix
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
      integer   iop, ipsw, iwr
      integer   i, j, i1(6), i2(6)
      real*8    detF, F(3,3), PK2(6), Cmt(6,6), Sig(6), Cst(6,6)
      real*8    detFi, Tmat(6,6)
C
      include 'include/feninc/const.h'
C
      data      i1 /1,2,3,1,2,3/
      data      i2 /1,2,3,2,3,1/
C
C         Entry section
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'ENTERING SUBROUTINE PUSH3D'
          if (ipsw .gt. 2)                        then
              write(iwr,9010) 'WITH INPUT ARGUMENTS'
              write(iwr,9020) 'iop    =', iop
              write(iwr,9030) 'detF   =', detF
              call rprin0(F  , 3, 3,'F     ',iwr)
              if (iop .eq. 0 .or. iop .eq. 1)     then
                 call rprin0(PK2, 1, 6,'PK2   ',iwr)
              end if
              if (iop .eq. 0 .or. iop .eq. 2)     then
                 call rprin0(Cmt, 6, 6,'Cmt   ',iwr)
              end if
          end if
      end if
C
      detFi = one / detF
C
      if (iop .eq. 0 .or. iop .eq. 1)             then
C
C         Push forward material stress (PK2) to current configuration
C
C         Sig(i,j) = F(i,k)*PK2(k,l)*F(j,l)/detf
C
C         Tmat = F^T * PK2
C
         do i = 1,3
            Tmat(i,1) = F(i,1)*PK2(1) + F(i,2)*PK2(4) + F(i,3)*PK2(6)
            Tmat(i,2) = F(i,1)*PK2(4) + F(i,2)*PK2(2) + F(i,3)*PK2(5)
            Tmat(i,3) = F(i,1)*PK2(6) + F(i,2)*PK2(5) + F(i,3)*PK2(3)
         end do ! i
C
C         Sig = Tmat*F
C
         Sig(1) = Tmat(1,1)*F(1,1) + Tmat(1,2)*F(1,2) + Tmat(1,3)*F(1,3)
         Sig(2) = Tmat(2,1)*F(2,1) + Tmat(2,2)*F(2,2) + Tmat(2,3)*F(2,3)
         Sig(3) = Tmat(3,1)*F(3,1) + Tmat(3,2)*F(3,2) + Tmat(3,3)*F(3,3)
         Sig(4) = Tmat(1,1)*F(2,1) + Tmat(1,2)*F(2,2) + Tmat(1,3)*F(2,3)
         Sig(5) = Tmat(2,1)*F(3,1) + Tmat(2,2)*F(3,2) + Tmat(2,3)*F(3,3)
         Sig(6) = Tmat(3,1)*F(1,1) + Tmat(3,2)*F(1,2) + Tmat(3,3)*F(1,3)
C
         Sig   = detFi * Sig
C
      end if
C
      if (iop .eq. 0 .or. iop .eq. 2)             then
C
C         Form transformation array for a 4th rank tensor in matrix form
C
C         Tmat(a,b) = F(i,I)*F(j,J) : a -> I,J ; b -> i,j
C              a,b  |  1    2    3    4    5    6
C             ------+-----------------------------
C             (I,J) | 1,1  2,2  3,3  1,2  2,3  3,1
C          or (i,j) |                2,1  3,2  1,3
C
         do i = 1,3
            do j = 1,3
               Tmat(i,j) =  F(i1(j),i1(i))*F(i2(j),i2(i))
            end do ! j
            do j = 4,6
               Tmat(i,j) = (F(i1(j),i1(i))*F(i2(j),i2(i))
     &                   +  F(i2(j),i2(i))*F(i1(j),i1(i)))*one2
            end do ! j
         end do ! i
C
         do i = 4,6
            do j = 1,3
               Tmat(i,j) =  F(i1(j),i1(i))*F(i2(j),i2(i))
     &                   +  F(i2(j),i2(i))*F(i1(j),i1(i))
            end do ! j
            do j = 4,6
               Tmat(i,j) = (F(i1(j),i1(i))*F(i2(j),i2(i))
     &                   +  F(i2(j),i1(i))*F(i1(j),i2(i))
     &                   +  F(i1(j),i2(i))*F(i2(j),i1(i))
     &                   +  F(i2(j),i2(i))*F(i1(j),i1(i)))*one2
            end do ! j
         end do ! i
C
C         Compute matrix product: Cst = Tmat^t * Cmt * Tmat
C
         Cst = matmul(transpose(Tmat),matmul(Cmt,Tmat))
C
         Cst = detFi * Cst
C
      end if
C
C         Closing section
C
      if (ipsw .gt. 0)                            then
          write(iwr,9010) 'LEAVING SUBROUTINE PUSH3D'
          if (ipsw .gt. 3)                        then
              write(iwr,9010) 'WITH OUTPUT ARGUMENTS'
              if (iop .eq. 0 .or. iop .eq. 1)     then
                 call rprin0(Sig, 1, 6,'Sig   ',iwr)
              end if
              if (iop .eq. 0 .or. iop .eq. 2)     then
                 call rprin0(Cmt, 6, 6,'Cmt   ',iwr)
              end if
          endif
      endif
C
      return
C
 9010 format(3X,A)
 9020 format(3X,A,I6)
 9030 format(3X,A,E15.6)
C
      end
