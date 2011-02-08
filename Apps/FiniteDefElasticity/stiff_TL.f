C $Id: stiff_TL.f,v 1.1 2010-12-18 20:03:54 kmo Exp $
C=======================================================================
C!
C! \file stiff_TL.f
C!
C! \date Jun 8 2010
C!
C! \author Knut Morten Okstad / SINTEF
C!
C! \brief Routines for Total Lagrangian material stiffness integration.
C!
C=======================================================================
C
      subroutine stiff_TL2D (NENOD,detJW,dNdX,F,Cmat,EK)
C
      implicit none
C
      integer          NENOD, a, b, i, j, k, l, m, n
      double precision detJW, dNdX(NENOD,2), F(2,2), Cmat(3,3)
      double precision EK(2,NENOD,2,NENOD), C(2,2,2,2), km
C
      do 200 i = 1, 2
       do 200 j = i, 2
         if (i == j) then
            a = i
         else
            a = 3
         end if
         do 100 k = 1, 2
          do 100 l = k, 2
            if (k == l) then
               b = k
            else
               b = 3
            end if
            C(i,j,k,l) = Cmat(a,b)
            C(j,i,k,l) = Cmat(a,b)
            C(i,j,l,k) = Cmat(a,b)
            C(j,i,l,k) = Cmat(a,b)
  100    continue
  200 continue
C
      do 450 a = 1, NENOD
       do 450 b = a, NENOD
         do 400 m = 1, 2
          do 400 n = 1, 2
            km = 0.0D0
            do 350 i = 1, 2
             do 350 j = 1, 2
              do 300 k = 1, 2
               do 300 l = 1, 2
                  km = km + dNdX(a,i)*F(m,j)*C(i,j,k,l)*F(n,k)*dNdX(b,l)
  300          continue
  350       continue
            EK(m,a,n,b) = EK(m,a,n,b) + km*detJW
            if (b > a) then
               EK(n,b,m,a) = EK(n,b,m,a) + km*detJW
            end if
  400    continue
  450 continue
C
      end
C
C
      subroutine stiff_TL3D (NENOD,detJW,dNdX,F,Cmat,EK)
C
      implicit none
C
      integer          NENOD, a, b, i, j, k, l, m, n
      double precision detJW, DNDX(NENOD,3), F(3,3), Cmat(6,6)
      double precision EK(3,NENOD,3,NENOD), C(3,3,3,3), km
C
      do 200 i = 1, 3
       do 200 j = i, 3
         if (i == j) then
            a = i
         else if (i == 1) then
            a = j+2
         else
            a = 6
         end if
         do 100 k = 1, 3
          do 100 l = k, 3
            if (k == l) then
               b = k
            else if (k == 1) then
               b = l+2
            else
               b = 6
            end if
            C(i,j,k,l) = Cmat(a,b)
            C(j,i,k,l) = Cmat(a,b)
            C(i,j,l,k) = Cmat(a,b)
            C(j,i,l,k) = Cmat(a,b)
  100    continue
  200 continue
C
      do 450 a = 1, NENOD
       do 450 b = a, NENOD
         do 400 m = 1, 3
          do 400 n = 1, 3
            km = 0.0D0
            do 350 i = 1, 3
             do 350 j = 1, 3
              do 300 k = 1, 3
               do 300 l = 1, 3
                  km = km + dNdX(a,i)*F(m,j)*C(i,j,k,l)*F(n,k)*dNdX(b,l)
  300          continue
  350       continue
            EK(m,a,n,b) = EK(m,a,n,b) + km*detJW
            if (b > a) then
               EK(n,b,m,a) = EK(n,b,m,a) + km*detJW
            end if
  400    continue
  450 continue
C
      end
C
C
      subroutine stiff_TL2D_isoel (NENOD,detJW,dNdX,F,C1,C2,C3,EK)
C
C     Use this subroutine in 2D, but only for isotrophic materials,
C     where the fourth-order constitutive tensor is assumed to only
C     contain three different values, C1, C2, C3 (and zeroes):
C
C     C1 = C(1,1,1,1) = C(2,2,2,2)
C     C2 = C(1,1,2,2) = C(2,2,1,1
C     C3 = C(1,2,1,2) = C(1,2,2,1) = C(2,1,1,2) = C(2,1,2,1)
C
      implicit none
C
      integer          NENOD, a, b, m, n
      double precision detJW, DNDX(NENOD,2), F(2,2), C1, C2, C3
      double precision EK(2,NENOD,2,NENOD), ek1, ek2, ek3
C
      do 350 a = 1, NENOD
       do 350 b = a, NENOD
         do 300 m = 1, 2
          do 300 n = 1, 2
            ek1 = 0.0D0
            ek2 = 0.0D0
            ek3 = 0.0D0
            ek1 = dNdX(a,1)*F(m,1)*F(n,1)*dNdX(b,1)
     +          + dNdX(a,2)*F(m,2)*F(n,2)*dNdX(b,2)
            ek2 = dNdX(a,1)*F(m,1)*F(n,2)*dNdX(b,2)
     +          + dNdX(a,2)*F(m,2)*F(n,1)*dNdX(b,1)
            ek3 = dNdX(a,1)*F(m,2)*(F(n,1)*dNdX(b,2) +
     +                              F(n,2)*dNdX(b,1))
     +          + dNdX(a,2)*F(m,1)*(F(n,2)*dNdX(b,1) +
     +                              F(n,1)*dNdX(b,2))
            EK(m,a,n,b) = EK(m,a,n,b) + (ek1*C1+ek2*C2+ek3*C3)*detJW
            if (b > a) then
               EK(n,b,m,a) = EK(n,b,m,a) + (ek1*C1+ek2*C2+ek3*C3)*detJW
            end if
  300    continue
  350 continue
C
      end
C
C
      subroutine stiff_TL3D_isoel (NENOD,detJW,dNdX,F,C1,C2,C3,EK)
C
C     Use this subroutine in 3D, but only for isotrophic materials,
C     where the fourth-order constitutive tensor is assumed to only
C     contain three different values, C1, C2, C3 (and zeroes):
C
C     C1 = C(1,1,1,1) = C(2,2,2,2) = C(3,3,3,3)
C     C2 = C(1,1,2,2) = C(1,1,3,3) = C(2,2,3,3)
C        = C(2,2,1,1) = C(3,3,1,1) = C(3,3,2,2)
C     C3 = C(1,2,1,2) = C(1,2,2,1) = C(2,1,1,2) = C(2,1,2,1)
C        = C(1,3,1,3) = C(1,3,3,1) = C(3,1,1,3) = C(3,1,3,1)
C        = C(2,3,2,3) = C(2,3,3,2) = C(3,2,2,3) = C(3,2,3,2)
C
      implicit none
C
      integer          NENOD, a, b, i, j, m, n
      double precision detJW, DNDX(NENOD,3), F(3,3), C1, C2, C3
      double precision EK(3,NENOD,3,NENOD), ek1, ek2, ek3
C
      do 350 a = 1, NENOD
       do 350 b = a, NENOD
         do 300 m = 1, 3
          do 300 n = 1, 3
            ek1 = 0.0D0
            ek2 = 0.0D0
            ek3 = 0.0D0
            do 200 i = 1, 3
               ek1 = ek1 + dNdX(a,i)*F(m,i)*F(n,i)*dNdX(b,i)
               do 100 j = i+1, 3
                  ek2 = ek2 + dNdX(a,i)*F(m,i)*F(n,j)*dNdX(b,j)
     +                      + dNdX(a,j)*F(m,j)*F(n,i)*dNdX(b,i)
                  ek3 = ek3 + dNdX(a,i)*F(m,j)*(F(n,i)*dNdX(b,j) +
     +                                          F(n,j)*dNdX(b,i))
     +                      + dNdX(a,j)*F(m,i)*(F(n,j)*dNdX(b,i) +
     +                                          F(n,i)*dNdX(b,j))
  100          continue
  200       continue
            EK(m,a,n,b) = EK(m,a,n,b) + (ek1*C1+ek2*C2+ek3*C3)*detJW
            if (b > a) then
               EK(n,b,m,a) = EK(n,b,m,a) + (ek1*C1+ek2*C2+ek3*C3)*detJW
            end if
  300    continue
  350 continue
C
      end
