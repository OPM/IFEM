C $Id: accKM.f,v 1.1 2010-12-29 18:51:53 kmo Exp $
C=======================================================================
C!
C! \file accKM.f
C!
C! \date Dec 18 2010
C!
C! \author Knut Morten Okstad / SINTEF
C!
C! \brief Material stiffness accumulation routines.
C!
C=======================================================================
C
      subroutine accKM2D (NENOD,Shpf,D,EKt)
C
      implicit none
C
      integer          NENOD
      double precision Shpf(NENOD,2), D(7,7)
      double precision EKt(2,NENOD,2,NENOD)
C
      integer          a, b, i, j
      double precision BBC(2,4)
C
      do a = 1, NENOD
C
         do j = 1, 4
            BBC(1,j) = Shpf(a,1)*D(1,j)
     &           +     Shpf(a,2)*D(4,j)
            BBC(2,j) = Shpf(a,1)*D(4,j)
     &           +     Shpf(a,2)*D(2,j)
         end do
C
         do b = 1, NENOD
C
            EKt(1,a,1,b) = EKt(1,a,1,b)
     &           +         BBC(1,1) * Shpf(b,1)
     &           +         BBC(1,4) * Shpf(b,2)
C
            EKt(1,a,2,b) = EKt(1,a,2,b)
     &           +         BBC(1,4) * Shpf(b,1)
     &           +         BBC(1,2) * Shpf(b,2)
C
            EKt(2,a,1,b) = EKt(2,a,1,b)
     &           +         BBC(2,1) * Shpf(b,1)
     &           +         BBC(2,4) * Shpf(b,2)
C
            EKt(2,a,2,b) = EKt(2,a,2,b)
     &           +         BBC(2,4) * Shpf(b,1)
     &           +         BBC(2,2) * Shpf(b,2)
C
         end do
C
      end do
C
      end
C
C
      subroutine accKM3D (NENOD,Shpf,D,EKt)
C
      implicit none
C
      integer          NENOD
      double precision Shpf(NENOD,3), D(7,7)
      double precision EKt(3,NENOD,3,NENOD)
C
      integer          a, b, i, j
      double precision BBC(3,6)
C
      do a = 1, NENOD
C
         do j = 1, 6
            BBC(1,j) = Shpf(a,1)*D(1,j)
     &           +     Shpf(a,2)*D(4,j)
     &           +     Shpf(a,3)*D(6,j)
            BBC(2,j) = Shpf(a,1)*D(4,j)
     &           +     Shpf(a,2)*D(2,j)
     &           +     Shpf(a,3)*D(5,j)
            BBC(3,j) = Shpf(a,1)*D(6,j)
     &           +     Shpf(a,2)*D(5,j)
     &           +     Shpf(a,3)*D(3,j)
         end do
C
         do b = 1, NENOD
C
            EKt(1,a,1,b) = EKt(1,a,1,b)
     &           +         BBC(1,1) * Shpf(b,1)
     &           +         BBC(1,4) * Shpf(b,2)
     &           +         BBC(1,6) * Shpf(b,3)
C
            EKt(1,a,2,b) = EKt(1,a,2,b)
     &           +         BBC(1,4) * Shpf(b,1)
     &           +         BBC(1,2) * Shpf(b,2)
     &           +         BBC(1,5) * Shpf(b,3)
C
            EKt(1,a,3,b) = EKt(1,a,3,b)
     &           +         BBC(1,6) * Shpf(b,1)
     &           +         BBC(1,5) * Shpf(b,2)
     &           +         BBC(1,3) * Shpf(b,3)
C
            EKt(2,a,1,b) = EKt(2,a,1,b)
     &           +         BBC(2,1) * Shpf(b,1)
     &           +         BBC(2,4) * Shpf(b,2)
     &           +         BBC(2,6) * Shpf(b,3)
C
            EKt(2,a,2,b) = EKt(2,a,2,b)
     &           +         BBC(2,4) * Shpf(b,1)
     &           +         BBC(2,2) * Shpf(b,2)
     &           +         BBC(2,5) * Shpf(b,3)
C
            EKt(2,a,3,b) = EKt(2,a,3,b)
     &           +         BBC(2,6) * Shpf(b,1)
     &           +         BBC(2,5) * Shpf(b,2)
     &           +         BBC(2,3) * Shpf(b,3)
C
            EKt(3,a,1,b) = EKt(3,a,1,b)
     &           +         BBC(3,1) * Shpf(b,1)
     &           +         BBC(3,4) * Shpf(b,2)
     &           +         BBC(3,6) * Shpf(b,3)
C
            EKt(3,a,2,b) = EKt(3,a,2,b)
     &           +         BBC(3,4) * Shpf(b,1)
     &           +         BBC(3,2) * Shpf(b,2)
     &           +         BBC(3,5) * Shpf(b,3)
C
            EKt(3,a,3,b) = EKt(3,a,3,b)
     &           +         BBC(3,6) * Shpf(b,1)
     &           +         BBC(3,5) * Shpf(b,2)
     &           +         BBC(3,3) * Shpf(b,3)
C
         end do
C
      end do
C
      end
