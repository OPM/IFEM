C $Id: accKMx.f,v 1.2 2010-12-20 17:01:21 kmo Exp $
C=======================================================================
C!
C! \file accKMx.f
C!
C! \date Dec 18 2010
C!
C! \author Knut Morten Okstad / SINTEF
C!
C! \brief Material stiffness accumulation routines for B-bar methods.
C!
C=======================================================================
C
      subroutine accKMx2D (axiSymm,NENOD,Shpr,Shpf,Shpbar,D,EKt)
C
      implicit none
C
      integer          axiSymm, NENOD
      double precision Shpr(NENOD), Shpf(NENOD,2)
      double precision Shpbar(NENOD,2), D(7,7)
      double precision EKt(2,NENOD,2,NENOD)
C
      integer          a, b, i, j
      double precision BBC(2,7)
C
      do a = 1, NENOD
C
         do j = 1, 7
            BBC(1,j) = Shpf(a,1)*D(1,j)
     &           +     Shpf(a,2)*D(4,j)
     &           +   Shpbar(a,1)*D(7,j)
            BBC(2,j) = Shpf(a,1)*D(4,j)
     &           +     Shpf(a,2)*D(2,j)
     &           +   Shpbar(a,2)*D(7,j)
         end do
         if (axiSymm .eq. 1) then
            BBC(1,:) = BBC(1,:) + Shpr(a)*D(3,:)
         end if
C
         do b = 1, NENOD
C
            EKt(1,a,1,b) = EKt(1,a,1,b)
     &           +         BBC(1,1) *   Shpf(b,1)
     &           +         BBC(1,4) *   Shpf(b,2)
     &           +         BBC(1,7) * Shpbar(b,1)
C
            EKt(1,a,2,b) = EKt(1,a,2,b)
     &           +         BBC(1,4) *   Shpf(b,1)
     &           +         BBC(1,2) *   Shpf(b,2)
     &           +         BBC(1,7) * Shpbar(b,2)
C
            EKt(2,a,1,b) = EKt(2,a,1,b)
     &           +         BBC(2,1) *   Shpf(b,1)
     &           +         BBC(2,4) *   Shpf(b,2)
     &           +         BBC(2,7) * Shpbar(b,1)
C
            EKt(2,a,2,b) = EKt(2,a,2,b)
     &           +         BBC(2,4) *   Shpf(b,1)
     &           +         BBC(2,2) *   Shpf(b,2)
     &           +         BBC(2,7) * Shpbar(b,2)
C
            if (axiSymm .eq. 1) then
               EKt(1,a,1,b) = EKt(1,a,1,b) + BBC(1,3)*Shpr(b)
               EKt(2,a,1,b) = EKt(2,a,1,b) + BBC(2,3)*Shpr(b)
            end if
         end do
C
      end do
C
      end
C
C
      subroutine accKMx3D (NENOD,Shpf,Shpbar,D,EKt)
C
      implicit none
C
      integer          NENOD
      double precision Shpf(NENOD,3), Shpbar(NENOD,3), D(7,7)
      double precision EKt(3,NENOD,3,NENOD)
C
      integer          a, b, i, j
      double precision BBC(3,7)
C
      do a = 1, NENOD
C
         do j = 1, 7
            BBC(1,j) = Shpf(a,1)*D(1,j)
     &           +     Shpf(a,2)*D(4,j)
     &           +     Shpf(a,3)*D(6,j)
     &           +   Shpbar(a,1)*D(7,j)
            BBC(2,j) = Shpf(a,1)*D(4,j)
     &           +     Shpf(a,2)*D(2,j)
     &           +     Shpf(a,3)*D(5,j)
     &           +   Shpbar(a,2)*D(7,j)
            BBC(3,j) = Shpf(a,1)*D(6,j)
     &           +     Shpf(a,2)*D(5,j)
     &           +     Shpf(a,3)*D(3,j)
     &           +   Shpbar(a,3)*D(7,j)
         end do
C
         do b = 1, NENOD
C
            EKt(1,a,1,b) = EKt(1,a,1,b)
     &           +         BBC(1,1) *   Shpf(b,1)
     &           +         BBC(1,4) *   Shpf(b,2)
     &           +         BBC(1,6) *   Shpf(b,3)
     &           +         BBC(1,7) * Shpbar(b,1)
C
            EKt(1,a,2,b) = EKt(1,a,2,b)
     &           +         BBC(1,4) *   Shpf(b,1)
     &           +         BBC(1,2) *   Shpf(b,2)
     &           +         BBC(1,5) *   Shpf(b,3)
     &           +         BBC(1,7) * Shpbar(b,2)
C
            EKt(1,a,3,b) = EKt(1,a,3,b)
     &           +         BBC(1,6) *   Shpf(b,1)
     &           +         BBC(1,5) *   Shpf(b,2)
     &           +         BBC(1,3) *   Shpf(b,3)
     &           +         BBC(1,7) * Shpbar(b,3)
C
            EKt(2,a,1,b) = EKt(2,a,1,b)
     &           +         BBC(2,1) *   Shpf(b,1)
     &           +         BBC(2,4) *   Shpf(b,2)
     &           +         BBC(2,6) *   Shpf(b,3)
     &           +         BBC(2,7) * Shpbar(b,1)
C
            EKt(2,a,2,b) = EKt(2,a,2,b)
     &           +         BBC(2,4) *   Shpf(b,1)
     &           +         BBC(2,2) *   Shpf(b,2)
     &           +         BBC(2,5) *   Shpf(b,3)
     &           +         BBC(2,7) * Shpbar(b,2)
C
            EKt(2,a,3,b) = EKt(2,a,3,b)
     &           +         BBC(2,6) *   Shpf(b,1)
     &           +         BBC(2,5) *   Shpf(b,2)
     &           +         BBC(2,3) *   Shpf(b,3)
     &           +         BBC(2,7) * Shpbar(b,3)
C
            EKt(3,a,1,b) = EKt(3,a,1,b)
     &           +         BBC(3,1) *   Shpf(b,1)
     &           +         BBC(3,4) *   Shpf(b,2)
     &           +         BBC(3,6) *   Shpf(b,3)
     &           +         BBC(3,7) * Shpbar(b,1)
C
            EKt(3,a,2,b) = EKt(3,a,2,b)
     &           +         BBC(3,4) *   Shpf(b,1)
     &           +         BBC(3,2) *   Shpf(b,2)
     &           +         BBC(3,5) *   Shpf(b,3)
     &           +         BBC(3,7) * Shpbar(b,2)
C
            EKt(3,a,3,b) = EKt(3,a,3,b)
     &           +         BBC(3,6) *   Shpf(b,1)
     &           +         BBC(3,5) *   Shpf(b,2)
     &           +         BBC(3,3) *   Shpf(b,3)
     &           +         BBC(3,7) * Shpbar(b,3)
C
         end do
C
      end do
C
      end
