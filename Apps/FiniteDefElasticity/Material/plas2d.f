      subroutine plas2d (ipsw,iwr,iter,lfirst,ntm,pmat,detF,
     &                   ndf,Fn1,Fn,be,Epp,Epl,Sig,Cst,ierr)
C
C     PLAS2D: 2D wrapper on the 3D plastic material of FENRIS.
C
      implicit none
C
      integer  ipsw, iwr, iter, lfirst, ntm, ndf, ierr
      real*8   detF, Fn1(*), Fn(*), pmat(*), be(*)
      real*8   Epp, Epl(*), Sig(3), Cst(3,3)
      real*8   F3D(3,3), F3P(3,3), S3D(6), C3D(6,6)
C
      if (ndf .eq. 2) then
         F3D      = 0.0D0
         F3D(1,1) = Fn(1)
         F3D(2,1) = Fn(2)
         F3D(1,2) = Fn(3)
         F3D(2,2) = Fn(4)
         F3D(3,3) = 1.0D0
         F3P      = 0.0D0
         F3P(1,1) = Fn1(1)
         F3P(2,1) = Fn1(2)
         F3P(1,2) = Fn1(3)
         F3P(2,2) = Fn1(4)
         F3P(3,3) = 1.0D0
         call plas3d (ipsw,iwr,iter,lfirst,ntm,pmat,detF,
     &                F3P,F3D,be,Epp,Epl,S3D,C3D,ierr)
      else
         call plas3d (ipsw,iwr,iter,lfirst,ntm,pmat,detF,
     &                Fn,Fn1,be,Epp,Epl,S3D,C3D,ierr)
      end if
C
      Cst(1:2,1:2) = C3D(1:2,1:2)
      Cst(1:2,3)   = C3D(1:2,4)
      Cst(3,1:2)   = C3D(4,1:2)
      Cst(3,3)     = C3D(4,4)
      Sig(1:2)     = S3D(1:2)
      Sig(3)       = S3D(4)
C
      return
      end
