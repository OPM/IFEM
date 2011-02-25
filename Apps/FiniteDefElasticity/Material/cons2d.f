      subroutine cons2d (ipsw,iter,iwr,lfirst,mTYP,mVER,
     &                   detF,F,Fp,pmat,HV,U,Sig,Cst,ierr)
C
C     CONS2D: 2D wrapper on the 3D material of FENRIS.
C
      implicit none
C
      integer  ipsw, iter, iwr, lfirst, mTYP, mVER, ierr
      real*8   detF, F(4), Fp(4), pmat(*), HV, U, Sig(3), Cst(3,3)
      real*8   F3D(3,3), F3P(3,3), S3D(6), C3D(6,6)
C
      F3D      = 0.0D0
      F3D(1,1) = F(1)
      F3D(2,1) = F(2)
      F3D(1,2) = F(3)
      F3D(2,2) = F(4)
      F3D(3,3) = 1.0D0
      F3P      = 0.0D0
      F3P(1,1) = Fp(1)
      F3P(2,1) = Fp(2)
      F3P(1,2) = Fp(3)
      F3P(2,2) = Fp(4)
      F3P(3,3) = 1.0D0
C
      call cons3d (ipsw,iter,iwr,lfirst,mTYP,mVER,
     &             detF,F3D,F3P,pmat,HV,U,S3D,C3D,ierr)
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
