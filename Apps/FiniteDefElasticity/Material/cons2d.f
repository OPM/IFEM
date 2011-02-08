      subroutine cons2d (mTYP,mVER,detF,F,pmat,U,Sig,Cst,ipsw,iwr,ierr)
C
C     CONS2D: 2D wrapper on the 3D material of FENRIS.
C
      implicit none
C
      integer  mTYP, mVER, ipsw, iwr, ierr
      real*8   detF, U, F(4), pmat(*), Sig(3), Cst(3,3)
      real*8   F3D(3,3), S3D(6), C3D(6,6)
C
      F3D      = 0.0D0
      F3D(1,1) = F(1)
      F3D(2,1) = F(2)
      F3D(1,2) = F(3)
      F3D(2,2) = F(4)
      F3D(3,3) = 1.0D0

      call cons3d (mTYP,mVER,detF,F3D,pmat,U,S3D,C3D,ipsw,iwr,ierr)

      Cst(1:2,1:2) = C3D(1:2,1:2)
      Cst(1:2,3)   = C3D(1:2,4)
      Cst(3,1:2)   = C3D(4,1:2)
      Cst(3,3)     = C3D(4,4)
      Sig(1:2)     = S3D(1:2)
      Sig(3)       = S3D(4)
C
      return
      end

