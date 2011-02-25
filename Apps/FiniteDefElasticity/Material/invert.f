      subroutine invert (lda,n,a,ierr)
C
      implicit none
C
      integer, intent(in)    :: lda,n
      integer, intent(out)   :: ierr
      real*8 , intent(inout) :: a(lda,n)
C
      integer, allocatable :: ipiv(:)
      real*8 , allocatable :: work(:)
      real*8               :: rwork
C
      allocate(ipiv(n))
      call dgetrf (n,n,a,lda,ipiv,ierr)
      if (ierr /= 0) then
         deallocate(ipiv)
         return
      end if
C
      call dgetri (n,a,lda,ipiv,rwork,-1,ierr)
      allocate(work(int(rwork)))
      call dgetri (n,a,lda,ipiv,work,size(work),ierr)
      deallocate(ipiv,work)
C
      return
      end
