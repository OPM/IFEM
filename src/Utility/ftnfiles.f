      integer function openftnfile (fname)
      implicit none
      character(*) fname
      integer iUnit, iErr
      logical exist, named
      iUnit = 7
      iErr  = 0
      named = .true.
      exist = .false.
      do while (named)
         iUnit = iUnit + 1
         inquire(iUnit,EXIST=exist,NAMED=named)
      end do
      open(iUnit,FILE=fname,IOSTAT=iErr)
      if (iErr .eq. 0) then
         openftnfile = iUnit
      else if (iErr .gt. 0) then
         openftnfile = -iErr
      else
         openftnfile = 0
      end if
      end

      subroutine closeftnfile (iUnit)
      implicit none
      integer iUnit
      close(iUnit)
      end
