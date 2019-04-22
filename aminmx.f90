module aminmx_mod


contains

      subroutine aminmx(array,ifirst,ilast,istride,amin,amax, &
        indmin,indmax)
      implicit integer (i-n), real*8 (a-h,o-z)

!     compute max and min with indices

      dimension array(ilast)
!
      amin = array(ifirst)
      amax = array(ifirst)
      indmin = ifirst
      indmax = ifirst
      do i=ifirst+istride,ilast,istride
        if (array(i) .lt. amin) then
          amin = array(i)
          indmin = i
        end if
        if (array(i) .gt. amax) then
          amax = array(i)
          indmax = i
        end if
      end do
!
      return
      end
end module aminmx_mod