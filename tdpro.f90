module tdpro_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      real(c_double) function tdpro(psi,rplasm,acoef)
      implicit integer (i-n), real(c_double) (a-h,o-z)
!..................................................................
!     Calculate ASDEX YAG1 type profiles.
!     acoef(i),i=1,4  must be specified.
!..................................................................

      dimension acoef(4)

      x = rplasm*psi
      x2=x*x
      arg = acoef(4)
      do 10  i=1,3
        arg = arg*x2 + acoef(4-i)
 10   continue
      tdpro = exp(arg)
      return
      end
end module tdpro_mod
