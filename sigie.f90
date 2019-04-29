module sigie_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      real(c_double) function sigie(en)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
      data a1/-3.17385e+01/
      data a2/1.143818e+01/
      data a3/-3.8339981/
      data a4/7.046692e-01/
      data a5/-7.431486e-02/
      data a6/4.153749e-03/
      data a7/-9.486967e-05/
!
!mnt  electron impact ionization. en in kev.
!
      zte=1.e+03*2./3.*en
!990131      zx=alog(zte)
      zx=log(zte)
      zlog=(((((a7*zx+a6)*zx+a5)*zx+a4)*zx+a3)*zx+a2)*zx+a1
      sigie=exp(zlog)
      return
      end
end module sigie_mod
