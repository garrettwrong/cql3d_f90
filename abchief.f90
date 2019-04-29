module abchief_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  use tdchief_mod, only : tdchief

  !---END USE
  !
  !

contains

  subroutine abchief
    use param_mod
    use comm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save

    !
    !      cputime=second()

    !.......................................................................
    !     2d code controlled by achief1, called at beginning of tdchief
    !.......................................................................
    call tdchief
    return
  end subroutine abchief

end module abchief_mod
