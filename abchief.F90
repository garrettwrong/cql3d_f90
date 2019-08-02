module abchief_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use aclear_mod, only : aclear
  use tdchief_mod, only : tdchief

  !---END USE
  !
  !

contains

  subroutine abchief(nml_file)
    use param_mod
    use cqlcomm_mod
    implicit integer (i-n), real(c_double) (a-h,o-z)
    save

    character(len=*), intent(in), optional :: nml_file

    !
    !      cputime=second()

    !.......................................................................
    !     2d code controlled by achief1, called at beginning of tdchief
    !.......................................................................
    call tdchief(nml_file)
    call aclear ! reset for STEP mode
    return
  end subroutine abchief

end module abchief_mod
