module abchief_mod

  !---BEGIN USE

  use tdchief_mod, only : tdchief

  !---END USE
  !
  !

contains

  subroutine abchief
    use param_mod
    use comm_mod
    implicit integer (i-n), real*8 (a-h,o-z)
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
