module dskin_idl_mod

  !---BEGIN USE

  use dskin_mod, only : dskin

  !---END USE



contains

      subroutine dskin_idl(argc,argv)  !Called by IDL
      implicit integer (i-n), real*8 (a-h,o-z)
      integer*4 argc,argv(*)
      integer*4 argv_len
      include 'dskin.h'
!
!     IDL setup to read distn function data
!
!     iya=80
!     jxa=80
!     lrza=15
!     ngena=1
!     x=dblarr(jxa)
!     y=dblarr(iya)
!     rovera=dblarr(lrza)
!     elecfld=dblarr(lrza)
!     bthr=dblarr(lrza)
!     btoru=dblarr(lrza)
!     bthr0=dblarr(lrza)
!     btor0=dblarr(lrza)
!     reden=dblarr(lrza,ngena)
!     temp=dblarr(lrza,ngena)
!     radmin=0d0
!     vnorm=0d0
!     vmaxdvt=0d0
!     eovedd=0d0
!     f=dblarr(iya,jxa,lrza,ngena)
!
!     I can't see any way to get the data from the dskin.h
!       arrays to the idl variables:  The idl variables
!       are set up with specific storage.  The dskin.h variables
!       are set up with different storage locations.  Somehow, the
!       dskin.h variables need to be copied to the idl
!       variables.
!
!     If we abandon the dskin.h common variables, and pass the needed
!       storage through subroutine arguments as in the standard
!       IDL/fortran examples, of course the problem will be
!       solved.
!
!     I'd just like to know how to set up new variables at the
!       the assigned locations.   Maybe with pointers.
!       Need to try pgf77 here.     WORKS!
!
!


      argv_len=loc(argc)         !Obtains number of arguments (argc)
                          !Because argc is passed by value,
                          !and fortran thinks is was a reference.

      open(unit=5,file='dskin_op',status='unknown')

!     The following line gives a segmentation fault
!     write(*,*) 'argc: ',argc
!     This is not too suprising, since from the fortran
!     point of view, argc is an address.  De-Referencing the address
!     (as occurs in the write statement) is probably illegal.
!     All we can do is ask what the address is: i.e., loc(argc).
!


      write(*,*) 'argv_len:', argv_len
      write(*,*) '(argv(j),j=1,argv_len): ',(argv(j),j=1,argv_len)
      write(5,*) 'argv_len:',argv_len
      write(5,*) '(argv(j),j=1,argv_len): ',(argv(j),j=1,argv_len)

!     NONE OF FOLLOWING USE OF VAL() IS LEGAL WITH PGF77.
!     ALSO CAN'T WRITE(*,*) %VAL WITH PGF77
!     (However %VAL works with g77).
!     Check can get value of initial
!      initial=val(argv(14))
!      energy=val(argv(15))
!      pitch=val(arg(16))
!      rho=val(argv(17))
!
!      write(*,*) initial,energy,pitch,rho
!      write(*,*) %val(argv(14)),%val(argv(15)),
!     +           %val(argv(16)),%val(argv(17))


      call dskin(%VAL(argv(14)),%VAL(argv(15)), &
                 %VAL(argv(16)),%VAL(argv(17)))

      write(*,*) rovera,radmin,vnorm
      write(5,*) rovera,radmin,vnorm

      close(unit=5)


      Setting pointers to addresses for data storage set up in IDL:
      fptr=argv(13)
      roveraptr=argv(3)
      radminptr=argv(9)

!*********************************************************************
!
!     GREAT:  The following works.  Can read f, rovera, radmin in IDL.
!
!*********************************************************************
      do i=1,iya
         do j=1,jxa
            do l=1,lrza
               f_idl(i,j,l,1)=f(i,j,l,1)
            enddo
         enddo
      enddo
      do l=1,lrza
         rovera_idl(l)=rovera(l)
      enddo
      radmin_idl=radmin


      return
      end

end module dskin_idl_mod
