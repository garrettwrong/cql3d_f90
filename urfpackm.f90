module urfpackm_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains


end module urfpackm_mod

        subroutine pack(ip, nbits, iu, nw)
!
!       PACK - Compresses stored data
!       This version prepared for Absoft compiled data
!       (BobH, 990201).
!       SYNOPSIS
!           CALL PACK(ip, nbits, iu, nw)
!           PACK compresses stored integer data.  The following is
!           a list of valid arguments for this routine.
!           ip        On exit, vector of packed data in range -128:127,
!                         the range of integer*1 for Absoft fortran.
!           nbits     Number of rightmost bits of data in each
!                     partial word; must be  8
!           iu        Vector of integers in range 0:255 (beware),
!                         to be compressed.
!           nw        Number of integers to be compressed.
!
!           PACK takes the  8 rightmost bits
!           of several partial words and concatenates them
!           into INTEGER*1 words.
!           This implements a Cray builtin subroutine.
!           For alignment, it will be best to take nw to
!           be an integral multiple of 8, for 64-bit words.

        integer iu(nw),iu0,iu255
        integer*1 ip(nw),ip0,ip255

        do i=1,nw
           ip(i)=iu(i)-128
        enddo

        return
        end



!
!
        subroutine unpack(ip, nbits, iu, nw)
!       NAME
!           UNPACK - Compresses stored data
!      SYNOPSIS
!           CALL UNPACK(p, nbits, u, nw)
!           UNPACK compresses stored integer data.  The following is
!           a list of valid arguments for this routine.
!           ip        On exit, vector of packed data.
!           nbits     Number of rightmost bits of data in each
!                     partial word; must be  8
!           iu        Vector of partial words to be uncompressed.
!           nw        Number of partial words to be uncompressed.
!
!           This is a modification of a Cray Subroutine.
!           For alignment, it will be best to take nw to
!           be an integral multiple of 8, for 64-bit words.

        integer iu(nw)
        integer*1 :: ip(nw)

        do i=1,nw
           iu(i)=ip(i)+128
        enddo

        return
        end


!
!
        subroutine pack16(ip, nbits, iu, nw)
!
!       PACK - Compresses stored data
!       This version prepared for Absoft compiled data
!       (BobH, 990201).  (lf95, BobH, 050809)
!       SYNOPSIS
!           CALL PACK(ip, nbits, iu, nw)
!           PACK compresses stored integer data.  The following is
!           a list of valid arguments for this routine.
!           ip        On exit, vector of packed data in range -32768:32767,
!                         the range of integer*2 for lf95 fortran.
!           nbits     Number of rightmost bits of data in each
!                     partial word; must be  8
!           iu        Vector of integers in range 0:65535 (beware),
!                         to be compressed.
!           nw        Number of integers to be compressed.
!
!           PACK takes the  8 rightmost bits
!           of several partial words and concatenates them
!           into INTEGER*2 words.
!           This implements a Cray builtin subroutine.
!           For alignment, it will be best to take nw to
!           be an integral multiple of 8, for 64-bit words.

        integer iu(nw),iu0,iu255
        integer*2 ip(nw),ip0,ip255

        do i=1,nw
           ip(i)=iu(i)-32768
        enddo

        return
        end



!
!
        subroutine unpack16(ip, nbits, iu, nw)
!       NAME
!           UNPACK - Compresses stored data
!      SYNOPSIS
!           CALL UNPACK(p, nbits, u, nw)
!           UNPACK compresses stored integer data.  The following is
!           a list of valid arguments for this routine.
!           ip        On exit, vector of packed data.
!           nbits     Number of rightmost bits of data in each
!                     partial word; must be  8
!           iu        Vector of partial words to be uncompressed.
!           nw        Number of partial words to be uncompressed.
!
!           This is a modification of a Cray Subroutine.
!           For alignment, it will be best to take nw to
!           be an integral multiple of 8, for 64-bit words.

        integer iu(nw)
        integer*2 ip(nw)

        do i=1,nw
           iu(i)=ip(i)+32768
        enddo

        return
        end
