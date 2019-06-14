!     frname_decl.h
!
!.......................................................................
!     This file will hold all declarations of the namelist variables
!     for the FREYA module given in frname.h, type and dimensions.
!.......................................................................
      
      integer nprim,nimp
      common /numbrs_/ nprim,nimp
      
      integer naptr,ialign20,mf,nsourc,nbeams,npart,npskip,ne_tk
      common /nub_/ naptr,ialign20,mf,nsourc,nbeams,npart,npskip,ne_tk
      
      character*8 ashape,bshape
      common /nub_/ ashape(nap,kb),bshape(kb)
      
      !scalars:
      real(c_double) :: beamon,btime,ds_tk,fe_tk,ranseed,relnub
      real(c_double) :: fionx,hdepsmth,psicutoff
      common /nub_/ beamon,btime,ds_tk,fe_tk,ranseed,relnub
      common /nub_/ fionx,hdepsmth,psicutoff
      
      !vectors ------------------------------
      real(c_double) :: anglev(kb),angleh(kb),bleni(kb),blenp(kb)
      real(c_double) :: bcur(kb),bhofset(kb),bvofset(kb),bptor(kb),bheigh(kb),bwidth(kb)
      real(c_double) :: bhfoc(kb),bvfoc(kb)
      real(c_double) :: bhdiv(kb),bvdiv(kb)
      real(c_double) :: ebkev(kb)
      real(c_double) :: sfrac1(kb),rpivot(kb),zpivot(kb)
      real(c_double) :: timbplt(5)
      common /nub_/ anglev,angleh,bleni,blenp
      common /nub_/ bcur,bhofset,bvofset,bptor,bheigh,bwidth
      common /nub_/ bhfoc,bvfoc
      common /nub_/ bhdiv,bvdiv
      common /nub_/ ebkev
      common /nub_/ sfrac1,rpivot,zpivot
      common /nub_/ timbplt
      !vectors ------------------------------
      
      ! 2D arrays ---------------------------------
      real(c_double) :: aheigh(nap,kb)
      real(c_double) :: awidth(nap,kb)
      real(c_double) :: alen(nap,kb)
      real(c_double) :: fbcur(ke,kb)
      common /nub_/ aheigh
      common /nub_/ awidth
      common /nub_/ alen
      common /nub_/ fbcur
      ! 2D arrays ---------------------------------


!................................................................

      character*8    frplt,frmod,fr_gyro,beamplse,multiply
      common /nub2_/ frplt,frmod,fr_gyro,beamplse,multiply

      integer ibcur,ibcx,iborb,ibslow,inubpat,npat(2)
      integer itrapfi,itrapech
      common /nub2_/ ibcur,ibcx,iborb,ibslow,inubpat,npat
      common /nub2_/ itrapfi,itrapech
     
!     ONETWO DIVERGENCE
      integer        nfrplt,multiplyn
      common /nub2_/ nfrplt,multiplyn
      real(c_double) :: smooth,bmsprd,beampon,beampoff
      common /nub2_/    smooth,bmsprd,beampon,beampoff

!................................................................

      integer iexcit,ilorent,mstate
      integer kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      integer izstrp(kimp)
      common /nub3_/      iexcit,ilorent,mstate,kdene
      common /nub3_/ kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh
      common /nub3_/ ngl
      common /nub3_/ izstrp

!................................................................
      integer nouthx
      real(c_double) :: xdebug(20)
      common /io_/  xdebug,nouthx

!................................................................

!.......................................................................

      character*8 read_birth_pts
      character*128, dimension(24) ::  birth_pts_files
      integer nbirth_pts_files, nbirth_pts
      common /nubeam_list/ nbirth_pts_files, nbirth_pts
      common /nubeam_list/ read_birth_pts, birth_pts_files

!.......................................................................

      character*8 nameb
      common /ions_/ nameb
