!     frname_decl.h
!
!.......................................................................
!     This file will hold all declarations of the namelist variables
!     for the FREYA module given in frname.h, type and dimensions.
!.......................................................................

      common /numbrs_/   nprim,nimp

      character*8 ashape,bshape
      common /nub_/      anglev(kb),angleh(kb),naptr,ialign20 &
        ,ashape(nap,kb),aheigh(nap,kb),awidth(nap,kb),bcur(kb) &
        ,alen(nap,kb),bleni(kb),blenp(kb) &
        ,bhofset(kb),bvofset(kb),bptor(kb) &
        ,bshape(kb),bheigh(kb),bwidth(kb) &
        ,bhfoc(kb),bvfoc(kb) &
        ,bhdiv(kb),bvdiv(kb) &
        ,beamon,btime,ebkev(kb) &
        ,fbcur(ke,kb) &
        ,mf,nbeams,npart,npskip,ne_tk &
        ,nsourc,ds_tk,fe_tk,ranseed,relnub &
        ,sfrac1(kb),rpivot(kb),zpivot(kb) &
        ,timbplt(5) &
        ,fionx,hdepsmth &
!     Remove NBI source at psi outside of psicutoff:
        ,psicutoff


!................................................................

      character*8 frplt,frmod,fr_gyro,beamplse,multiply

      common /nub2_/ ibcur,ibcx,iborb,ibslow,inubpat,npat(2),itrapfi,itrapech &
!     ONETWO DIVERGENCE
        ,smooth,multiply,bmsprd,frmod,fr_gyro,beamplse,frplt,nfrplt &
        ,multiplyn,beampon,beampoff

!................................................................

      common /nub3_/      iexcit,ilorent,mstate,kdene &
        ,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh &
        ,ngl &
        ,izstrp(kimp)

!................................................................

      common /io_/  xdebug(20),nouthx

!................................................................

!.......................................................................

      character*8 read_birth_pts
      character*128, dimension(24) ::  birth_pts_files
      integer nbirth_pts_files, nbirth_pts
      common /nubeam_list/ nbirth_pts_files, nbirth_pts, &
        read_birth_pts, birth_pts_files

!.......................................................................

      character*8 nameb
      common /ions_/ nameb
