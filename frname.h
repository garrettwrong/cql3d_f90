!     frname.h
!***********************************************************************
!     Namelist variables for the NFREYA fr.... routines
!***********************************************************************

!................................................................
!     Freya namelist variables.   It is not safe to put these
!     into main cql3d namelist files name.h, as there could
!     be duplication of variable names.  If you want to combine
!     with namelist in name.h, then check all the following variables.
!     Might have to change a couple of variable names.
!     frname.h is used included into frcomm.h, which is used by
!     freya.f,frinitl.f,frnfreya.f
! # YuP[2019-05-31] frhexdrv.f is removed - not used   
!................................................................
      namelist/frsetup/ timbplt,beamon,btime,nameb,relnub &
	,anglev,angleh,ashape,aheigh,awidth,bcur,bptor,blenp,bshape &
	,bleni,bheigh,bwidth,bhfoc,bvfoc,bhdiv,bvdiv,ebkev,fbcur,nbeams &
	,naptr,alen,bvofset,bhofset,nsourc,sfrac1,mf,npart,npskip &
	,rpivot,zpivot,ranseed,fionx,nbinject &
	,xdebug &
	,a1rf,a2rf,wrfe,wrfi,nrfzon,rfzone,ichmod,betalm &
	,beamplse,beampon,beampoff &
!     ONETWO DIVERGENCE
	,nimp,nprim,frmod,fr_gyro,smooth,multiply,multiplyn,bmsprd &
	,frplt,nfrplt,inubpat,npat &
	,ibcur,ibcx,ibslow,iborb,iyoka,ishot,itime &
	,itrapfi &
	,iexcit,ilorent,mstate,izstrp,kdene &
	,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl,nouthx &
	,hdepsmth &
	,birth_pts_files,nbirth_pts_files,nbirth_pts,read_birth_pts &
	,ne_tk,ds_tk,fe_tk &
!     Remove NBI source at psi outside of psicutoff:
	,psicutoff


!Removed namelist, not used, BH070308:
! tfusbb,iddcal,fdbeam,nbinject,rfmode,rfon,rftime,rfpow,
! idamp,isrc,zsrc,fpsrc,phisrc,
! freq,rfrad1,rfrad2,wrfo,wrfx,rnormin,njqin,qine,qini,
! a1rf,a2rf,wrfe,wrfi,nrfzon,rfzone,ichmod,betalm
! relrf,nprf,iside,xkpar,nhigh,ykperp,navg
! nampel,pelrad,vpel,nbgpel,timpel
! ifus,iaslow,wtifus
! itrapech,
! nshell,nraysh,thetsh,phish,praysh
! wdelt,wgam,wohm,nqrad,qradr,qradin,refrad
! rnp,irfcur,ifb,rfcur,lmode,ifbprof,alphaf
