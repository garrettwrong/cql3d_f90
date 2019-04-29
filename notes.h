!
!
!        Notes File.  (Try starting this as of June 21, 1994).
!
!.......................................................................
!     Things to watch for in old cqlinput NAMELIST FILES:
!     (BobH, 990702)
!.......................................................................
!
!
!      -To make namelist files maximally compatible
!       across platforms (Joe Freeman experience),
!       begin and end all namelist sections with "$",
!       (rather than "&" or "/".
!      -Replace all "-quotes with '-quotes.
!
!      -On CRAY J90/f90 (only), need to remove non-namelist data
!       such as "**" between lines, header data before
!       first namelist section, and comments following ";".
!      -Fix eegy(,,) ==> eegy(,,,)
!      -Fix asor(,)  ==> asor(,,)
!      -Add noplots="enabled" to first invocation of "setup"
!       namelist, until all the dead-end references to GRAFLIB
!       routines are completely replaced by PGPLOT
!      -Check nconteq,izeff and multiply are character*8
!       (integer values now put in through nconteqn,multiplyn).
!      -Check fpld(,) is type real
!      -Careful that nrfstep1(3),pwrscale(3),wdscale(3)
!       not set if nmodsa.le.2,
!       nonrf(2),noffrf(2) not set if ngena=1.
!      -Set bsign=-1. in eqsetup namelist, if the toroidal field
!       in the eqdsk is negative.
!      -Need to set gamaset to nonzero value for ion simulation.
!       (Upgrade this?).
!
!
!.......................................................................
!     Modifications to 32-bit code to give code running on Cray J90
!     64-bit environment:                       990701 (BobH)
!.......................................................................
!
!      -Changes are indicated in the code by comment lines.
!       c_cray, as compared to c_pc. In further detail:
!      -Change machinea=2==>1 in param.h, for cray
!      -Reverse use of dgbtrf,dgbtrs to sgbtrf,sgbtrs (can use
!       replace_sgbtrf script).
!      -malloc,free ==> hpalloc,hpdeallc
!      -Reverse use of drand(iflag) and ranf()
!      -Reverse use of call system and call ishell
!      -Need to set environmenatal variables
!       PGPLOT_DIR=(pgplot directory) ,PGPLOT_DEV=/XWINDOW, e.g.
!      -makefile:remove references to r8lsode and urfpackm.
!      -cqlinput namelist file cannot have header of characters
!       between the namelist sections (such as "**").
!      -use netcdf.inc compatible with installed netCDF.
!      -note on netcdf: the NCDOUBLE designation for storage
!       of netcdf variables is netCDF based, rather than on
!       the architecture of the computer.
!       Thus, real(c_double) is NCDOUBLE, regardless.
!
!
!.......................................................................
!
!  f_code, indicating distribution functions calculated in the
!     code, is equal to f*vnorm**3,
!     where  f=distribution function in particles per cm**3, and
!              per velocity**3 (cgs).  For the relatavistic case
!              velocity is actually (momentum/rest_mass).
!            vnorm=maximum velocity of the grid (maximum
!              momentum/rest_mass, for the relativistic case).
!
!
!
!  c..................................................................
!  c     generate radial (rho) mesh. rya will be the normalized mesh.
!  c     rrz will be the intermediate (mid) mesh points. Later rz
!  c     will be defined to be the non-normalized actual radial
!  c     toroidal flux mesh.
!  c..................................................................
!  c
!  c
!  c.......................................................................
!  c     determine array rovera
!  c.......................................................................
!
!        do 30 ll=1,lrzmax
!          rovera(ll)=rya(ll)
!          if (0.lt.rovera(ll) .and. rovera(ll).lt.1.e-8)
!       +    rovera(ll)=1.e-8
!   30   continue
!
!  c..................................................................
!  c     fill in input arrays between rya=0. and rya=1.
!  c     This determines density, temperature, source, etc profiles as
!  c     functions of the normalized radial coordinate rho. If density
!  c     and temperatures are to be specified as functions of bbpsi
!  c     poloidal flux coordinate these arrays will be overwritten
!  c     in subroutine eqcoord
!  c     If eqsource="tsc" (running with tsc):
!  c     iprone=iprote=iproti="disabled" (set in tdtscinp).
!  c..................................................................
!  c
!  c
!  c..................................................................
!  c..................................................................
!  c       cosz(1:iy,1:lzmax,1:lrzmax)  is cosine of pitch angle
!  c         corresponding to y(i,lmdlpln_), as given by cnsts of motion.
!  c..................................................................
!          if (ilcqlp .eq. 0) then
!            do 140 i=1,iiyh
!              cosz(i,l,lr_)=cvmgt((1.-bbpsi(l,lr_)*(sinn(i,lmdpln_)**2)),
!       1        1.e-90,1.-bbpsi(l,lr_)*(sinn(i,lmdpln_)**2) .gt. 1.e-11)
!              cosz(i,l,lr_)=sqrt(cosz(i,l,lr_))
!   140      continue
!          else
!            do 141 i=1,iiyh
!              cosz(i,l,lr_)=coss(i,indxls(l))
!   141      continue
!          endif
!  c..................................................................
!  c     imax(l,lr_) is the highest i such that a particle with a midplane
!  c     pitch angle y(i,lmdpln_) attains and passes z(l,lr_) before turning.
!  c     exception: imax(lz,lr_)=itl.
!  c..................................................................
!          do 150 i=1,iiyh
!            iii=iiy+1-i
!   150    cosz(iii,l,lr_)=-cosz(i,l,lr_)
!
!          imax(l,lr_)=1
!          if (ilcqlp .eq. 0) then
!            do 160 i=2,iyh_(lmdpln_)
!              if (cosz(i,l,lr_) .le. 1.e-10) go to 170
!   160      continue
!            imax(l,lr_)=iyh_(lmdpln_)
!            go to 180
!   170      imax(l,lr_)=i-1
!   180      continue
!          else
!            imax(l,lr_)=iiyh
!          endif
!
!  c..................................................................
!  c     sinz and tanz are defined analogously to cosz (above)
!  c..................................................................
!  c
!  c
!  c
!  c
!  c
!  c..................................................................
!  c     lmax(i,lr_) is the highest l such theta particle with midplane
!  c     pitch angle y(i,l_) attains and passes z(l,lr_). exception lmax(itl
!  c..................................................................
!
!  c
!  c
!  c
!  c
!  c
!  c
!  c..................................................................
!  c     zboun(i,lr_) is the point along a orbit at which a trapped
!  c     particle bounces. Passing particles are assigned zmax(lr_).
!  c     Particles at the passed trapped boundary are temporarily
!  c     assigned zstar. See subroutine baviorbt.
!  c..................................................................
!  c
!  c
!  c
!
!
!
!  There are (at least) 2 meshes along the magnetic field line:
!  (1)Subroutine eqfndpsi works with 2D arrays with arguments
!     l=1:lorbit(lr_),lr_=1,lrzmax.  lorbit(lr_) takes values
!     up to lfielda (=250, presently), and gives a fine mesh of
!     points along a flux surface:
!        es(1:lorbit(lr_),lr_) is distance along the field line,
!          measured from the outer minimum B-field point (i.e., the
!          equatorial plane.
!        solr(1:lorbit(lr_),lr_) corresonding major radii
!        solz(             ,   )              verical height
!        eqbpol(               )   B_poloidal (T)
!        bpsi(                 )   B(es)/B(0.)
!        thtpol(               )   atan(solz/(rmag-solr))
!        eqdell(               )   delta(poloidal distance between
!                                  mesh points).
!  (2)Subroutine micxiniz works with a courser mesh along the field
!     line (requiring less storage and computations), but derived
!     from the above mesh:
!     z(1:lz,1:lrzmax) distance along field line,
!                      where  namelist lz.le.lza=80 presently.
!                      Spacing determined by tfacz.
!     pol(1:lz,lrzmax) from interpolation of thtpol at z( , )
!                           (psimodel.eq."spline")
!     bbpsi(           )                     bpsi at z( , )
!     solrz(  ,      )                       solr at z( , )
!
