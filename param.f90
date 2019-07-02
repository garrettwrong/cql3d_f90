  !     param.f90
  !**********************************************************************
  !**********************************************************************

  !     PARAMETERS CHOSEN BY USER FOLLOW

  !..................................................................
  !     version is simply a setup0%mnemonic printed out to indicate
  !     the version of the source being run.  Further below, the non-
  !     user set parameter precursr indicates from which author-released
  !     version this source is descended.
  !
  !_cray machinea is =1 for 64-bit integers (e.g., Cray), and
  !_pc               =2 for 32-bit integers (e.g., PCs).
  !BH081218:  Present usage, machinea=2 works with 32- and 64-bit machines
  !
  !     ngena is the maximum number of general (time advanced) species.
  !
  !     nmaxa is the maximum number of background species which are
  !     not time advanced and which are assumed to be Maxwellian at a
  !     fixed temperature and density. These species exist only
  !     to contribute to the collision integral.
  !
  !     mx is the order of Legendre terms used in Rosenbluth
  !     expansions of the collision integral.
  !BH180901/YuP: Following restriction on mx has been removed.
  !     NOTE:   If relativ="enabled" mx => 2,
  !             BUT, 2*jx*(mx+3)*(mx+3) <= (iyp1+1)*(jxp1+1) to
  !                  avoid overwrite by tamt1 and tamt2.
  !
  !     lza is the number of poloidal mesh points for purposes of
  !     computing bounce averages (number of "z" mesh points).
  !     (if setup0%cqlpmod="enabled", lza=lsa=number of mesh points
  !     along the magnetic field line)
  !
  !     jx is the number of speed (momentum) mesh points.
  !
  !     iy is the number of pitch angle mesh points.
  !
  !     noncha is the number of elements in arrays carrying time plot
  !     information.
  !
  !     nbctimea is max number of points in arrays supplying time
  !     dependent profile information.
  !
  !     ndtr1a is maximum number of time step intervals dtr1().
  !
  !     nplota is maximum number of plot times for 2d and 3d plots.
  !
  !     nefitera is the maximum number of iterations permitted for
  !     the electric field per time step (to obtain a target current).
  !
  !..................................................................

  !     PARAMETERS CHOSEN BY USER FOLLOW

module param_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE


  character(len=64), parameter :: version="cql3d_git_190531.0"
  character(len=64), parameter :: precursr="cql3d_git_190309.0_GW"
  parameter(machinea=2)
  !BH081218:  Present usage, machinea=2 works with 32- and 64-bit machines
  parameter(ngena=4)
  !BH120727:  Moved nmaxa up to 8, for SWIM Plasma State Abridged species
  !BH120727:  list + elec.
  parameter(nmaxa=8)
  !        parameter(mxa=3)     YuP-101216: Not used anymore
  !        parameter(jxa=300)   YuP-101216: Not used anymore
  !        parameter(iya=160)   YuP-101216: Not used anymore
  !        parameter(noncha=2000) YuP-101221: Not used anymore
  parameter(nbctimea=101)
  parameter(ndtr1a=10)
  parameter(nplota=10)
  parameter(nsavea=10)   !Max number of distn time steps saves
  !as specified through n.eq.nsave()
  parameter(nefitera=10)
  parameter(ntavga=16)

  !*******************************************************************
  !     BEGIN parameters for SOURCE routines
  !*******************************************************************

  !..................................................................
  !     nsoa is the maximum number of sources per species allowed.
  !..................................................................

  parameter (nsoa=3)
  !..................................................................
  !     jfla is maximum length of grid for the reduced parallel
  !     distribution function used for the Besedin-Pankratov
  !     knock-on source. Use odd number.
  !..................................................................
  !        parameter (jfla=jxa)    YuP-101216: Not used anymore


  !*******************************************************************
  !     Parameter for knock-on source soureko routine
  !*******************************************************************

  !..................................................................
  !
  !     i0param is a table length in pitch angle, for translation
  !     of sin(theta) a given poloidal position into the corresponding
  !     equatorial pitch angle.
  !
  !..................................................................

  !BH180720      parameter (i0param=1001)
  parameter (i0param=2001)


  !*******************************************************************
  !     BEGIN parameters for 3-D (td...) routines
  !*******************************************************************

  !..................................................................
  !     lrza is the number of radial mesh points
  !
  !     lsa is the number of parallel (spatial) mesh points
  !     lsa should be .ge.lrza [BH081106:  perhaps this the usual
  !       case before big memory and bigger lrza cases.   Its
  !       contravention is causing an out-of-bounds ref for
  !       setup0%cqlpmod.ne.enabled.   Could check code out further
  !       to see how essential it is.]
  !
  !     Note: Limitations on relative sizes of lrza, lsa, and lza!
  !
  !     nbanda is the maximum bandwidth of the parallel transport equation
  !
  !     msxr is the order of the highest polynomial used in the Legendre
  !     expansion for computing SXR spectra (It is set to mx in aindflt.f).
  !
  !     nena is the maximum number of energies in the SXR spectra calc.
  !
  !     nva is the maximum number of SXR sightlines.
  !
  !     ntrmdera is the maximum of node for computing df/ds
  !
  !     njenea is the number of nodes for spline profiles
  !     (ryain, tein, tiin, enein, elecin)
  !..................................................................

  parameter(lrza=128)
  parameter(lsa=max(128,lrza), lsa1=lsa+1)
  !     lza should be .ge. lsa. if setup0%cqlpmod, lz is set equal to setup0%ls
  parameter(lza=lsa)
  !     lrorsa should be equal to max of lrza and lsa
  parameter(lrorsa=lsa)

  parameter(nbanda=5)

  parameter(nena=60,nva=100)
  parameter(ntrmdera=4)
  parameter(njenea=256)

  !*******************************************************************
  !     BEGIN parameters for the NON-CIRCULAR CROSS SECTION (eq..)
  !     routines.
  !*******************************************************************

  !..................................................................
  !     nnra (nnza) is the number of R (Z) meshpoints used in
  !     MHD equilibrium calculation
  !
  !     lfielda is the maximum number of poloidal points allowed
  !     in the orbit calculation. (250 is a good number)
  !..................................................................

  !      parameter(nnra=201)
  !      parameter(nnza=201)
  parameter(nnra=257)
  parameter(nnza=257)
  parameter(lfielda=250)
  !$$$P-> Grid size for storing the values
  !$$$c     of equilibrium field (Beqr,Beqz,Beqphi, etc.)
  !$$$c     on req(1:nreqa),zeq(1:nzeqa) grid
  !$$$c     for finite-orbit-width (FOW) tracing.
  !$$$      parameter(nreqa=128, nzeqa=196)
  !$$$c-YuP-> Could be set to (nnra,nnza) ?


  !*******************************************************************
  !     BEGIN parameters for the NFREYA (beam deposition) (fr..) code.
  !*******************************************************************

  !..................................................................
  !     Parameters for the "fr" routines (NFREYA code) as it exists in
  !     the transport code ONETWO. Parameters and common blocks
  !     appearing will be defined as in ONETWO.
  !     Where possible the coding has been lifted from ONETWO.
  !
  !..................................................................

  !        parameter(maxp=1500)  YuP-101211: maxp is not used anymore.

  !*******************************************************************
  !     BEGIN parameters for LH, FW,  and ECH (urf..) routines
  !*******************************************************************

  !..................................................................
  !     nraya is the maximum allowable number of LH rays
  !
  !     nrayelta is the maximum allowable number of ray elements
  !     per ray. A ray element is a small subsection of a ray
  !     characterized by a position, a length, to which the
  !     ray tracing code assigns a number of wave characteristics.
  !
  !BH060314  nmodsa is the maximum number of either number of wave modes,
  !BH060314  or number of harmonics for a single wave mode (presently
  !BH060314  there is some hardwiring to values .le.3).(and some hardwiring
  !BH060314  to values .ge. 3. E.G., nrfstep1,powrf,powurf,powurfi).
  !BH060314  Need to check this out before changing nmodsa.ne.3!
  !BH060314  (BobH, 990702).  This has been upgraded as follows.

  !     nmodsa is maximum of the sum over{ wave types (mrf) times
  !     harmonics for each wave type}.  (BobH, 060314).
  !
  !     nharma is maximum harmonic number for cyclotron interactions.
  !..................................................................

  !      parameter (nraya=1)
  !      parameter (nrayelta=1)
  !     YuP 101122: nraya and nrayelta are not used anymore.
  !     YuP 101122: Instead, nrayn and nrayelts are determined
  !     YuP 101122: in urfsetup by reading rays\' data files.

  parameter (nmodsa=155)  !Suggest not using .le.3, unless check
  !cqlinput that some vestigial inputs
  !are not set for index larger than nmodsa.
  !        parameter (nharma=1,nharm2a=nharma+2)
  !     YuP 101208: nharma is not used anymore.
  !
  !..................................................................
  !     rdcmod related
  !..................................................................
  parameter(nrdca=10)   !Max number of diffusion coeff files,
  !for rdcmod="format1"
  !
  !..................................................................
  !     NPA related:
  parameter(npaproca=5)
  !..................................................................


  !*******************************************************************
  !     Parameters for sigma-v routines
  !*******************************************************************

  !..................................................................
  !     mtaba is a table length for passage of sigma values.
  !
  !     mmsv is the order of the highest polynomial used in the Legendre
  !     expansion for computing sigma-v (set to mx in aindflt.f)
  !
  !..................................................................

  !BH150620  Found that energy increment delegy is much too large
  !BH150620  in a case with general species e,d,t,alpha.
  !BH150620  Suspect table controlled by electrons.
  !BH150620  Quick fix may be to increase mtaba by !mp/me
  !HB150620  With mtaba=1000000, still only get 31 nonzero entries in svtab()
  !BH150620  Needs further investigation/coding!
  !BH150620      parameter (mtaba=1000)
  parameter (mtaba=1000000)

  !        parameter (mmsva=mxa)    YuP-101216: Not used anymore


  !     END OF PARAMETER SET WHICH IS NORMALLY USER CHOSEN FOLLOW
  !***********************************************************************



  !..................................................................
  !     Parameters defined just below should not be altered
  !     (with the possible exception of negyrga (used in plots)).
  !..................................................................

  !BH160911      parameter(negyrga=3)
  parameter(negyrga=4)
  parameter(mbeta=10)
  !        parameter(mxp1a=mxa+1)  YuP-101216: Not used anymore
  !        parameter(ngenpa=ngena+1)  Not used anymore
  !        parameter(iyjxa=iya*jxa) YuP-101216: Not used anymore
  !        parameter(iyjx2a=(iya+2)*(jxa+2)) YuP-101216: Not used anymore
  !        parameter(iyjxnga=iyjxa*ngena) YuP-101216: Not used anymore
  !        parameter(jxp1a=jxa+1) YuP-101216: Not used anymore
  !        parameter(iyp1a=iya+1) YuP-101216: Not used anymore
  parameter(ntotala=ngena+nmaxa)
  parameter(ift07a=01,ift21a=01)

  !        parameter(ipxya=iya)    YuP-101216: Not used anymore
  !        parameter(jpxya=jxa+1)  YuP-101216: Not used anymore
  !        parameter(iyjxua=iya*(jxa+1)) YuP-101216: Not used anymore
  !        parameter(iyujxa=(iya+1)*jxa) YuP-101216: Not used anymore
  !        parameter(miyjxa=6*iyjxa)     YuP-101216: Not used ?
  parameter(incza=301,inczpa=incza+1)

  integer tlfld1a
  parameter(tlfld1a=3*lfielda+1)
  !     parameter(nconteqa=nnra)
  parameter(nnrab2=(nnra+1)/2)
  parameter(nconteqa=nnrab2)
  parameter(nrz3p1a=3*(nnza+nnra)+1)


  parameter(ki=nnra,kix2=2,kj=nnza,kjx2=2)
  parameter(kikj=ki*kj,kwork=3*(kj-1)+kj)

  !..................................................................
  !     kb (presently =8) is the maximum number of neutral injectors
  !     ke (=3) is the number of beam energy states
  !     kf is the maximum number of flux zones (or volumes) used
  !       by FREYA, and HEX** routines (reaction rate routines).
  !     kz=kf-1
  !     nap is the number of source aperatures.
  !     kprim is the maximum allowable number of primary species
  !     kimp is the maximum number of impurities.
  !     kion is the maximum number of ionic species
  !..................................................................

  parameter(kprim=ntotala)
  parameter(kimp=ntotala)
  parameter(kion=ntotala)
  parameter(kkq=kion+3,kbctim=11,kb=16,ke=3,kf=nconteqa+3, &
       kz=kf-1,kzm1=kz-1)

  !..................................................................
  !     ONETWO DIVERGENCE and source of confusion. Cray32 uses kj
  !     as the number of R mesh points, cray30 and the rest of the
  !     ONETWO code uses kj as he number of radial mesh points.
  !     We will retain the meaning as in cray32. In some of the
  !     common blocks below which are lifted from ONETWO kj
  !     has been changed to k_. These arrays are not used in CQL3d,
  !     but are retained to maintain continuity with ONETWO.
  !..................................................................

  parameter(k_=3,k_m1=k_-1)
  !     WARNING: DO NOT ALTER nap UNLESS IT IS ALTERED ALSO IN SUB ROTATE
  parameter(nap=10)
  !        parameter(kjp=maxp+1) YuP-101211: Not used?

  !     ibytes is bytes per integer word, for 64 or 32-bit integers.
  !     (1 byte contains 8 bits).
  !     It is used for packing data for urf subroutines.
  !     jjxa is 1st multiple of 8 greater than jxa.
  parameter(ibytes=8/machinea)
  !        parameter(jjxa=((jxa-1)/ibytes)*ibytes+ibytes)
  !        parameter(ipacka=jjxa/ibytes*nrayelta*nraya) !moved to urfalloc

  !     Set up new dimensioning for ifct1_,ifct2_ (from previous
  !     1 byte per word to 2 bytes per word, improving accuracy.
  !     BH, 050812).
  !     ibytes16 is number of 16-bit words per integer word.
  !     ipack16a is number of integer words required to store 1 set
  !     of ray data (*jjxa) in the ifct1_,ifct2_ 16-bit-word arrays.
  parameter(ibytes16=ibytes/2)
  !        parameter(ipack16a=jjxa/ibytes16*nrayelta*nraya)!->to urfalloc

  !BH070118      parameter (nrada=129,ninta=8,nint1a=ninta+1)
  parameter (nrada=nnra,ninta=8,nint1a=ninta+1)
  !
  !.......................................................................
  !     JK - new for addition of ADAS
  !.......................................................................
  integer kcm, kcmp1 ! kcmp1 used in getsgxn, nbsgold, nbsgxn,
  ! wrap_xboley
  integer kbe, ksge
  parameter(kcm=3,kcmp1=kcm+1,kbe=kb*ke,ksge=20)
  !.......................................................................
  !     maximum number of options in tdoutput (for nlotp1,.. arrays)
  !.......................................................................

  parameter(noutpta=10)

  ! constants from random places
  real(c_double), parameter :: ep100 = 1.d+100
  real(c_double), parameter :: zero=0.d0
  real(c_double), parameter :: one=1.d0
  !      pi=3.141592653589793d0
  real(c_double), parameter :: pi=atan2(zero,-one)


end module param_mod
