c
!     a_change.h
!
!
!***********************************************************************
!
!     This file documents changes in the code 
!
!***********************************************************************


![290] Added namelist variable bctimescal, which scales the bctime(),
![290] for background plasma time-variation.  [BH190102, YuP190309].

![289] version="cql3d_cswim_180101.3"
![289] For nlrestrt="ncdfdist" or "ncregrid", added resetting
![289] of elecfld(0:lrz) and elecfldb from distrfunc.nc.
![289] Also, disallowed reading of distrfunc.nc with multiple f(,,,)
![289] since there are improper netcdfshort settings. [BH181025].

![288] Added new option for saving distr. function 
![288] into mnemonic.nc file, netcdfshort="lngshrtf",
![288] which saves f at selected time steps only;
![288] The steps are set by specifying nplot()
![288] values in cqlinput.
![288] The value of time is recorded into 
![288] 'time_select' variable in the nc file.
![288] YuP[2018-09-28], BH.

![287] Corrected netcdfrw2.f for storage of arrays
![287] denra, ucrit, eoe0.  YuP[2018-09-24]

![286] Noticed that sometimes Te is getting to 
![286] negative values at plasma edge.
![286] It is related to interpolation by subroutine tdinterp().
![286] When tdinterp() is called by profiles.f, 
![286] at the last (upper) point in ryain-coordinate, 
![286] it uses iup="free" boundary condition, which means 
![286] the first derivative is calculated 
![286] by fitting a cubic to 4 last points
![286] in original arrays. 
![286] This procedure occasionally gives an interpolated value
![286] of Te() at the last point in rya() which is zero or negative value.
![286] To avoid such condition, an additional option is added 
![286] for calling tdinterp() with iup="linear", 
![286] which means that the derivative 
![286] is set from the last two (outer radius) points in original array,
![286]   cd2(nold)= (y(nold)-y(nold-1))/(x(nold)-x(nold-1)), 
![286] where x() corresponds to ryain() array, nold==njene.
![286] Now this subr. is called with tdinterp("zero","linear",...)
![286] for all the spline plasma profiles (in subr. profiles and tdinitl).
![286] YuP[2018-09-19]


![285] Improved the definition of reden() through zeff and 
![285] density of other species. In original definition, 
![285] reden could get a small negative value,
![285] because of a rounding error. Added lower limit =0.d0,
![285] reden(k,ll)=max(reden(k,ll),zero) .
![285] YuP[2018-09-18]

![284] Fixed bug in diagscal pertaining to scaling density of species k
![284] to specified time-dependent values, for spline-t and 
![284] lbdry(k)=scale or consscal. [YuP,BH180918]

![283] Added enerkev and uoc, corresponding to x(:) array to
![283] .nc output file (for convenience), in runaway plots.
![283] Bug fix: jfl reset to jx for sourceko='enabled' affecting fle,
![183] but with no noticed effect except proper termination of run.
![183] [BH180706]

![282] version=cql3d_cswim_180101.2
![282] For rdcmod="format1"/"aorsa", increased the number of input
![282] RF diffusion coefficient files that can be read, 
![282] and applied to the same or separate general species.  This
![282] modification also enables rdcmod with multiple general species.
![282] Also adjusted code for multiple general species restart, 
![282] using nlrestrt='ncdfdist',nlwritf='ncdfdist'.
![282] The rdcmod is not (yet) introduced for rdcmod="format2", but
![282] is readily extended to this case as need arises.
![282] Namelist nrdc.ge.1 is the number of rdcmod files.  See
![282] cqlinput_help for additional new namelist variables,
![282] rdcfile(1:nrdc),nrdcspecies(1:nrdc),rdcscale(1:nrdc),
![282] rdc_plot,rdc_netcdf.   [BH180608].

![281] version=cql3d_cswim_180101.1
![281] Added capability to make plots in color, particularly 
![281] contour plots for distr.function, change of distr.func.
![281] over time step, urfb coeffs, and source.
![281] Main changes are in pltcont.f.
![281] The color can be added by setting
![281] pltd='color', or pltd='df_color', or pltso='color',
![281] pltso='first_cl', or plturfb='color'. 
![281] BH,YuP[2018-02-08].

![280] YuP[2018-01-08] Revised some of plt*.f files, also souplt.f,
![280] to correct the functionality of pltlim plotting options.

![279] version=cql3d_cswim_180101.0
![279] YuP[2018-01-05] Added resetting of vth() array in profiles.f.
![279] vth is the thermal velocity =sqrt(T/m) (at t=0 defined in ainpla).
![279] But, T==temp(k,lr) was changed above, 
![279] in case of iprote (or iproti) equal to "prbola-t" or "spline-t".
![279] Then, reset the values of vth() [used for plots, and also 
![279] in lossegy, efield, restvty, sourceko, tdnpa, tdtrdfus,
![279] tdoutput, and vlf*,vlh* subroutines].
![279] Similarly - for the energy() array.
![279] The energy() array, initially set as 1.5*temp() [see ainpla.f]
![279] is re-defined at each time step in subr. diaggnde
![279] as the m*v^2/2 average over distr.function in vel. space.
![279] BUT, it is done for k=1:ngen species only!
![279] So, in profiles.f, we re-define the energy() array  
![279] for the Maxwellian species only.

![278] YuP[2018-01-04] Adjusted the printout of variables 
![278] that are specified in cqlinput into *.ps file.
![278] Now the long text lines are wrapped around 
![278] to additional lines, instead of just cutting them.
![278] See ainplt.f, lines ~~180-240.

![278] YuP[2018-01-04] Adjusted subr.pltstrml, 
![278] to make plots of streamlines of phase flow.
![278] Also - small corrections in pltcont.f, for the plots 
![278] of v=vnorm and v=vth lines in different units (v/vnorm or v/c).

![278] YuP[2018-01-02] Corrected plots of pitch angle avg source vs u.
![278] In subr.pltsofvv (file souplt.f), added (similar to pltcont):
![278] Setup different types of horizontal axis, depending on 
![278] namelist settings of pltlim = 'u/c' or 'energy' or 'x'.

![277] BH,YuP[2018-01-02] Added resetting of ipro*** values
![277] to "prbola-t", in case when the run
![277] is done with nbctime>0 (time-dependent profiles),
![277] but ipro*** values are set to "parabola".
![277] See ainsetva.f.
  
![276] version=171124.1
![276] Adding time-dependent scale factor, difus_io_drrscale(,k), 
![276] difus_io_drscale(,k) at times  difus_io_t(1:ndifus_io_t). 
![276] This is applied only in difus_io(k)="drrin" and "drrdrin" cases.
![276] Fixed several out-of-bounds bugs affecting the lrz=1 cases.
![276] [BH171231].

![275] version=171124.0
![275] Added capablility to output and input the radial diffusion 
![275] coeffs d_rr and pinch arrays as function of theta,u,species,
![275] rho to/from a netCDF file.  This enables outputing a
![275] template file, and inputing numerical d_rr and d_r.
![275] Additional namelists in the trsetup namelist are difus_io(1:ngen) 
![275] and difus_io_file. Added for MST simulations. [BH171122].
![275] Installed neutral t-dependent splines (ennin_t) from FOW cql3d.
![275] [YuP171127].
![275] Also, added some modifications added to facilitate TRANSP
![275] compile system, labelled CMG (Marina Gorelenkova), from YuP work.
![275] Added improvement to eqfndpsi.f according to cql3d-fow, which
![275] fixed a problem in setting up a psi(rho) grid for some eqdsk
![275] cases with small psi gradient near the magnetic axis.[BH171211]
![275] Changed the default of nconteq="psigrid" to nconteq="disabled",
![275] and set nconteqn=50.  See cqlinput_help.   [BH171211]

![274] YuP[2017-11-17] Migrated eqrhopsi from the CQL3D-FOW version. 
![274] It has corrections for eqpsi(1:nconteqn) definition, see do_11 loop.
![274] In cql3d-fow (or in "mirror" version this adjustment 
![274] is done on [03-2016; with correction on 2017-12-05]. 
![274] This correction gives an effect on 
![274] NBI deposition points (very little, but visible - 
![274] in coords of points in R,Z plane).

![274] YuP[2017-11-21]: In definition of sorpw_nbi() [Neutral Beam power]
![274] added if(gone()) check: ONLY COUNT not-lost particles.
![274] It is done similar to the CQL3D-FOW version (added on 03/23/2015).
![274] The array sorpw_nbi is only used 
![274] for plots and NetCDF file. 
![274] It has no effect on the source operator or solution.
![274] The consequence of adding if(gone..) condition is
![274] the drop of NBI power profile at plasma edge 
![274] (if lossmode is not disabled) because of banana+Larm.rad losses;
![274] see tdpltjop, plot 'FSA SOURCE POWER DEN: ...';
![274] Also, 'Power from NBI', and 'Power integrated over radius'
![274] become smaller (values are printed in *.ps file). 

![273] version="cql3d_cswim_170101.6"
![273] Adjusted cql3d to accomodate inclusion in TRANSP.
![273] Routine to detect if &fsetup namelist being used instead of
![273] the standard %setup0, and provide for reading it.
![273] (&fsetup was an accomodation required by the pathscale
![273] compiler, since only was sensitive to 5 characters in a
![273] namelist section name.)
![273] Explicit deallocation of arrays before end of cql3d,
![273] rather than just depending on the operating system to deallocate
![273] at the end of a run, adding subroutine de_alloc in it3daloc.f.
![273] Additional zero initializations.   [Petrov, Nov/2017].

![272] Migrated some of recent modifications from FOW version:
![272] 1. Updated freya subroutines in cql3d to calculate NBI
![272] stopping according to ADAS data in GA freya code[Kinsey, 130821].
![272] 2. Capability to read a list of neutral beam injected
![272] ion starting points from NUBEAM.  List is generated by Deyong
![272] from an NSTX TRANSP run for a modulated beam case.  
![272] New frsetup namelist: read_birth_pts, nbirth_pts_files,
![272] birth_pts_files, nbirth_pts [BH130304].
![272] 3. NBI pulse square-wave modulation.  New namelist variables:
![272] beamplse,beampon,beampoff.  Output of time-averaged 4D (f4d)
![272] distribution function.  New namelists: tavg, tavg1, and setting
![272] f4d_out=tavg.  First use is NSTX NBI+HHFW [BH120608].  
![272] 4. Correction for B-field signs by introducing cursign and bsign
![272] factors  in freyasou.f, tdnpa0.f, and tdsxr0.f.  It may not
![272] have a large effect for cases where the toroidal field is dominant,
![272] and not effect if cursign/bsign are positive. [YuP and BH, 110815].
![272] 5. Ampere-Faraday solution option (see 'ampfmod').


![271] version="cql3d_cswim_170101.3.1". Young-soon Bae, KSTAR, found
![271] problem with updated ray files generated by cql3d, when using
![271] multiple ray input files.  Second file overwritten with first 
![271] file data, for a urfmod, netcdfshort='longer_f' case.  
![271] Fixed by including a krf dependency in output of the ray data 
![271] in netcdfr3d.f.  [BH, 170920].

![270.5] Found that the unphysical looking ledge in the distribution
![270.5] in wh80.ps knockon/runaway electron case was due to a bug fix in 
![270.5] sourceko.f on 1/23/2000. (Shown by reverting the fix.) Need to check
![270.5] out further the reason for the "ledge".  [BH170927].

![270] version="cql3d_cswim_170101.3".  Restoring original 
![270] lossmode(k)="simplban" functionality.  lossmode(1)= 'simplbn1'
![270] needs work.  lossmode(k)='simplbn2' gives same results as 'simplban',
![270] but plan to update it to use subroutine deltar first order calc
![270] of radial orbit shifts.  deltar only works with ngen=1 at this time.
![270] Fixed bug in tdrmshst.f which was giving incorrect area/vol
![270] mesh increments, resulting in incorrect total currents and powers.
![270] Also, removed fow_orbits.f from ZOW code.  [BH and YuP, 170717].

![269] version="cql3d_cswim_170101.3"
![269] Re-setting namelist lbdry(k)="scale" to lbdry(k)="consscal" in
![269] in ainsetva.f. The "consscal" model for density rescaling
![269] is more straight-forward to interpret.
![269] This adjustment fixes a problem in lbdry(k)="scale" which arose
![269] due a fix of a "nipple" forming on the distn at v=0.  The "nipple"
![269] let to constant density rise in a lbdry(k)="scale" modeling
![269] of LHCD in C-Mod.  Problem identified by Syun\'ichi Shiraiwa.
![269] [BH and YuP, 170703].

![268] version="cql3d_cswim_170101.2"
![268] YuP[04-2017] Adjusted files in /mpi/ subfolder
![268] to match the recent changes in urf* files.

![267] Following the message from John Wright[04-2017, 04/03/2017]:
![267] Changed index 'n' in statement functions in subr.cfpcoefr
![267] to 'ng'.  Index 'n' is also in comm.h, which causes 
![267] a compilation error (bug) in new Intel compiler (2017).

![266] version="cql3d_cswim_170101.1"
![266] Regression tests of a ngen=2 (2 ion general species) test with
![266] and NBI heating and multiharmonic FW heating of both general
![266] species, showed that results had significantly deviated in an
![266] unphysical way from prior 130124 results.  The changes were tracked
![266] down to problematic changes in losscone.f and a change in the
![266] treatment of a resonance contribution (rrr= in the code) in the
![266] QL diffusion (urfb0.f) and damping (urfdamp1.f).  Coding was
![266] reverted to the long-standing approaches, pending further 
![266] studies of the changes.  [YuP and BH, 170105].

![265] version="cql3d_cswim_161220.1"
![265] YuP provides cql3d_cswim_160605_update1.zip, updating
![265] cql3d_cswim_160602 version  from cql3d_160720_mirror 
![265] which incorporates changes for the cql3d-FOW version in
![265] urfb0.f etc.  This is a fix for a problem found by Young-Soon Bae
![265] with damping on 3rd harmonic EC, for top launch of mainly
![265] 2nd harmonic EC in Kstar.  Updated sigsetup.f according to 
![265] YuP changes in Aug/2016 [YuP, with BH, 161220]
![265] More details on the updates:
![265] -- Removed storage arrays for <E> and <F> QL coeffs, 
![265] as they are expressed through <B> and <C> coeffs.
![265] (For urf modules only).
![265] -- Upgraded vlf and vlh modules: now they are valid for multi-surface
![265] runs (arrays cqlb and other cql* now include additional radial index).
![265] -- Value of symm is defined in aingeom.f now, and then used in other
![265]  files (originally - defined in many places/files). 
![265] -- Value of truncd is defined in micxinit.f now, and then used in 
![265] other files (originally - defined in many places/files). 
![265] -- Plots of QL coeffs in urfbplt/pltcont: the all-modes plot 
![265] is modified to sum-up all modes and then plot contours of the 
![265] total array (originally - overlapping contours of separate modes in one plot).
![265] -- In mixcinit: Check that there are sufficient number of v-grid points  
![265] over thermal part of the distribution (at least 3 points).
![265] If not, print warning [see if(j_thermal.lt.3) ...].



![263] version="cql3d_cswim_160602"
![263] Upgraded lossmode(1)='simplban' so banana (+gyro) losses are based
![263] on the first order deltarho orbit shift.  Kept the old pre-160403
![263] model (='simplbn2'), and an intermediate model (='simplbn1') for
![263] backward compatibility and cross-checking.  Main changes are in 
![263] losscome.f.  [See file log.1_simplban_extended.gz for comparisons.]
![263] Should use lossmode(1)='simplban' in combination with ndeltarho=
![263] ='freya', to prevent unphysical scrapeoff of beams. [BH, 160602].

![262] Fixed inadaquecy in tdnpa0.f which was causing NaN for
![262] some cases for efluxt.  Similar fix for tdsxr0.f [BH160508].

![212] The total integrated powers - same as with original version.
![212] The older expression for pwrrf is still present, commented out.


![211] version="cql3d_cswim_160312.0"
![211] Bug fixed in tdrmsst.f which was causing EC absorption to
![211] to be (impossibly) radially inboard of the closest approach
![211] of the input EC rays to the magnetic axis, for a eqsym='none'
![211] case [YuP and BH 160312].

![210] Attemping to get cqlinput_mirror_ions_IC.10 case going:
![210] Uses Gary Smith b-field model, single flux surface, added couple 
![210] of logic lines (See cBH151202) but then realized more variables
![210] need setting for ICRF vlfmod to work in this environment.
![210] Eventually, single flux surface, particle source + RF QL 
![210] worked cqlinput_mirror_ions_IC.4.8.4) (diffusion(BH151202).

![209] Enabled various capitalizations of species namep/namei, for 
![209] operation with Plasma State software [YuP, 150620].

![208] YuP[2015/05/03] version=cql3d_fullFOW_140425.1
![208] Several adjustments in subr.eqorbit().
![208] 1. Set a limitation/resetting for the value of epsicon_
![208]   which is the input argument for eqorbit().
![208]   If it is too small (outside of LCFS), reset with warning:
![208]     epsicon_= psilim + 0.001*(psimag-psilim)
![208]   The shifting fraction, 0.001, could be changed.
![208]   The starting point for tracing the surface
![208]   should not be too close to the LCFS (psi=psilim),
![208]   otherwise the "orbit" of the surface may diverge 
![208]   outside of the LCFS.
![208] 2. Modified the adjustment of 1st/last point:
![208]   At small rho, the "orbit" usually "drifts" towards magn.axis,
![208]   so that solr_(nn)<solr_(1).
![208]   But at outer surfaces, the "orbit" diverges to the outside,
![208]   so that solr_(nn)>solr_(1).
![208]   In such a case, the orbit can get outside of LCFS,
![208]   if the starting point solr_(1) is close to the LCFS.
![208]   So, setting solr_(1)=max(solr_(1),solr_(nn)) would have
![208]   a bad result. Simply move the last point to the starting point:
![208]   solr_(nn)=solr_(1)
![208] 3. Sometimes, at ending point, in case of eqsym='none',
![208]   an orbit/surface cannot get to the starting point
![208]   because of accuracy of the orbit integrator.
![208]   Orbit can get to a smaller R, where the value of bval/bmidplne_
![208]   is a bit smaller than at the starting point, say 0.9999
![208]   instead of 1.d0. Although the ending point is forced 
![208]   to be the same as the starting point (later in the subr),
![208]   the previous point (lorbit_-1) may still have the value of
![208]   bpsi_(lorbit_-1) = 0.9999 or so.
![208]   In any case, bpsi_(l) cannot be less than 1.d0 - 
![208]   the code cannot handle such equilibria.
![208]   The adjustment that is made only works for eqsym=none
![208]   (during 2nd half of surface):
![208]   Adjust (R,Z) points by setting them in a straight line 
![208]   connecting the first point where bpsi_<1 is detected 
![208]   (designated l=lbad) to the end/starting point [to be exact, 
![208]   connecting the point(lbad-1) to the point(lorbit_)]
![208] There are few more adjustments for the last two points 
![208] in surface tracing (mostly for the case of eqsym='none').
![208] Search for 2015/05/03 in eqorbit.f.

![207] YuP[2015/05/03] Added plotting of the LCFS in plots
![207] with flux surfaces, in case of eqsym='none'.

![206] YuP[2015/05/03] Bug fix in eqrhopsi: 
![206]   ez2=ez(jmag+1) ! was (imag+1)
![206] It was affecting the tracing of small flux surfaces
![206] near magn.axis in some cases of nonsymmetrical 
![206] up-dn equilibria (eqsym='none'), because of error in 
![206] determination of (supposedly more accurate) value of zmag.

![205] YuP[2015/05/02] Added a more accurate def. of volume and area
![205] near magnetic axis.  This definition 
![205] is not dependent/sensitive to tracing the surface.
![205] Gives better convergence in eqfndpsi in cases for 
![205] radcoord.eq."sqtorflx", when 'rhonew' is proportional to 
![205] sqrt of 'volum' (near magn.axis).

![204] YuP got rid of the nipple at region of f(v=0) which was
![204] appearing in restart runs due to that [YuP believes] f(j=1) is 
![204] not found as the solution, but taken from the previous time step.
![204] YuP[04/10/2015].

![203] Slight correction of simplban loss model: use local bmod0
![203] rather than edge bmod0, for gyro-losses.
![203] Enabled use of density rescaling along with lbdry(k)="conserv",
![203] using new option lbdry(k)="consscal". Previously, density
![203] rescaling only occured with lbdry(k)="fixed". [YuP+BH 140806].

![202] version="cql3d_cswim_140804"
![202] Fixed bug in calc of cfpgamma (Coul log) for diagnostic
![202] purposes in tdoutput.f (tdnflxs not called in l=1,lrz loop).
![202] (Also, lr_ improperly used in a do loop variable.) BH140804.

![201] version="cql3d_cswim_131030"
![201] Reduced size of  lossfile, netcdfnm, eqdsk, rffile() from 
![201] charcter*512 to character*256, to avoid possible limitation
![201] with a character read of variable greater than ~300 with intel
![201] compiler.
![201] Fixed bug in input time-dependent parabolic profiles (tdxin13d.f).
![201] Fixed bug in netcdf file writes for wpar,wper, and longer_f
![201] save of distributions (netcdfwr2.f)  [BH131026-30].

![200] version="cql3d_cswim_130928_mpi"
![200] Added delim='apostrophe' to open (iunwrif,... in tdwritef.f.
![200] This is to address a concern/problem about long character variable
![200] when writing/reading them as namelist elements with various
![200] compilers.  [gfortran has this delim as default, but this 
![200] convention is not universal].  Checked namelist variables
![200] lossfile, netcdfnm, eqdsk, rffile(), and all are set to 
![200] character*512  (enabling long path specs). Also, checked that
![200]  all open() for writing namelist contain the delim statement. 
![200] Could reduce variable length per Bonoli, but trying this first.
![200] [BH130928].

![199] Restored printout of rfpwr for rdcmod.ne.disabled cases
![199] (fixing oversight in ~2011 when added printout of more than
![199] 3 harmonics of rfpwr for urfmod.eq.'enabled') [BH130611].

![198] version="cql3d_cswim_121120_mpi"
![198] Fix error exit hang for MPI, discovered by John Wright.
![198] [YuP, 121117].

![197] Modified tdsxr0.f and tdnpa0.f to handle a sxr/npa sightline
![197] which occurs very near to the magnetic axis, when the
![197] magnetic axis psi value is slightly underestimated [BH121115].

![196] version="cql3d_cswim_120816_mpi"
![196] Added powrft (rf power summed over modes) to use by SWIM IPS.
![196] Some cleanup of netcdfrw2.f. Additional storage in tdsxr0.f.
![196] Reverted time step counter update to previous use [YuP+BH, 
![196] an MPI related mod, now not necessary]  [BH120816].

![195] Added option to specify radial diffusion profile as tein, nein, 
![195] etc, using a new namelist variable: difin
![195] If sum(abs(difin))>0, then difin(ryain) is used to define rshape 
![195] in tdtrdfus.f instead of difus_rshape(1-7)
![195] This is useful since it has been found that a practical profile
![195] is D0 * chi_e(r). So the chi_e profile from power balance
![195] can be given in difin and D0 in difusr. TCV experiments showed 
![195] that D0~0.2 works for many cases [OS 20120809]

![194] version="cql3d_cswim_120807_mpi"
![194] Fixed bug preventing power from more than two neutral beams to be
![194] to be injected [BH120807].

![193] version="cql3d_cswim_120315_mpi"
![193] Reworked and debugged nuclear fusion rates.  The DD-rates remain
![193] unchanged, but D-T, D-He3 were not working.  <sigma-v> results
![193] now agree within a few percent with standard Bosch-Hale results,
![193] for the four nuclear fusion rates in the code.  Logic was adjusted
![193] in subroutines sigsetup, sigmax, sigv5d.  Namelist variable 
![193] isigsgv2 was made inoperative, as has no apparent physics use
![193] [BH120315].

![192] Add netcdf ouput of FSA current parallel curr.  Fixed minor bug 
![192] re printout of Fisch normalized RFCD [changing value].  
![192] Fixed bug in output of rfpwr to netcdf file. Power density
![192] in individual modes, at total (not sorpw/sorpwi) were
![192] erroneously divided by dvol() [BH120221-3].

![191] version="cql3d_cswim_120202_mpi"
![191] Added required blas and lapack routines to r8subs.f, so
![191] no further need for loading blas and lapack libraries.
![191] Updated recent makefiles.  Added some makefiles for PPPL.
![191] Add powurf and powurfc to netcdf .nc file for all mrfn urf modes.
![191] Corrected defn of rfpwr in .nc file, and added data to all "modes"
![191] [BH120202].

![190] version="cql3d_cswim_120124_mpi"
![190] For multiple general species, added capability to damp same urf
![190] wave data sets on multiple species.  Previously, only could do
![190] this for one general species (ngen=1)  [BH120125].

![189] version="cql3d_cswim_120122_mpi"
![192] noplots="enabled1" inhibits all plots.  Previously, 
![192] noplots="enabled" inhibited all GRAFLIB plots, but over a
![192] period of ~ 10 years, all such routines had been converted to
![192] PGPLOT library routines, so it has no effect.  The purpose of
![192] noplots="enabled1" is to enable compilation without the pgplot 
![192] library [along with using, e.g.,  the -Wl,-noinhibit-exec gnu 
![192] loader (ld) option in gfortran which enables loading of the code 
![192] in the presence of unsatisfied externals] [BH121122].
![192] Also increased max number of NBs, kb=3 ==> 8 (for D3D+).

![191] Increased iotalim from 750 to 6400, for use with larger
![191] values of jx [BH120104].

![190] Fixed an apparent bug in symmetrization in the calculation of
![190] collisional contribution to C,D, and F in the trapped region. 
![190] [See cfpcoefn.] The effect on distributions for a test case
![190] seemed small,  but more checking is warranted [YuP111202].
![190] Had no effect in an aorsa-dc-cql3d simulation [BH].

![189] Fixed write of two distributions to the .nc file, for
![189] ngen=2 cases [BH111103].

![189] version="cql3d_cswim_110804_mpi"
![189] Small mods, bug fixes, improvements:
![189] -option added for output of freya NBI birth points to ascii file
![189] -clipping option to smooth external diffusion coeffs from DC
![189]  or AORSA added in rdc_multi.
![189] -removed some parameter-setting restrictions in netcdfrw2.f,
![189]  persuant to dynamic dimensioning modification of cql3d.
![189] -Removed some redundant "save" statements in several subroutines
![189]  and added a character*8 ntitle,dat in equilib.f, following
![189]  Jacob Urban compilation with g95.
![189] -Added deltarho, first order orbit shift to the netcdf output.
![189] -diagimpd.f changes to trunc-d option to improved treatment of
![189]  NBI source into prompt loss region.
![189] -maxp maximum number of freya NBI source particles increased
![189]  from 150000 to 1500000   [BH and YuP, 110804].

![188] Updated NPA recombination rate in tdnpa.f from Path and Wolf
![188] formula to Hutchinson, Eq. 6.3.5, per Aaron Bader. Not much
![188] change in C-Mod NPA spectra, since recomb. is small [BH110804].

![187] Inconsistency in distribution function extrapolation when
![187] increasing enorm, was fixed.  Had caused neg density of
![187] regridded density, and code to halt, is some cases [BH110525].

![186] rdc_clipping functionality added.  Imported rf diffusion coeffs
![186] such as from the DC code (via rdcmod=enabled), can have have
![186] unphysical spikes near the trapped passing boundary.  These can
![186] be removed using a running median filter. See code for refs
![186] [BH110430].

![185] Slight addition of ineg='enabled1' option, removing spikes
![185] in the distribution occuring when source points are in
![185] prompt loss regions of velocity space [BH 110408].

![184] Add printout of freya NBI birth point positions and velocities
![184] [BH110406].

![183] Code version='cql3d_cswim_110401_mpi'
![183] Added mpi capabilities for X-ray, neutron, and npa diagnostics
![183] [YuP110404].

![182] Enabled that number of Legendre polynomials in the approximation
![182] of f and sigma for xray calc, msxr, can be set independently
![182] of number of Legenedre polynomials in the FP collision
![182] coefficients [previously msxr.le.mx]. mx=3 is probably OK
![182] for collsion coeffs, but msxr~8 is required for accurate calc
![182] of XR spectra due to electron runaway distributions. 
![182] Similarly, the number of Legendre polynomials used in fusion
![182] rate calcs, mmsv, in no longer limited to .le.mx.  Changing
![182] from mmsv=3 to 8, though, made little difference in neutron
![182] rate, for a NSTX NBI+HHFW case.  [BH110401].

![182] Enabled MPI version of code to work for transp="enabled" and
![182] the default soln_method='direct'.  (Additional work is required
![182] for parallelization of the full, soln_method='it3drv' case
![182] [YuP and BH 110329].

![181] Added coding to include multi-species data in the output
![181] netCDF file, for cql3d runs with two or more general species.
![181] Output now includes multiple distribution functions, and
![181] time-dependent wpar,wperp, specific power densities and
![181] and currents [BH110321].

![180] Adjusted coding for larger number of sxr/NPA view cord
![180] plot points.  A problem appeared due to dynamic dimensioning
![180] of some temporary storage, previous set to a larger
![180] dimension by fixed parameter values [BH110318].

![179] Code version="cql3d_cswim_110315_mpi". Fixed two
![179] bugs in the freya NBI modules: (1) For single source
![179] NBI cases, using a source axis centered aperature, s-xxxx,
![179] resulted in undefined aperature widowing (iap not defined
![179] in zfreya.f/rotate). Modified ashape namelist defn slightly
![179] so in nsourc=1 case, can enter centered aperature dimension
![179] with ashape(,)=s-xxxx (source centered) or s-xxxx, equivalently;
![179] (2) namelist smooth not passed properly from frcomm.h to comm.h.
![179] Changed name to smooth1 in comm.h and elsewhere.
![179] Affected smoothing of NB source profile.  [BH and YuP, 110314].

![178] Code version="cql3d_cswim_110308_mpi". Should have same 
![178] functionality as serial version cql3d_cswim_110308.
![178] The line  if (noplots.eq."enabled") return 
![178] in pltprppr.f is commented out, so that
![178] the parallel distribution function will be plotted 
![178] unless pltprpp="disabled" [YuP 03-07-2011].

![177] Achief1.f (used for a single flux surface calculations)
![177] is not used anymore. Now the solution is found through 
![177] the call of achiefn(0) in ll-loop in tdchief, same as for
![177] multi-flux-surfaces calculations. 
![177] The probable reason for using achief1 was the problem with 
![177] tdxinitl, called from tdinitl. Now tdxinitl, which "fills in" 
![177] input data parabolically and computes the normalized radial mesh,
![177] is only called for lrzmax>1 [YuP 03-07-2011].

![176] Modification made after [03-2011] to ineg="enabled" or "trunc_d": 
![176] (See diagimpd.f)
![176] If distr.function f(i,j) is negative or zero 
![176] at some (i,j)-point in vel. space,
![176] it is set to zero for ALL jj.ge.j 
![176] (at fixed i-index of pitch angle),
![176] for energies greater than the maximum in the source term
![176] (for example from NBI).  For lower energies, neg f is set = 0.
![176] Before this modification, it was set to zero 
![176] for only this specific (i,j)-point where f(i,j)<0.
![176] So, no more "islands" in distribution function 
![176] remaining beyond such (i,j)-points, except when they are 
![176] caused by sources like NBI.


!================================================
![175] MPI-enabled version.   February 2011 [YuP].
![175] The parallelization capabilities of the code are upgraded. 
![175] The MPI executable is produced by launching
![175] make -f makefile_mpi.franklin 
![175] (similarly for Hopper NERSC machines)
![175] The makefile will run Python script doparallel.py 
![175] (in /mpi/ subdirectory) that converts source files 
![175] by inserting MPI commands from mpins_par.f.
![175] The location of such commands in original source files  
![175] can be tracked by searching all phrazes that start with CMPI.
![175] The executable can be launched as a batch job (example):
![175] qsub -q regular franklin_batchscript_mpi128
![175] The file franklin_batchscript_mpi128 should contain
![175] #PBS -N xcql3d_mpi.franklin
![175] #PBS -q regular
![175] #PBS -A m876
![175] #PBS -m e
![175] #PBS -j oe
![175] #PBS -l walltime=00:05:00
![175] #PBS -l mppwidth=128
![175] cd $PBS_O_WORKDIR
![175] aprun -n 128  /???/xcql3d_mpi.franklin
![175] (Specify your working directory instead of /???/)
![175] 
![175] The parallelization is done for:
![175] 1. Impavnc0 solver, together with collisional coefficients 
![175] generator; see achiefn.f. Parallel run for each flux surface.
![175] 2. Energy transfer diagnostics; see diagimpd.f, 
![175] call diagentr(lefct,k). Parallel run for lefct=-1:3,5:8,11,12.
![175] 3. Calculation of diff. coefficients for the ray-tracing;
![175] see urfb0.f. Parallel run for the number (mrfn) 
![175] of excitation wave modes. 
![175] 4. Calculation of damping of rays; see urfdamp1.f and urfdamp2.f.
![175] Parallel run for combined (number of rays)x(number of wave modes)
![175] which is nrayn*mrfn (these two numbers are determined from 
![175] reading the rays data file).
![175] From above, the optimal number of cores (ranks, or mppwidth) is
![175] the largest of: lrz+1, 11+1, or nrayn*mrfn+1 (if rays are used);
![175] [+1 because rank=0 does not perform calculations, only collects
![175] data from other ranks].  More cores can be requested, but
![175] extra cores will be idling.
![175]
![175] The parallelization of the impavnc0 solver is only done 
![175] for soln_method="direct".  In other cases of soln_method, there 
![175] will be not much speed-up because the solver will be running 
![175] on rank=0, for all flux surfaces; some speed-up can still occur
![175] if many rays are used, due to parallelization of urf-routines.
![175] For soln_method="direct", the speed-up can be ~10 times.

![175a] Changes made in course of MPI-upgrade:
![175a] Call_diagscal is moved out of diagimpd.f. 
![175a] Updating of dtr, dtreff, dttr, n_(l_) is moved out of achiefn.f
![175a] to tdchief.f. 
![175a] Call_profiles is moved out of achiefn.f to tdchiefn.f, 
![175a] before ll-loop.
![175a] Definition of gfu,gft and gfi functions is moved out of advnce.h
![175a] to the end of diagentr.f as function subprograms, 
![175a] to avoid cvmgt-construct.


![174] Removed reading of restart namelist data from distrfunc file
![174] for restarts with nlrestrt='ncdfdist' or 'ncregrid'. BH110201.

![173] Added restart option nlrestrt='ncregrid' enabling change in
![173] cqlinput of enorm, jx and xfac from values in the restart
![173] .nc file.  The restart distribution is extrapolated and
![173] interpolated onto the code momentum grid. BH110201.

![172] [January 2011]
![172] Adjusted plotting limits for flux surfaces in tdplteq.f.
![172] In some cases of up/dn symmetrical surfaces, when only half of a
![172] surface is plotted, the surface can be in negative-Z area. YuP

![171] Dynamically allocated local arrays in tdinlegw.f - 
![171] allows using large iy now.  YuP

![170] Changed default value nonboot=1 to nonboot=2  in aindflt.f 
![170] (turn on computational bootstrap at n=nonboot).
![170] The radial derivative of distr. function 
![170] for bootstrap current calculations uses f() at neighbouring
![170] flux surfaces; in general, it is known from previous time step;
![170] during n=1 time step the distr.function is still zero on many
![170] flux surfaces until completion of the time step.
![170] So it\'s better to start bootcalc at 2nd time step. YuP

![169] Problems in tests with bootcalc problem: need high resolution
![169] in theta (pitch-angle). This can only be achieved by setting 
![169] nii=25 in baviorbt (baviorbt is used when taunew='enabled').
![169] nii is the number of pitch angle subintervals used in the
![169] calculation of dtau and tau. 
![169] Using nii=1 was ok for most tests, except bootstrap calculations.
![169] Since higher nii does not add much of computational time, 
![169] nii=25 will be the default value.
![169] With nii set to 25 in baviorbt (and iy=160 theta grid), 
![169] the bootstrap current profile is smooth.  
![169] Also restored factor of 2 for the first node of poloidal-grid:
![169] if(l.eq.1) dtau(i,l,lr_)=dtau(i,l,lr_)*2.0 
![169] !-YuP Should be *2 because first node is over half-interval.
![169] Result: Better shape of tau at pitch-angle theta~pi/2. YuP

![168] Moved z00 function from advnce.h  back to the end of impavnc0.f.
![168] Added k-index into z00, because 
![168] di(i,j,k,l_) and dj(i,j,k,l_) have k-index.
![168] This modification allows to reduce the usage of cvmgt function -
![168] better optimization. See note [162] below.  Results are the same.
![168] YuP

![167] Fixed a bug in bsl(jj,kk,ll) and bsu(jj,kk,ll) functions:
![167] x(jj) in bsl and bsu is only defined for jj=1:jx,
![167] while bsl and bsu are called with jj=0:jx+1.
![167] Made this fix:
![167]   jjj=max(jj,1)   ! to limit from below, for x(jjj)
![167]   jjj=min(jjj,jx) ! to limit from above, for x(jjj)
![167] Alternatively, could set bsl and bsu to zero for jj=0 and jx+1.
![167] Almost no effect on results, but prevents out-of-bounds error. YuP

![166] Fixed a bug in reading of Complex*16 array into a Real*8  
![166] dummy array in urfread_i.  YuP

![165] Possible bug in reading ray##/text data file.
![165] The formats in GENRAY (write3d.f) for saving data:
![165] 1    format (2i5,1pe16.9)
![165] 2    format (5(1pe16.9))
![165] 3    format (9i5)
![165] 4    format (9(' ',i4))
![165] But in CQL3D (urfread_.f) for reading data:
![165] 1    format(2i5,1pe16.9)
![165] 2    format(5(1pe16.9))
![165] 3    format(9i6)
![165] 4    format(9i5)
![165] Should we make the last two formats in CQL3D as in GENRAY?
![165] No changes for now, just keep in mind.

![164] Update 110107:
![164] Corrected error in advnce.h in bsu(j,k,l), function fpj0.
![164] Should be bsu(j,k,l_).  YuP, BH

![163] Fixed a bug with denpar and temppar - they are dimensioned as 
![163] (ntotala,0:lza+1), but in clear.f they are zeroed over lsa+1.
![163] These two arrays have dual usage: 
![163] either have radial coord. dependence, 
![163] or along-field-line dependence (when cqlpmod="enabled").
![163] The problem is fixed by forcing  
![163] parameter(lza=lsa)  in param.h.  YuP.

![162] Modifications to address the problem of optimization on Franklin
![162] and gfortran compilers. The compilers could not perform 
![162] full optimization because of functions cvmgt or similar
![162] (not an intrinsic, but a declared function in r8subs).
![162] During invoking of such function, both 1st and 2nd arguments 
![162] are evaluated, although only one is needed.
![162] This function is heavily used in definition of
![162] other functions, which are called in nested i,j,k,l loops.
![162] Tried to use intrinsic function merge() - no improvement.
![162] Reduced the usage of cvmg# functions to a minimum
![162] by using if-else statements instead.
![162] For some cases the code now runs ~7 times faster on Franklin. YuP.

![161] Update YuP-101230:
![161] Dynamic dimensioning of sounor(ngen,nsoa,lz,lrz) 
![161] to reduce memory usage.
![161] Changes lza->lz in many dynamically allocated arrays.

![160] YuP-101228: Corrected error in diagimpd.f (do 410 loop)
![160] related to ineg='enabled'. Now the negative parts of 
![160] distr.function are really set to zero. 
![160] The effect from ineg bug fix is very small  -
![160] in 3rd-4th digit. 

![159] Added dyi()=1./dy() and dxi()=1./dx(), 
![159] made changes in advnce.h and diagentr.f to re-arrange terms, 
![159] to make the code run faster.

![158] YuP-101224: Corrected error in netcdfrf.f:
![158] In pack21(y,1,iy,1,lrors,urftmp,iymax,lrors),
![158] replaced urftmp by tem1.  In some cases the size of urftmp 
![158] is smaller than size of y.

![156] YuP-101221: Eliminated parameter noncha.
![156] Arrays that depended on noncha are allocated now 
![156] using nonch which is set to nonch=nstop+1 in ainsetpa.

![155] YuP-101220: 
![155] Reduced size of many arrays dependent on ngena; 
![155] now they depend on ngen.
![155] (Also, possible bug in tdxinitl: changed tauegyz to tauegy). 
![155] Similarly - for nmodsa. Large arrays that depended on nmodsa are
![155] cqlb...cqlf and wcqlb...wcqlf. 
![155] Now nmodsa is replaced by mrfn in these arrays.
![155] Moved allocation of cqlb...cqlf to vlh.f and vlf.f
![155] where mrfn is determined. 
![155] Allocation of wcqlb...wcqlf is moved to vlf.f

![154] YuP-2010 December 08-17
![154] Other parameters eliminated: maxp,jxa,iya,mxa,jfla,
![154] and related to them.  Many changes through the code -  
![154] new version cql3d_101217.   
![154] Code runs ~2 times faster; smaller memory footprint.
![154] urfdamp2.f is re-organized to make it run faster.

![153] YuP-101208: parameter nharma is not needed anymore.

![152] YuP-101207: added in tdchief.f, just before call achiefn(1)  :
![152]   do k=1,ngen  ! Compute density gains and losses, and powers.
![152]      call coefstup(k) ! To define da...df coeffs, gon(i,j), etc
![152]      call coefmidv(da,temp1,1,vptb(1,lr_))
![152]      call coefmidv(db,temp1,2,vptb(1,lr_))
![152]      call coefmidv(dc,temp1,3,vptb(1,lr_))
![152]      call coefmidt(dd,temp1,1)
![152]      call coefmidt(de,temp1,2)
![152]      call coefmidt(df,temp1,3)
![152]      call coefwtj(k)
![152]      call coefwti(k)
![152]      call diagimpd(k) ! to calculate sgain(1:8,k)
![152]   enddo ! k
![152] Needed to compute sgain(1:8,k)

![151] YuP-101207:  Modified definition of ipack and ipack16:
![151] ipack16= jjx/ibytes16*nrayelts*nrayn +1 
![151] ipack= jjx/ibytes*nrayelts*nrayn +1
![151] No need to multiply by mrfn, 
![151] because ifct1,2_(ipack16,mrfn) include mrfn,
![151] and ilowp(ipack,mrfn), iupp(ipack,mrfn) include mrfn.
![151] Saves considerable amount of memory!

![151] Reverted following, in favor of implementation in the full FOW 
![151] version of cql3d [BH+YuP170506].
![151] YuP:  Added fow_orbits.f to compilation !!!!!!!!!!!!!!
![151] It contains subroutines for Finite Orbit Width option.
![151] THE WORK is in PROGRESS. The routines are NOT used by default.
![151] Called from tdinitl:
![151] if(fow .eq. 'enabled') then
![151]   call fow_alloc ! Allocate arrays for orbit tracing and com_map.
![151]   call fow_setb ! Setup rectangular grid in (R,Z) for storing the
![151]              ! values of equilibrium m.field and its derivatives.
![151]              ! Setup Beq*(ir,iz), psieq(ir,iz), and derivatives
![151]              ! on the req(ir),zeq(iz) grid.
![151]              ! Store in common/Beq/ and common/BeqGrad/
![151]              ! Needed for finite-orbit-width calculations.
![151]   call com_map ! Setup 3D grid for storing a map 
![151]             ! (U,mu,Pfi)->Rmidplane 
![151]             ! In other words, setup a lookup table
![151]             ! which gives the midplane value(s) of R 
![151]             ! for FOW orbits ("leg's" major radius at midplane)
![151]             ! as a function of three indices corresponding to COM
![151] endif 


![150] YuP: In coefwti:  Added to prevent jchang=0 :
![150] if(jchang(k,l_).le.0) jchang(k,l_)=1  


![149] YuP: In ilut:  Re-defined cutlo=dsqrt(EPSILON(one)) 
![149] (the smallest number - machine dependent)


![148] YuP: Allocation of some arrays is adjusted to reduce memory load
![148] (changed from jxa to jx, iya for iy, ngena to ngen):
![148] cal, cbl, ccl, cdl, cel, cfl, eal, ebl, ...


![147] YuP-101122: nraya and nrayelta are not used anymore.
![147] Instead, nrayn and nrayelts are determined 
![147] in urfsetup by reading ray data files.
![147] In urfalloc:
![147] Added if(istat.eq.0) in front of call bcast() or ibcast()
![147] If istat=41 (not enough memory), cannot call bcast()
![147] because array is not allocated => results in Seg.Fault.
![147] The arrays ifct1_, ifct2_ are quite large 
![147] and may cause memory allocation problem.


![146] YuP: Corrections in tdchief. "bug" affected diagnostics output:
![146] Added
![146]     call coefstup(k) ! To define da...df coeffs, gon(i,j), etc
![146]     call coefmidv(da,temp1,1,vptb(1,lr_))
![146]     call coefmidv(db,temp1,2,vptb(1,lr_))
![146]     call coefmidv(dc,temp1,3,vptb(1,lr_))
![146]     call coefmidt(dd,temp1,1)
![146]     call coefmidt(de,temp1,2)
![146]     call coefmidt(df,temp1,3)
![146]     call coefwtj(k)
![146]     call coefwti(k)
![146] just before
![146]     call diagimpd(k) 


![145] Added option (netcdfshort=long_urf) to output all urf or 
![145] rdc coefficients at the last time step [BH100917].

![144] Added option to output specific current currv(u,r) and 
![144] rf power pwrrf(u,r) at each time step, rather than just
![144] on the last. Enabled by netcdfshort="long_jp"
![144] [BH100912].

![143] Fixed bug in urfalloc.f where insufficient space could
![143] by allocated for ilowp/iupp and ifct1_/ifct2_ for compressed
![143] urf ray data storage for 32 bit integer machines.   Evidently, 
![143] this usually did not cause a problem  [BH100903].

![142] Added some checking of memory allocation status in ainalloc
![142] and urfalloc.  However, this didnt properly catch a case of
![142] too large memory request, which led to a Segmentation fault.
![142] So, added print out of when allocation routines are entered
![142] and left (in ainalloc,urfalloc,eqalloc,sigalloc,vlfalloc,wpalloc,
![142] freyasou,impavnc0,urfbes,tdtraloc, and vlfsetup [BH100330].

![141] version='cql3d_merge_100829'
![141] Added multiple SXR and NPA synthetic diagnostic sightlines
![141] and starting postions to the netcdf file.  Further debugged
![141] NPA.  [BH100829].

![140] Fixed up NPA plotting in .ps file, and added NPA data to
![140] the .nc netCDF file. Added ennscal(1:npaproc) scale factors
![140] for corresponding density profiles [BH100815].

![138] Added NPA related namelist variables, npaproc, npa_process(),
![138] giving access to CX with boron+4 and electron recombination.
![138] Modified ennl/ennb from scalars to arrays (1:npaproc),
![138] enabling separate density profiles for related CX species.
![138] ennin modified from 1D to 2D array, ennin(1:njene,1:npaproc).
![138] rd_npa,thetd_npa,x_npa,z_npa modified from scalar to
![138] arrays (1:nv_npa), enabling separate specification of
![138] detector locations [while maintaining backwards compatibility.
![138] [BH100720].

![137] version="cql3d_merge_100608".  This is a major modification.
![137] It combines several separate branches of cql3d, and includes
![137] fully-implicit 3d radial transport (soln_method=it3drv), 
![137] URF routines applied to multiple general species(multiURF),
![137] full non-updown symmetric equilibria (eqsym=non),
![137] NPA diagnostics, deltar first-order finite-orbit-width
![137] shift (and output to .nc file).  
![137] The 1st order orbit shift is not yet integrated into 
![137] the calculation of RF diffusion coefficients in urf-routines 
![137] or into diagnostics). [Yup, BH, 100608].

![136] Modified method for calculation of NPA, removing
![136] use of Legendre polynomial expansion of distribution
![136] function (such as used with SXR), since CX cross-section
![136] is much simpler (assuming CX ions are much faster than
![136] the background neutrals.  Taking neutrals to have zero temp,
![136] then CX of FI to neutral occurs without change in energy of dirn.
![136] Removed m_npa from NPA related namelist.  [BH, 100521]

![135] nrstrt namelist variable functionality removed.  It was
![135] generally not used anymore, and can take on only default
![135] value, nrstrt=1. Reason for removal:  Was causing difficult
![135] logic in achiefn.f, related to new soln_method=it3dv and it3drv.
![135] [YuP and BH, 100517].

![134] Fixed code bomb for lrz=1,nmlstout='trnscrib', in which
![134] cqlinput was open in ain_transcribe, when trying to open it.
![134] Changes to achief1.f, frset.f.  [BH100517].

![133] Fixed an indexing problem with the equation scaling in tdtranspn
![133] and impavnc0. This problem arose for the new soln_method='it3drv'
![133] and 'it3dv' functionality.  For lbdry(k)='scale' or 'fixed', 
![133] the v=0 boundary point (j=1 index) gets a special treatment - 
![133] it is not included into the sparse matrix. In tdtranspn, 
![133] for j=1 & lbdry(k).ne."conserv", set radial elements to zero. 
![133] In impavnc0, the j=1 point  is not updated from solution 
![133] matrix (rhs);  update rhs=>f is performed from jstart=2.
![133] Also, the re-scaling of f 
![133] (call diagimpd -> calls diagscal -> renorm f if lbdry.eq."scale")
![133] is moved from impavnc0 to tdchief.
![133] Tdchief now has two loops in radial index.
![133] The first loop calls achiefn(0)-->impavnc0, 
![133] which defines matrix coefficients and finds new f. 
![133] For soln_method=it3drv, the solution is only found
![133] when the loop reaches the last radial index, 
![133] so it is important to postpone with diagnostics or
![133] re-scaling of f until the end of the first loop in radial index.
![133] The second loop re-scales f if needed, and computes diagnostics.
![133] These changes resulted in total current to be different by ~15% 
![133] in a MST radial transport (soln_method='it3drv')test case. 
![133] The value of current  obtained with lbdry(k)='scale' is now 
![133] different by that from  lbdry(k)='conserv' by less than 1%  
![133] (before corrections: 15%). The shapes of f now are also 
![133] almost same in runs with lbdry(k)='scale' and 'conserv'. 
![133] No more problem with v=0 point in f. [YuP and BH, 100517].

![132] De-equivalenced ca(:,:), cb(:,:), etc., storage from da(:,:),
![132] db(:,:), etc.  There was an error in re-scaling (call dscal) of
![132] ca, cb, etc. in cfpcoefn.f/cfpcoefr.f; not all of coefficients
![132] could have been properly re-scaled because of presence of zero
![132] index in ca...cf. Now the size of ca...cf is set to (1:iya,1:jxa), 
![132] no zero index. The difference between the corrected and the old
![132] 100217 version is ~3.3% in value of total current. This might have
![132] affected only soln_method='it3dv' and 'it3drv'.  Had no effect
![132] on a rdcmod='enabled', soln_method='direct' case [YuP, 100513].

![131] Multi-species QL diffusion (ngen.ge.2) capability was added
![131] (previously only general species 1 could be QL diffused).
![131] This was work in Sept'08 and Oct-Nov'09 by BH, under GA contract.
![131] This work merged into mainline cql3d_cswim_svn cql3d version,
![131] including debugging [Merge mostly YuP, with BH, Apr-May,2010].

![130] A fully-implicit 3d (2d-in-vel, 1d-in-radius) solution of
![130] the FP equation was introduced into CQL3D using conjugate
![130] gradient sparse-matrix solve techniques, including 
![130] "drop tolerance".  This uses SPARSEKIT2 (Yousef Saad, U Miss.)
![130] routines. [BH, see c[100]].
![130] Merged this code into mainline cql3d_cswim_svn cql3d version,
![130] including debugging [Merge mostly YuP, with BH, Feb-Mar,2010].


![129] Modified initial posn of soft xray detector from scalar to
![129] array dimensioned [1:nv].  This enables multiple detector
![129] positions, as needed by MST.  Mod by Mike Kaufman, UW,
![129] added to mainline cql3d, 100510.

![128] Fixed code bomb for lrz=1,nmlstout='trnscrib', in which
![128] cqlinput was open in ain_transcribe, when trying to open it.
![128] Changes to achief1.f, frset.f.  [BH100517].

![127] version='cql3d_100420'.
![127] Petrov upgraded some netCDF coding from netCDF2 and netCDF3.
![127] BH modified NPA routines to plot energy spectra.
![127] BH modified tdpltmne for profile plotting at n=1:
![127] Incorrect background species energy vs radius were
![127] being plotted for nstop=1 cases, as occur in CQL3D/AORSA
![127] coupled runs [BH100420].

![126] version="cql3d_100126"
![126] enescal/tescal/tiscal applied to all input profiles, static
![126] and time-dependent.  Useful for adjusting units, e.g.
![126] [BH100126].

![125] version="cql3d_100116"
![125] Added rdcmod="format2" read capability to rdc_multi.f, to
![125] read multiple du0u0 files for each radius from DC. [YuP, 100101].

![124] Located bug in tdtranspn.f which was causing unphysical bulge
![124] in tail of fe in EC test case, and added 3-radial-point smoothing
![124] of density profile and distn function to get working 
![124] soln_method='it3drv' (fully implicit 3D interative solve)
![124] [YuP and BH, Dec. 9, 2009].

![124a] Fixed bug in setting density profiles for multiple ion species
![124a] with same charge (e.g., D-T)  [BH091121].

![123] version="cql3d_091118".
![123] Reworked urfread_ read of ray data to accomdate possibility of
![123] added data (complex cnper from toray) [BH091116].

![122] Fixed bug to properly reset urfmod='disabled' in ainsetva, 
![122] in case no urf modules are are setup [BH091116].

![121] Re-arranged expressions for derivatives of gamma^m * alpha^-n,
![121] to increase accuracy, etc.  The code is somewhat more stable 
![121] for relativ=fully.  Can use mx=5, if higher-m coeffs at j<33 
![121] are zero-ed, see cfpleg.  But not as good as hoped - still a 
![121] noise at v~0 starts to develop.  The major part of cfpcoefr 
![121] should be re-written to resolve the problem at v~0, if
![121] in MeV range need be used.  Elsewise, the quasi-relativistic
![121] relativ='enabled' approximation is quite sufficient, see
![121] report CompX-2009-1_Fully-Rel.pdf  [Yuri Petrov, 091016].

![120] The Intel ifort compiler on viz.pppl.gov differed in compiling 
![120] a comparison between a real*8 variable and 0.0, so changed all
![120] .ne.0. and .eq.0. in the code to be comparisons with real*8
![120] zero=0.d0  (in about 35 source files). [BH090904].

![119] Found error in velocity theta flux expression in advnce.h,
![119] a "k" rather than "i" in indexing of distribution function.
![119] Affected plot of velocity-space flux vectors, e.g., efld vectors
![119] not purely parallel, and not much (perhaps nothing) else
![119] [BH090827].

![118] Found bounds check violation in choose(,) in fully relativistic
![118] FP collision coeff calculation, and fixed coding in comm.h, 
![118] micxinil.f and fpcoefr.f to agree with Franz\'s thesis
![118] [YUP+BH, 090809].

![117] Added additional conversions of old GRAFLIB plotting to PGPLOT
![117] library, removing all references to graflib [YUP+BH, 090807].

![116] Added pwrrf (rf power versus x=u/unorm) to netcdf file, modifying
![116] storage [BH090611].

![115] Modified urffflx to assign ray elements which are outside the
![115] LCFS to the outer radial bin inside the LCFS (urfmod="enabled").  
![115] This removed on out-of-bounds reference, resulting from the
![115] new GENRAY capability for ray tracing in the SOL [BH090602].

![114] Updated fully-implicit 3D interative solve from mainline
![114] cql3d_f90ptr_090514 version.  This is fully f90 pointered
![114] version rather than cray pointers, facilitating debugging
![114] [bobh, 090514].

![114a] Added option to output soft x-ray fluxes to the netcdf file at 
![114a] each time step (see nml softxry); previously only at first and 
![114a] last time step.
![114a] Similarly, an option was added to output the distribution function 
![114a] f and the QL coeff urfb (through nml netcdfshort) at each time
![114a] step.  Updated to version="cql3d_f90ptr_090514" [BH090514].

![113] Adding namelist rdc_upar_sign to reverse u_par order of rdc
![113] diffusion coeffs (for rdcmod.eq."enabled"), appropriate
![113] DC originated coeffs when eqdsk toroidal magnetic field is neg.
![113] Reverted scaling of rdc diffusion coeffs for cases with
![113] when coeffs originated on a grid with different normalization
![113] energy [BH090424].

![112] Fixing read a time-dependent generalized parabolic density/
![112] temperature/zeff/elecfld/vphi.  There were errors when an
![112] impurity was automatically added, showed up for a runaway
![112] knockon test case [BH090312].

![111] Fixing/adjusting plotting to remove glitches in 64-bit plotting
![111] at very low ordinate values (e.g. power 10**-43 watts).
![111] Added setup0 namelist variable lnwidth, specifying plotting
![111] line width (useful for publication) [BH090220].

![110] Add namelist &setup0 variable special_calls, which if set
![110] to "disabled" will avoid system calls (which are not enabled
![110] for some systems (e.g., franklin.nersc.gov, jaguar.ornl.gov).
![110] default="enabled" [BH081216].

![109] Switch entirely to f90 pointers rather than cray pointers.
![109] This makes debugging easier for gfortran/gdb (and probably
![109] debugging systems) BH081216.

![108] Generalized rdc_multi to input data velocity normalization
![108] different from vnorm (may be less, or greater), and added
![108] option to read file specifying prompt loss.  This is for
![108] coupling to the DC diffusion coeff calculation [BH081201].

![108] Dimensioned namelist variables iurfcoll and iurfl 1:nmodsa,
![108] in support of multiple ray type and multiple general species
![108] application of rf [BH081106].

![107] Fixed inadaquacy of calc of iprozeff = "parabola" or "spline"
![107] for multiple ion species with same bnumb plus ions being
![107] a general species [BH081031].

![106] Fixed bug (gfortran compiler?) wherein xpts(1) and rx(1), etc.,
![106] were offset by one real*8 memory position, contrary to equivalence
![106] statement in freya.f, by adding dummy integer in frcomm.h
![106] after npts.  The problem was causing freyasou bomb [BH,081029].


![107] Zeroed vflux before summing in diagimpd [BH081125].

![106] Modified tdsxray.f/tdsxr0.f so that sightline distance from
![106] the detector to the plasma is increased.  Otherwise, the
![106] calculated sightlines may not reaching the plasma [BH081106].

![105] Checked for possibility that lsa.lt.lrza.  Has not usually been
![105] the case in the past, but can cause memory overwrite [BH081106].

![104] Fixed inadaquacy of calc of iprozeff = "parabola" or "spline"
![104] for multiple ion species with same bnumb plus ions being
![104] a general species [BH081031].

![105] version="cql3d_biptr-mpi_080925"
![105] Brought mainline cql3d version up to date w.r.t. Aug 19-21,2006
![105] changes below.
![105] Added changes for multispecies (ngen.ge.2) urf and rdc
![105] QL diffusion. Add namelist .... [BH080918]

![104] version="cql3d_biptr-mpi_080909"
![104] Removed calculation of z00(i,j) [rhs of main equation] from
![104] advnce.h statement functions, putting the z00 directly
![104] into impavnc0 (impavnc) as an internal procedure function.
![104] This aided in debugging a problem for nso=1, when source
![104] parameters where set to zero.  Several traps agaist divide
![104] by zero under these conditions were added.
![104] Removed equivalences of dxp5/dxm5 and exp5/exm5, and
![104] set the variables in micxinit.f [BH,080909].

![127] version="cql3d_3d-implicit_f90ptr_080303".
![127] Combined fully-implict 3D eqn solve and mainline cql3d
![127] [bobh, 080303]

![103] Added option to restart variable, nlrestrt="ncdfdist", to
![103] restart using netcdf saved distribution, as a higher accuracy
![103] alternative to restart from distfunc text file [bobh, 080203].

![102] version="cql3d_biptr-mpi_080125".
![102] cql3d running on 64-bit machines such as Intel core 2
![102] quad processor, using gfortran64 and Lahey-Fujits lf64.
![102] Added SAVE statements to subprograms with data statements,
![102] removing some bugs related to differences in compilers.
![102] Also fixed bug in asor, giving false turnon
![102] of anaytic particle source [error introduced
![102] in [98], below.  Added new D1MACH() function, based on
![102] f90 calls.  Previous version had hex definitions which
![102] were probably incorrect for 64 bit machines[bobh, 080125].

 
![101] Made a few modifications, principally in tdxinitl.f and profiles.f     
![101] to enable iprozeff='parabola' or 'spline' to work with two ion          
![101] species with same bnumb() (e.g., H+,D+, and/or T+).                     
![101] Did a ngen=2 (D+,H+) test run (with ngena=2) in                 
![101] /home/bobh/cql3d/aorsa/D3D_test_case/122080.0/116MHZ/8th.1/tmp_ngen2   
![101] Density H+ is 10**-4 of D+.                                             
![101] Uses makefile_lf95.  Results differ in 4th sig fig from previous       
![101] ngen=1 D+ run, as expected.                                             
![101] Adjusted ainsetva.f to account for iprozeff=1 option for
![101] calc of ion densities from ene/zeff so treat different ions
![101] with same charge (e.g., H+,D+ and/or T+). [bobh, 060819].
![101]                                                                 
![101] Execution time increased from approx 8 minutes to 20, but               
![101] addition of the 2nd general species.  There is no QL diffusion          
![101] yet on the 2nd species.    [bobh, Aug 19-21, 2006]                      

![100] First results from new soln_method="it3drv" option using
![100] sparskit2 to iteratively solve full 3d (2V-1R)
![100] implicit  cql3d equations [bobh, 071108]
![100] Dimensioned wk ==> wk(*) in zcunix.f: coeffs to prevent
![100] tripping of array length checker [bobh, 070805].
![100] New iterative sparse-matrix solution capabilities for solving 
![100] the basic FP matrix equation are introduced, as specified
![100] by new namelist var soln_method, invoking SPARSKIT2 routines.  
![100] Additional new namelist vars are droptol and lfil. 
![100] The distinction between distributions for the velocity and
![100] radial steps has been removed, according to setting of the
![100] internal variable ipacktr=0 in tdtrmuy.f.  The resulting
![100] transp="enabled" solutions did not change substantially.
![100] This is prepatory to solution of fully-implicit 2D-1r 3D
![100] equations (with transport) by iterative sparse-matrix methods
![100] [bobh, 070424].

![99] Changed nlwritf and nlrestrt from logical to character*8,
![99] for more flexibility in specifying restart settings.
![99] Will need to check use of past cqlinput files. [bobh, 070414].

![98] Broke setup namelist up into two: setup0 and setup.
![98] This was necessary for SWIM IPS, involving namelist writes,
![98] and reads.  Backwards compatibility is maintained by checking
![98] the cqlinput file for existance of &setup0:  if not, then
![98] cqlinput is rewritten with first occurance of &setup changed
![98] to &setup0.  Initial (old) cqlinput is restored at the
![98] end of the run.  [bobh, 070414]

![97] Added runaway electron related variables (denra,curra,
![97] ucrit,knockon,eoe0,srckotot), also wpar,wperp to the 
![97] netcdf output file. [bobh, 070407]

![96] The namelist input system was restructured to facilitate
![96  coupling, and backward compatibility after further cql3d
![96] upgrades, with the SWIM Integrated Plasma Simulator (IPS):
![96] namelists, namelist type and dimensions (i.e. declarations),
![96] and subroutines setting defaults have been separated out from
![96] other code variables . Consequently, the IPS interface uses
![96] files from the mainline cql3d distributions:  the name.h,
![96] name_decl.h,frname.h,frname_decl.h include files, and sets 
![96] defaults using the subroutines in aindfpa.f,aindlft.f,eqindflt.f,
![96] urfindfl,frinitl.f (of the same root names) [BH070310].

![95] A 'if ...write(*,*)' [never executed] was added to zcunix.f
![95] to get around compiler problem on viz.pppl.gov SGI machine.

![94] Modified eqtopeol (reads Culham TOPEOL file in lieu of eqdsk)
![94] to integrate poloidal B-fields from the equatorial plane
![94] rather than minimum -Z plane.   This gets around a problem
![94] of inexact B-field/singularities outside last closed flux surface
![94] [Part of benchmarking effort with Saveliev] [bobh, 070116].

![93a] Vickie Lynch, ORNL found an ancient bug in freyasou (setting
![93a] of array with unset index) which could cause an overwrite,
![93a] but seems to have not usually been a problem. [bobh, 060824].

![93] Added Vincent Tang (MIT, next LLNL) coding for a passive NPA
![93] synthetic diagnostic.  This is a modification of the SXR
![93] synthetic diagnostic.  In future work, it is intended to
![93] convert active NPA synthetic diagnostic coding developed
![93] in Matlab by Vincent to a fortran module with cql3d. 
![93] [bobh, 060622].

![92] makefile_lf95_ptr is system developed by Nikolay Ershov     
![92] to produce a Fortran 90 pointered version of the code,      
![92] as alternate to the standard Cray pointered version obtained
![92] with makefile_lf95.  Comparison of these two makefiles shows the    
![92] only differences are in four lines referring to doptr.py    
![92] and tmpptr.f.                                            
![92] The source code modifications are carried out with the      
![92] python script doptr.py in ./ptr, and additional files       
![92] in that subdirectory.
![92] Future additions to code memory should be added in both the
![92] mainline cray pointered version, and in the f90 pointer mods in
![92] ./ptr.
![92] Please code consistent with this scheme for modification in code
![92] storage, so that the two pointering systems are carried forward.
![92] [N. Ershov; bobh, May\'06].     

![91] Removed limitation in urfmod and vlfmod options of no more 
![91] than 3 rf wave types, or, 3 harmonics for 1 wave type.  Now
![91] can have multi-wave types with multi-harmonics.  The limit
![91] on number of wave-types*harmonics is set by parameter 
![91] nmodsa [bobh, 060319].

![90] Changed passing of the distribution function to functions
![90] bsl and bsu by common block rather than as function argument.
![90] This GREATLY reduced execution time for f90 pointered version
![90] cql3d (Lahey-Fujitsu f95 compiler, V6.0)  [Ershov, 060316].

![89] Introduced parallelization of impavnc/impanv0 solver of the
![89] 2D bounce-averaged equations on each flux surface.
![89] See MPIreport.pdf in cql3d_manual directory.  Parallelization
![89] is obtained with CMPI comments/insertion points in source,
![89] activated by makefile_mpi, using python scripts in ./mpi
![89] [Ershov, 060114].

![88] Fixed bug in bsu.f which may effect bootstrap current calc
![88] at outer radius, rho(lrzmax) [bobh, 051101].
!
![87] Added NPA (neutral particle analyzer) diagnostic skeleton code
![87] based on first general ion species and added neutral profile. 
![87] See npa_diag and ipronn related namelist [Vincent Tang (MIT), 
![87] bobh, 051007].
!
![86] Added additional phenomenological radial diffusion drr from
![86] radius 0. out to normalized radius difus_rshape(8), to 
![86] simulate sawteeth effects [bobh, 050921].
!
![85] Substantial bug fix in radial dependence of velocity independent
![85] radial diffusion coefficient.  Radial dependence specified by
![85] difus_rshape() is implemented.  Previously, radial dependence
![85] was constant with r, for constant in velocity space cases 
![85] [bobh, 050921].
!
![84] Added namelist var difus_type with possible values of "specify"
![84] (default, giving older methods of specifying vel and r variation
![84] of drr, and new "neo_smpl" giving simple, velocity-independent
![84] neoclassical, bannana regime drr, and "neo_plus" which adds
![84] "specify" type drr to the "neo_smpl" drr. [bobh, 050919]
!
![83] Fixed bug in tdtravct.f. Target density for radial convection
![83] was not being set for colmodl.ne.0, preventing radial transport
![83] runs with usual colmodl=3.  Error introduced with colmodl=0,
![83] radial transport update in March, 2002. [bobh, 050913].
!
![82] Fixed bug for bsign=-1 cases (neg tor fld) which was preventing
![82] addition of salphal additional damping for iurfl="enabled".
![82] Bug existed for about last year [bobh, 050911].
!
![81] Small k_parallel-width ray elements were giving zero damping
![81] due to inaccuracy in interpolating the diffusion coefficient
![81] to neighboring velocity-space grid points.  This was improved
![81] by increasing storage of ifct1/ifct2 from 8-bit words to 
![81] 16-bit words [bobh, 050812].
!
![80] Fix bug: Sauter\'s unicity of f at j=1 (v=0) in impavcn0 was
![80] being  applied in lbdry(1).ne."conserv" cases ("fixed"
![80] and "scale"), which already applied unicity.  An out of bounds
![80] reference was generated in the FP coeff matrix. Effects
![80] are unknown [bobh, 050804].
!
![79] Moved collisional and linear damping calc out of urfdamp1
![79] and urfdamp2, to new subroutine urfdampa.  This fixes a bug in 
![79] which the additional damping was not calculated for ray elements
![79] for which the QL damping was not calculated [bobh, 050426].
!
![78] Changed flag indicating output of diagnostic damping rate
![78] from iurfl="damp_out" to iurfcoll="damp_out", thus enabling
![78] simultaneous input of additional damping in salphal and
![78] output of diagnostic damping in salphac [bobh, 050423].
!
![77] Included gyro-radius with banana width in lossmode()=simplban
![77] scrape-off loss model [bobh, 050217].
!
![76] Fix bug: Added bsign1() in urfread_.f to account for the
![76] different sign conventions in genray and toray when bsign
![76] .lt.0. This bug, introduced in [74] could give zero cyclotron
![76] damping for negative toroidal field in the eqdsk [bobh, 050206].
!
![76] rdcmod modification to read in externally computed RF
![76] diffusion coefficients for an array of flux surfaces
![76] [bobh, 041111].
!
![75] Fix bug: fixed eqdskin functionality so specification of full
![75] path for the eqdsk works, not just local eqdsk [bobh, 040902].
!
![74] Version designated as cvs_cql3d_v2-00_040820
![74] Added netcdf read of standard ray data input files, as
![74] alternative to reading text files. Modified netcdfrf.f. 
![74] Added namelist variables rffile(1:3),rfread
![74] Names of input netcdf files can be specified through rffile.
![74] Bug fixed: sdpwri was output to .nc file in netcdfrf.f,
![74] but was not dimensioned.  Unknown effect.
![73] [bobh, 040814].
!
![73] Added namelist variable nmlstout, to control namelist o/p to 
![73] stdout.  Various slight changes in write(*,*).
![73] Added vol int of FSA power densities entr ==> entrintr(,) and
![73] put result in .nc file at each time step.
![73] Fixed bug in zmincon,zmaxcon determination.
![73] [bobh, 040326].
!
![72] Added netcdf output giving fusion/neutron rates and powers,
![72] and RF damping due to salphal.
![72] Added time-dependent input power in urf rays (nurftime.gt.0).
![72] [bobh, 030618].
!
![71] Version designated cvs_cql3d_v1-12_030115. 
![71] CVS repository changed /usr/local/cvsroot.
![71] Added netcdf output giving ratio of theoretical electrical 
![71] conductivity by Connor and by Kim-Hirshman/Sigmar  [bobh, 021121].
![71] Added additional output to _"flux_" netcdf files [bobh, 030115].
!
![70] Added variable vlh_karney to enable Karney and Fisch
![70] u-variation of LH Dql (PF 28, 116 (1985)).  [bobh, 021118].
!
![69] Version designated cvs_cql3d_v1-11_021028.
![69] Added new capability to output vel space fluxes in x,theta-coords
![69] (see netcdfvec.. variables in cqlinput_help).
![69] Added netcdf output of specific current dens j(u,species,radius).
![69] [bobh, 021025]
!
![68] Version designated cvs_cql3d_v1-10_020914.
![68] Added (much) more accurate calculation of qsafety.  I presume
![68] that in a bounce-avg code, this won\'t have significant effects.
![68] It does affect zmax. [bobh, 020914]
!
![67] Version designated cvs_cql3d_v1-9_020607.
![67] Added some symmetrization options (see eqsym namelist variable).
![67] Added eqdskin namelist.  Added netCDF variables [bobh, 020607].
!
![66] Version designated cvs_cql3d_v1-8_020101.
![66] Added most remaining data in screen output  into netcdf file.
![66] Fixed some dimensioning bugs in netcdf o/p of flux vectors.
![66] Adjusted vector plots of vel space fluxes for PGPLOT.
![66] [bobh, 020101]
!
![65] Added parallel current and resistivities to netcdf file.
![65] [bobh, 010814]
!
![64] Added simple banana loss model removing particles with banana
![64] width greater than distance to the plasma edge [bobh, 011125].
!
![63] Added parallel current and resistivities to netcdf file.
![63] [bobh, 010814]
!
![62] Fixed bug in restvty.f: efswtchn.eq.ne0_hh ==> ne0_hh [bh,010812]
!
![61] Added X-ray detector specs and calculated fluxes to
![61] the netcdf output file (when softxry.en."disabled" [bobh,010526]
!
![60] Added a globally convergent Newton-Raphson iteration scheme
![60] to find radial transport velocities which maintain a given
![60] target radial density profile [bobh, 010426].
!
![59] Added a general radial and velocity-space profiles for
![59] the radial diffusion coefficient d_rr [bobh, 010420].
!
![58] Version designated cvs_BH20010411. [bobh, 010411]
!
![57] Output files now all pre-pended with contents of the
![57] character*40 variable, mnemonic. [bobh, 010411]
!
![56] Added netCDF output facility for velocity-space flux vectors,
![56] through netcdfmain.  This approach could be expanded to
![56] additional output, similar to pltmain facility. [bobh,010411]
!
![55] Switched to CVS code maintenance system at bobh.compxco.com.
![55] [bobh, 010319]
!
![54] Fixed pitch angle references in vlf.f. Found that there is
![54] an incompatibility between cosz,sinz,tanz setting and the
![54] maximum pitch angle given by imax, with taunew="enabled".
![54] Changed taunew default to "disabled".  This needs more attention.
![54] [bobh, 010316]
!
![53] Added target current calculation from eqdsk.  Fixed elecfld
![53] iteration to achieve target current, for multiflux-surface
![53] runs. [bobh, 010301].
!
![52] Added e-e Bremsstrahlung to existing e-i, for Xray calc.
![52] Debugged Xray calc [bobh, 001200].
!
![51] Incorporated numerical bootstrap current calc  [bobh, 990731].
![51] Added numerical calc of Hirshman'88 and Sauter'99 analytic
![51] bootstrap current, for comparison  [bobh, 990825].
!
![50] Added netCDF o/p for main plasma parameters and RF data
![50] [bobh, 990601]
!
![49] Added quasineutrality calc of E_ll to cqlpmod="enabled" model
![49] [bobh, 9904]
!
![48] Added relativistic QL diffusion coeffs (vlfmod="enabled") to 
![48] cqlpmod="enabled" (1D-in-dist-along-B,2V) model.  [bobh, 990402].
!
![47] Major modifications to run code on both 64- and 32-bit
![47] machines (machinea=1 or 2), and replacing GRAFLIB graphics
![47] with more public domain PGPLOT.  Code compiles
![47] on PC with Absoft, DEC Alpha with DEC For, and J90[bobh, 990501].
!
![46] Added RFP upgrades: radcoord for new radial coords, calc and
![46] print pol and tor current, -(epsi) for psilim.lt.psimax,
![46] minor mods for RFP equilibria (removing a kluge for
![46] non-monotonicity of B(s).   [bobh, 980919].
!
![45] Corrected bounce-average in sourceko.f, to get agreement
![45] with MNR.  Multi-flux surface enabled. bobh,980501
!
![44] Corrected cross-section for KO operator to agree with Heitler
![44] bobh, 980417
!
![43] Fixed bug of incomplete copy of complex*16 numbers in urfread.f
![43] bobh,980411
!
![42] Update vth and sptzr at each time step according to current
![42] temp.  Option for calc of elecfld using neoclassical resistivity
![42] plus runaway electron current. Fixed pdenra in aclear[bobh, 970922].
!
![41] Changed out impavnc.f and impavnc0.f for Olivier Sauter\'s
![41] new version which uses standard LAPACK routines sgbtrf sgbtrs, 
![41] rather than zband1 (which might have error under f90)[bobh,970901].
!
![40] Added time-dep profile option for (1+cos(pi*t))-variation 
![40] bobh, 970714].
!
![39] Added calc of ko source using pitch angle average primary
![39] distribution: fsa [970701], local in pol. angle [bobh, 970706].
!
![38] Many little changes in plotting, sourceko.f soucrit.f[bobh,970417].
!
![37] Restored equivalence of temp2 and fpn(xpar,xperp), added
![37] subregion for ploting distributions and fluxes [bobh,970331].
!
![36] Added coding to warn of potential overwrite by tamt1,tamt2,
![36] when mxa is too large [bobh, 970312].
!
![35] Added improved, fully numerical calc of tau and dtau,giving uniform  
![35] poloidal density for isotropic distributions [bobh, 970212].
!
![34] Added pol. angle variation for parallel distribution appearing 
![34] in knock-on source.  Using Rosenbluth formula [sc,bobh, 970122].
!
![33] Upgrade in cfpcoefn for high Zeff contrib. to Cee [bobh, 961217]
!
![32] Added plots of runaway current and density [SC].
![32] Added write and read of preload distribution [bobh, 960919].
!
![31] Added new "conservative" treatment of knock-on source, with new
![31] subroutines for reduced distn fle, and sourceko. (bobh, 960820).
!
![30] (Temporarily) increased max resolution of calc of reduced distn in 
![30] pltprppr to 201x10001, as fix for runaway mesh resolution problem 
![30] (bh 960724 & 0801).  May use unnecessary storage. Added xprpmax.
!{30} Gave fpn storage.  Added jpxy.le.jpxya, ipxy.le.ipxya to namelist.
![30] Don\'t try anything that uses xtail, xhead, ytail, yhead!
!
![29] Changed Bremsstrahlung radiation evaluation in lossegy. An enhancement
![29] factor is multiplied in the highest energy range to force
![29] containment of runaway electrons. Also changed the limits for
![29] energy for different formulas so that transition is smooth. (SCC960625)
!
![28] Introduced array of time step settings. (bobh 960620).
!
![27] Further changes and additions relevant to electron runaway problem:
![27]  renamed efield.f to coefefld.f, restinit.f to efield.f and added
![27]  control algorithm to calc. elec. fld. to give specified target
![27]  current.  Added bremmstahlung reaction force on electrons, and
![27]  enegy dependent Coulomb logarithm.  Also additional preloading
![27]  distributions (in finit.f), time-dependent profiles (new
![27]  subroutines in  profiles.f), plotting of coefficients.
![27]  Also, knock-on electron source (march 96). (bobh, circa 960521)
!
![26] Omitted damping in npar1*npar2.lt.0. case, in urfpack. This situation
![26]  can occur for near perp. injection, but needs more work on logic.
![26]  (bobh 960421).
!
![25] Simple search for accurate magnetic axis rather than Newton
![25]  iteration, in eqorbit.f (bobh, 960309)
!
![24] Minor change in tdoutput.f of synchrotron o/p.  Not sure correct
![24]  yet (bh 960305)
!
![23] Took into acct. like-like collision factor in neutron rates
![23]  for D-D particles (i.e., 0.5). Changed fusion cross-sections to 
![23]  Bosch & Hale, Nucl Fus, \'92. Benchmarks with ONETWO (bh 950307).
!
![22] Fixed over-write of nrayelt0, for nharms.gt.1, in urfdamp0.
![22] Extended background interpolation for freya from lrz to lrzmax.
![22] (bh 950221)
!
![21] Added toroidal rotation effect for ions into FREYA (bh 950109).
!
![20] Added calculation and plots of fusion power (bh 950102).
!
![19] Added zeff profile namelist input, thru zeffin (bh 941222).
!
![18] Added possibility to taper rf diffusion coefficient over
![18] the last 10 velocity points of the mesh (ineg=trunc_d), to 
![18] assist in obtaining solutions in ion/icrf problems where there
![18] is otherwise an rf-induced ion runaway (bh).

![17] Added pure perpendicular rf diffusion in vlh module (bh).

![16] Added plot of the theta-averaged distribution function (scchiu).

![15] Small change to impavnc0.f: In the case that lbdry=fixed or
![15] scale, lbdry0=enabled, j=1, i=iy, nvar set =1 rather than inew=1.
![15] This cleans up the "coefficient" matrix somewhat, possibly
![15] removes a potential error, but had no effect on test cases
![15] lrz=1, test case with EC.  (bh 94/08/16).
!
![14] Added vlf.. routines.  These are single flux surface versions
![14] of the urf... routines, and permit study of LH, FW, and multiple
![14] cyclotron harmonic effects on electrons or ions, specifying
![14] wave parameters and region of flux surface for QL diffusion,
![14] through simple namelist input.   (bh 94/06/21).
!
![13] Generalized single flux surface LH QL diffusion model (vlh...)
![13] to multiple resonance regions.  Also more deeply pointered
![13] the urfb, urfc,urfe,urff, so save significant memory. 
![13] Started new file (notes.h) keeping information on variables,
![13] etc., on the code.   (bh 94/06/21).
!
![12] Combining Sauter version CQLP_B_940425 and Harvey version to that
![12] below date. Indenting done of bh code according to [10] below, using
![12] Sauter-s emacs indenting setup files.  (bh 94/05/07).
!
![11] Change bcast(ca,6*..) and sscal(ca,..) in cfpcoefn and cfpcoefr
![11] to bcast(ca,3*iyjxua,..), bcast(cd,3*iyujxa,..) and so on, as
![11] cd(1,1) is not after cc(iya,jxa+1) in memory => bhos_940425
![11] This version is copied to CQLP_B_940425 and is the start for
![11] B versions of CQLP/CQL3D codes. (os 94/04/25)
!
![10]  Indented whole code according to Emacs-Fortran rules, with
![10]  2 columns shift after loops, ifs and continuation lines.
![10]  Had to insert new continue statements, as one should not have
![10]  2 loops using same label. It is also recommended not to use
![10]  "1" as do loops label. The default continuation character has been
![10]  set to "+". (os) => new version 940304
!
![9]  Fix gftp defn. in urfdamp2.f => 5% change in powurf. bh, 940301
!
![8]  Changes in urfpack and urfedge correcting treatment of 
![8]  weighting of edge of resonance region. bh, 940227.
![8]  Added namelist variable wdscale.
![8]  Should check nharm=0 and inside of resonance layer cases further.
!
![7]  Synthesized Sauter and Harvey versions incorporating following
![7]  code changes [6] into CQL3D version with FREYA, ion-cyclotron
![7]  damping, and multiple cyclotron damping on electrons or ions.
![7]  Results stored in CFS /u650/940224a/bh940224.tar  bh, 940224.
!
![6]  CQLP is built starting from the CQL3D E_930209 version, then
![6]  several new versions  followed. The last one before O.Sauter
![6]  left GA is the 930624 version. Then the version used for the
![6]  4th PET93 workshop  paper is the 930910 version saved by O.Sauter
![6]  on cfs in CQLP_A931103.taZ.         
![6]  Now this version has been cleaned up, and dcofleg modified to
![6]  represents the contribution between i-1/2 and i+1/2, thus f is
![6]  not interpolated at i+1/2 as before (in particular in cfpcoefn). 
![6]  This version, 940128, is given to B. Harvey.  Olivier S. 940128.

![5]  Fixed us some storage problems:  rovsp re-dimensioned, 
![5]  lndum adjusted accordingly, and eqdumptr put in common.
![5]  940215, BobH

![4]  Changed sign of n_par, if partner="lsc". BH 930131.

![3]  Modified prppr.f slightly to get rid of minor gliches
![3]  in interpolation to f(xpar,xperp) at emax. BH 940130

![2]  Fixed ainspec.f so it properly calculates maxwellian
![2]  species designators kionm(nionm) BH 930130

![1]  Added the namelist variables urfrstrt and urfwrray and
![1]  made associated changes in urf routines. Bob H. 940127

![0]  EQRHOPSI.F :eqsource="tsc" eqpsimin/max set as for eqdsk.
![0]  Began this record file.   (BobH, 931222)

![-1] The CQL3D code has been under development since about 1985
![-1] by Harvey and McCoy, as a 3D 2v-1r (v0,theta0,rho) FP
![-1] code based on (1) the cql 2D-in-v Kerbel and McCoy 
![-1] FP collision code at each flux surface, plus (2)
![-1] an added radial variable enabling accounting for rf
![-1] quasilinear of electrons and ions using ray tracing
![-1] input data, determination of self-consistent distributions
![-1] and ray energy damping, developed by Harvey and McCoy 
![-1] [First reported for LH heating and current drive in a
![-1] combined IAEA paper by Soeldner et al., Washington, D.C., 1990.
![-1] Radial dependent synthetic diagnostics have been added, such
![-1] as Bremsstrahlung xray emission.  ECE emission is calculated
![-1] from CQL3D distributions, in the HORACE code.
![-1] Diffusive radial transport and radial pinch is solved by
![-1] an implicit, alternating direction differencing scheme.
![-1] See http://www.compxco.com/cql3d_manual_110218.pdf for code
![-1] status in 1992.


