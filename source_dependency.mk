# Source line to be included in makefiles
# YuP[2019-05-31] frhexdrv.f is removed - not used

SOURCES  =   param.f90 cqlcomm.f90 r8subs.f90 advnce.f90 pltmain.f90 \
	pltdf.f90 eqrhopsi.f90 equilib.f90 \
	impavnc0.f90  frplteq.f90 tdeqdsk.f90  \
	bsl.f90  bsu.f90 \
	a_cqlp.f90  abchief.f90  achief1.f90  achiefn.f90  aclear.f90 ainalloc.f90  \
	aindflt.f90  aindflt1.f90 aindfpa.f90  aingeom.f90  ainitial.f90   \
	ainpla.f90  ainplt.f90  ainpltpa.f90   \
	ainsetpa.f90  ainsetva.f90  ainspec.f90  ainvnorm.f90 ampfar.f90 \
	aminmx.f90 bavdens.f90  bavgmax.f90  baviorbt.f90   \
	bcast.f90 cfpcoefc.f90  cfpcoefn.f90  \
	cfpcoefr.f90  cfpgamma.f90  cfpleg.f90     \
	cfpmodbe.f90  cfpsymt.f90  coefefad.f90  coefefld.f90  coefegad.f90  \
	coeffpad.f90  coefload.f90  coefmidt.f90  coefmidv.f90   \
	coefrfad.f90  coefstup.f90  coefsyad.f90  coefwti.f90  \
	coefwtj.f90  cqlconf.f90 diag.f90  diagcfac.f90   \
	diagdens.f90  diagdenz.f90  diagentr.f90  diagescl.f90  diaggnde.f90  \
	diaggnde2.f90  diagimpd.f90  diagscal.f90  diagwrng.f90  diagxswt.f90   \
	diagxswx.f90  dsk_gr.f90 dskout.f90  efield.f90  eflditer.f90 eqalloc.f90   \
	eqcoord.f90  eqelpse.f90  eqflxavg.f90  eqfn.f90  \
	eqfndpsi.f90  eqfninv.f90  eqfpsi.f90   \
	eqindflt.f90  eqinitl.f90  eqjac.f90  eqonovrp.f90  \
	eqorbit.f90  eqrhs.f90   \
	eqtopeol.f90  eqvolpsi.f90  eqwrng.f90  \
	esefld.f90  exlin.f90  exsweep.f90  exsweept.f90   \
	exsweepx.f90  finit.f90  firstdrv.f90  fle.f90 flxfn.f90  \
	freya.f90  freyasou.f90  \
	frinitl.f90  frinitz.f90  frnbdep2.f90  frnfreya.f90  \
	frset.f90 frsmooth.f90  frsplft.f90   \
	frstup.f90  frsubs.f  frsuppor.f90  frwrong.f90  \
	ilut.f90  impchk.f90  impnorm.f90  \
	lookup.f90  losscone.f90  lossegy.f90  lossorbm.f90  \
	losstor.f90  micgetr.f90   \
	micgmbnd.f90  micgnbnd.f90  micxinil.f90  micxinim.f90  micxinit.f90  \
	micxiniz.f90  \
	netcdfrf.f  netcdfrw2.f \
	ntdstore.f90   \
	ntloop.f90  pack21.f90 pltcycl.f90 \
	pltdnz.f90  pltelec.f90  pltendn.f90   \
	pltfvsv.f90  pltinit.f90  pltlosc.f90  \
	pltpower.f90  \
	pltrun.f90  pltvec.f90  pltvectr.f90  \
	pltvflux.f90  profaxis.f90  profiles.f90 prppr.f90   \
	prpprctr.f90  psif.f90 r8lsode.f  \
	rdc_multi.f90 rdc_bplt.f90 restcon.f90 resthks.f90  \
	restvty.f90  rf.f90  sigalloc.f90  siggy.f90  sigmax.f90  sigsetup.f90  \
	sigv5d.f90 sigfn.f90 sigie.f90  sigmaxwl.f90   sigv.f90  \
	soucrit.f90 sounorm.f90  soup.f90  soup0.f90  \
	sourc0.f90  sourcee.f90  sourcef.f90  sourceko.f90 sourcpwr.f90   \
	synchrad.f90  tdbootst.f90 tdboothi.f90  tdchief.f90  tddiag.f90  \
	tddsig.f90  tdfinterp.f90  tdinitl.f90  tdinlegw.f90  \
	tdinterp.f90  tdnflxs.f90  tdnpa.f90  tdnpadiag.f90  tdnpa0.f90  \
	tdnpacxcs.f90  tdnpalam.f90  tdnpabscs.f90  tdoutput.f90   \
	tdplteq.f90  tdpltjop.f90  tdpltmne.f90  tdpro.f90  \
	tdreadf.f90  tdrmshst.f90  tdsetnpa.f90  tdsetsxr.f90   \
	tdstin.f90  tdsxr.f90  tdsxr0.f90  tdsxray.f90  tdsxrplt.f90  \
	tdtloop.f90  tdtoarad.f90   tdtranspn.f90  \
	tdtoaray.f90  tdtraloc.f90  tdtransp.f90  tdtravct.f90 \
	tdtrchk.f90  tdtrchkd.f90  tdtrcon.f90   \
	tdtrdfus.f90  tdtrfcop.f90  tdtrflg.f90  tdtrflx.f90  \
	tdtrmuy.f90  tdtrrsou.f90  tdtrrtov.f90   \
	tdtrrtov2.f90  tdtrsavf.f90  tdtrsym.f90  tdtrvint.f90  \
	tdtrvsou.f90  tdtrvtor.f90  tdtrvtor2.f90   \
	tdtrvtor3.f90  tdtrwtl.f90  tdtry.f90  tdtscinp.f90  tdtscout.f90  \
	tdwrng.f90  tdwritef.f90  tdxin13d.f90   \
	tdxin23d.f90  tdxin33d.f90  tdxinitl.f90  urfalloc.f90  \
	urfavg.f90  urfb0.f90  urfbes.f90   \
	urfbplt.f90  urfchief.f90  urfdamp0.f90  urfdamp1.f90  \
	urfdamp2.f90  urfdampa.f90 urfdout.f90  urfedge.f90   \
	urffflx.f90  urfindfl.f90  urfinitl.f90  urfmidv.f90  \
	urfpack.f90  urfpackm.f90  urfrays.f90  urfread.f90   \
	urfread_.f90  urfsetup.f90  urfwrite.f90  urfwrite_.f90  \
	urfwrong.f90  urfwr0.f90  urfwr0c.f90 \
	vlf.f90 vlfalloc.f90 vlfbplt.f90 vlfsetup.f90 vlh.f90  vlhbplt.f90 vlhd.f90  \
	wpalloc.f90  wparsou.f90  wpavg.f90  wpbdry.f90  \
	wpcheck.f90  wpchgdy.f90 wpcthta.f90  wpelecf.f90   \
	wpinitl.f90  wploweq.f90  wpsavf.f90  wptrafx.f90  wptramu.f90  \
	wptrmuy.f90  wpvptb.f90  wpwrng.f90  wpmshchk.f90 \
	zblock.f90  zcunix.f90  zfreya.f   znonsym.f90
