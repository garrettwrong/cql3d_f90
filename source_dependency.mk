# Source line to be included in makefiles
# YuP[2019-05-31] frhexdrv.F is removed - not used

SOURCES  =   param.F90 cqlcomm.F90 r8subs.F90 advnce.F90 pltmain.F90 \
	pltdf.F90 eqrhopsi.F90 equilib.F90 \
	impavnc0.F90  frplteq.F90 tdeqdsk.F90  \
	bsl.F90  bsu.F90 \
	a_cqlp.F90  abchief.F90  achief1.F90  achiefn.F90  aclear.F90 ainalloc.F90  \
	aindflt1.F90 aindfpa.F90  aingeom.F90  ainitial.F90   \
	ainpla.F90  ainplt.F90  ainpltpa.F90   \
	ainsetpa.F90  ainsetva.F90  ainspec.F90  ainvnorm.F90 ampfar.F90 \
	aminmx.F90 bavdens.F90  bavgmax.F90  baviorbt.F90   \
	bcast.F90 cfpcoefc.F90  cfpcoefn.F90  \
	cfpcoefr.F90  cfpgamma.F90  cfpleg.F90     \
	cfpmodbe.F90  cfpsymt.F90  coefefad.F90  coefefld.F90  coefegad.F90  \
	coeffpad.F90  coefload.F90  coefmidt.F90  coefmidv.F90   \
	coefrfad.F90  coefstup.F90  coefsyad.F90  coefwti.F90  \
	coefwtj.F90  cqlconf.F90 diag.F90  diagcfac.F90   \
	diagdens.F90  diagdenz.F90  diagentr.F90  diagescl.F90  diaggnde.F90  \
	diaggnde2.F90  diagimpd.F90  diagscal.F90  diagwrng.F90  diagxswt.F90   \
	diagxswx.F90  dsk_gr.F90 dskout.F90  efield.F90  eflditer.F90 eqalloc.F90   \
	eqcoord.F90  eqelpse.F90  eqflxavg.F90  eqfn.F90  \
	eqfndpsi.F90  eqfninv.F90  eqfpsi.F90   \
	eqinitl.F90  eqjac.F90  eqonovrp.F90  \
	eqorbit.F90  eqrhs.F90   \
	eqtopeol.F90  eqvolpsi.F90  eqwrng.F90  \
	esefld.F90  exlin.F90  exsweep.F90  exsweept.F90   \
	exsweepx.F90  finit.F90  firstdrv.F90  fle.F90 flxfn.F90  \
	freya.F  freyasou.F  \
	frinitl.F  frinitz.F  frnbdep2.F  frnfreya.F  \
	frset.F frsmooth.F  frsplft.F   \
	frstup.F  frsubs.F  frsuppor.F  frwrong.F  \
	ilut.F90  impchk.F90  impnorm.F90  \
	lookup.F90  losscone.F90  lossegy.F90  lossorbm.F90  \
	losstor.F90  micgetr.F90   \
	micgmbnd.F90  micgnbnd.F90  micxinil.F90  micxinim.F90  micxinit.F90  \
	micxiniz.F90  \
	netcdfrf.F  netcdfrw2.F \
	ntdstore.F90   \
	ntloop.F90  pack21.F90 pltcycl.F90 \
	pltdnz.F90  pltelec.F90  pltendn.F90   \
	pltfvsv.F90  pltinit.F90  pltlosc.F90  \
	pltpower.F90  \
	pltrun.F90  pltvec.F90  pltvectr.F90  \
	pltvflux.F90  profaxis.F90  profiles.F90 prppr.F90   \
	prpprctr.F90  psif.F90 r8lsode.F  \
	rdc_multi.F90 rdc_bplt.F90 restcon.F90 resthks.F90  \
	restvty.F90  rf.F90  sigalloc.F90  siggy.F90  sigmax.F90  sigsetup.F90  \
	sigv5d.F90 sigfn.F90 sigie.F90  sigmaxwl.F90   sigv.F90  \
	soucrit.F90 sounorm.F90  soup.F90  soup0.F90  \
	sourc0.F90  sourcee.F90  sourcef.F90  sourceko.F90 sourcpwr.F90   \
	synchrad.F90  tdbootst.F90 tdboothi.F90  tdchief.F90  tddiag.F90  \
	tddsig.F90  tdfinterp.F90  tdinitl.F90  tdinlegw.F90  \
	tdinterp.F90  tdnflxs.F90  tdnpa.F90  tdnpadiag.F90  tdnpa0.F90  \
	tdnpacxcs.F90  tdnpalam.F90  tdnpabscs.F90  tdoutput.F90   \
	tdplteq.F90  tdpltjop.F90  tdpltmne.F90  tdpro.F90  \
	tdreadf.F90  tdrmshst.F90  tdsetnpa.F90  tdsetsxr.F90   \
	tdstin.F90  tdsxr.F90  tdsxr0.F90  tdsxray.F90  tdsxrplt.F90  \
	tdtloop.F90  tdtoarad.F90   tdtranspn.F90  \
	tdtoaray.F90  tdtraloc.F90  tdtransp.F90  tdtravct.F90 \
	tdtrchk.F90  tdtrchkd.F90  tdtrcon.F90   \
	tdtrdfus.F90  tdtrfcop.F90  tdtrflg.F90  tdtrflx.F90  \
	tdtrmuy.F90  tdtrrsou.F90  tdtrrtov.F90   \
	tdtrrtov2.F90  tdtrsavf.F90  tdtrsym.F90  tdtrvint.F90  \
	tdtrvsou.F90  tdtrvtor.F90  tdtrvtor2.F90   \
	tdtrvtor3.F90  tdtrwtl.F90  tdtry.F90  tdtscinp.F90  tdtscout.F90  \
	tdwrng.F90  tdwritef.F90  tdxin13d.F90   \
	tdxin23d.F90  tdxin33d.F90  tdxinitl.F90  urfalloc.F90  \
	urfavg.F90  urfb0.F90  urfbes.F90   \
	urfbplt.F90  urfchief.F90  urfdamp0.F90  urfdamp1.F90  \
	urfdamp2.F90  urfdampa.F90 urfdout.F90  urfedge.F90   \
	urffflx.F90  urfinitl.F90  urfmidv.F90  \
	urfpack.F90  urfpackm.F90  urfrays.F90  urfread.F90   \
	urfread_.F90  urfsetup.F90  urfwrite.F90  urfwrite_.F90  \
	urfwrong.F90  urfwr0.F90  urfwr0c.F90 \
	vlf.F90 vlfalloc.F90 vlfbplt.F90 vlfsetup.F90 vlh.F90  vlhbplt.F90 vlhd.F90  \
	wpalloc.F90  wparsou.F90  wpavg.F90  wpbdry.F90  \
	wpcheck.F90  wpchgdy.F90 wpcthta.F90  wpelecf.F90   \
	wpinitl.F90  wploweq.F90  wpsavf.F90  wptrafx.F90  wptramu.F90  \
	wptrmuy.F90  wpvptb.F90  wpwrng.F90  wpmshchk.F90 \
	zblock.F90  zcunix.F90  zfreya.F   znonsym.F90
