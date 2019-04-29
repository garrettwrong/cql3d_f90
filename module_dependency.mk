#------ Module generated dependencies, to be included in makefiles...
a_cqlp.o: abchief.o impavnc0.o 
abchief.o: comm.o param.o tdchief.o 
achief1.o: achiefn.o ainalloc.o aindflt1.o aindflt.o \
	aindfpa.o aingeom.o ainitial.o ainpla.o ainsetva.o ainspec.o \
	ainvnorm.o cfpgamma.o coefmidt.o coefmidv.o coefstup.o coefwti.o \
	coefwtj.o comm.o diagimpd.o diagscal.o eqindflt.o eqinitl.o \
	micxinit.o micxiniz.o ntloop.o param.o pltmain.o profiles.o \
	r8subs.o tdnflxs.o wploweq.o 
achiefn.o: cfpcoefc.o cfpgamma.o comm.o diag.o diaggnde.o dskout.o \
	efield.o exsweep.o finit.o impavnc0.o ntdstore.o ntloop.o param.o \
	pltmain.o pltrun.o restvty.o sourcee.o tdoutput.o 
aclear.o: bcast.o comm.o param.o 
advnce.o: bsl.o bsu.o comm.o r8subs.o 
ainalloc.o: bcast.o comm.o param.o 
aindflt.o: comm.o param.o 
aindflt1.o: comm.o param.o 
aindfpa.o: comm.o param.o 
aingeom.o: comm.o diagwrng.o eqalloc.o eqcoord.o eqelpse.o \
	equilib.o flxfn.o param.o zcunix.o 
ainitial.o: bavgmax.o coefefld.o comm.o diag.o diagentr.o \
	diaggnde.o efield.o finit.o losscone.o lossegy.o micxinim.o \
	ntdstore.o param.o pltmain.o restvty.o sourcee.o synchrad.o \
	tdoutput.o vlf.o vlh.o wpvptb.o 
ainpla.o: cfpmodbe.o comm.o param.o 
ainplt.o: comm.o param.o 
ainpltpa.o: comm.o param.o 
ainsetpa.o: comm.o param.o 
ainsetva.o: comm.o diagwrng.o param.o wpwrng.o 
ainspec.o: comm.o param.o 
ainvnorm.o: comm.o diagwrng.o param.o 
ampfar.o: bcast.o comm.o param.o r8subs.o 
bavdens.o: bcast.o comm.o exlin.o param.o psif.o r8subs.o 
bavgmax.o: bcast.o comm.o param.o psif.o 
baviorbt.o: bcast.o comm.o lookup.o micgnbnd.o param.o psif.o \
	r8subs.o zcunix.o 
bsl.o: bcast.o comm.o lookup.o param.o 
bsu.o: bcast.o comm.o lookup.o param.o 
cfpcoefc.o: cfpcoefn.o cfpcoefr.o comm.o param.o 
cfpcoefn.o: bavdens.o bcast.o cfpleg.o cfpmodbe.o cfpsymt.o comm.o \
	diagescl.o param.o r8subs.o 
cfpcoefr.o: bavdens.o bcast.o cfpleg.o cfpmodbe.o cfpsymt.o comm.o \
	diagescl.o param.o r8subs.o 
cfpgamma.o: comm.o param.o 
cfpleg.o: bcast.o comm.o param.o r8subs.o 
cfpsymt.o: comm.o param.o 
coefefad.o: comm.o param.o 
coefefld.o: bcast.o comm.o param.o 
coefegad.o: comm.o param.o 
coeffpad.o: comm.o param.o 
coefload.o: comm.o losstor.o param.o 
coefmidt.o: bcast.o comm.o param.o 
coefmidv.o: bcast.o comm.o param.o r8subs.o 
coefrfad.o: bcast.o comm.o param.o 
coefstup.o: bcast.o coefefad.o coefegad.o coeffpad.o coefload.o \
	coefrfad.o coefsyad.o comm.o param.o r8subs.o 
coefsyad.o: comm.o param.o 
coefwti.o: advnce.o bcast.o coefmidt.o comm.o param.o 
coefwtj.o: advnce.o bcast.o coefmidv.o comm.o param.o 
comm.o: param.o 
diag.o: comm.o param.o 
diagcfac.o: bcast.o zcunix.o 
diagdens.o: bcast.o comm.o param.o 
diagdenz.o: bcast.o cfpleg.o comm.o param.o r8subs.o 
diagentr.o: advnce.o bcast.o coefefad.o coefegad.o coefload.o \
	coefmidv.o coefrfad.o coefsyad.o comm.o diagwrng.o param.o \
	r8subs.o 
diagescl.o: bcast.o comm.o param.o r8subs.o 
diaggnde.o: bcast.o comm.o diagcfac.o diagdenz.o diaggnde2.o \
	diagwrng.o param.o r8subs.o soucrit.o 
diaggnde2.o: bcast.o cfpleg.o comm.o diagcfac.o diagdenz.o \
	diagwrng.o param.o r8subs.o 
diagimpd.o: advnce.o aminmx.o bcast.o comm.o diagdens.o diagentr.o \
	param.o r8subs.o 
diagscal.o: bcast.o comm.o param.o r8subs.o 
diagwrng.o: comm.o param.o 
diagxswt.o: comm.o diagdens.o diagentr.o param.o r8subs.o 
diagxswx.o: advnce.o bcast.o coefmidt.o coefmidv.o coefstup.o \
	comm.o diagdens.o diagentr.o param.o r8subs.o 
dsk_gr.o: comm.o param.o 
dskout.o: comm.o param.o 
efield.o: cfpgamma.o comm.o param.o restcon.o resthks.o 
eflditer.o: bcast.o coefefad.o comm.o param.o r8subs.o restcon.o \
	resthks.o soucrit.o 
eqalloc.o: bcast.o comm.o param.o 
eqcoord.o: aminmx.o comm.o eqfndpsi.o eqrhopsi.o param.o 
eqelpse.o: comm.o eqfn.o param.o 
eqflxavg.o: comm.o eqorbit.o param.o 
eqfn.o: comm.o eqwrng.o param.o 
eqfndpsi.o: comm.o eqflxavg.o eqfpsi.o eqonovrp.o eqorbit.o \
	eqvolpsi.o eqwrng.o param.o 
eqfninv.o: comm.o eqwrng.o param.o 
eqfpsi.o: comm.o eqwrng.o param.o zcunix.o 
eqindflt.o: comm.o param.o 
eqinitl.o: comm.o param.o 
eqjac.o: comm.o param.o zcunix.o 
eqonovrp.o: comm.o eqflxavg.o eqorbit.o param.o 
eqorbit.o: aminmx.o bcast.o comm.o eqfpsi.o eqjac.o eqrhs.o \
	eqwrng.o exlin.o param.o zcunix.o 
eqrhopsi.o: aminmx.o comm.o eqflxavg.o eqfpsi.o eqonovrp.o \
	eqorbit.o eqvolpsi.o eqwrng.o exlin.o param.o zcunix.o 
eqrhs.o: comm.o param.o zcunix.o 
eqtopeol.o: comm.o param.o r8subs.o zcunix.o 
equilib.o: comm.o eqtopeol.o eqwrng.o param.o r8subs.o zcunix.o 
eqvolpsi.o: comm.o eqorbit.o param.o 
eqwrng.o: param.o 
esefld.o: param.o 
exsweep.o: advnce.o coefmidt.o coefmidv.o coefstup.o coefwti.o \
	coefwtj.o comm.o diagxswt.o diagxswx.o exsweept.o exsweepx.o \
	param.o r8subs.o 
exsweept.o: advnce.o comm.o diagwrng.o param.o 
exsweepx.o: advnce.o comm.o diagentr.o diagwrng.o param.o 
finit.o: bcast.o comm.o lossorbm.o param.o tdreadf.o 
freyasou.o: bcast.o comm.o eqfpsi.o param.o r8subs.o sourcpwr.o tdnflxs.o \
	tdtoaray.o urfb0.o zcunix.o
frinitz.o: comm.o param.o
frnfreya.o: comm.o frplteq.o param.o
frsmooth.o: comm.o param.o
frstup.o: bcast.o comm.o param.o
frwrong.o: comm.o param.o
fle.o: bcast.o comm.o diagwrng.o micgetr.o param.o r8subs.o 
flxfn.o: comm.o diagwrng.o param.o 
frplteq.o: comm.o param.o r8subs.o tdnflxs.o 
hpalloc0.o: param.o 
ilut.o: r8subs.o 
impavnc.o: advnce.o bcast.o coefmidt.o coefmidv.o coefstup.o \
	coefwti.o coefwtj.o comm.o diagimpd.o esefld.o impchk.o impnorm.o \
	param.o r8subs.o tdtrvsou.o 
impavnc0.o: advnce.o bcast.o bsl.o bsu.o coefmidt.o coefmidv.o \
	coefstup.o coefwti.o coefwtj.o comm.o esefld.o ilut.o impchk.o \
	impnorm.o param.o r8subs.o tdtranspn.o tdtrvsou.o 
impchk.o: advnce.o bcast.o comm.o diagentr.o diagwrng.o param.o 
lookup.o: urfb0.o 
losscone.o: bcast.o comm.o lossorbm.o param.o r8subs.o zcunix.o 
lossegy.o: bcast.o comm.o param.o 
lossorbm.o: comm.o param.o r8subs.o 
losstor.o: comm.o param.o 
micgmbnd.o: comm.o param.o psif.o 
micgnbnd.o: comm.o micgmbnd.o param.o 
micxinil.o: bcast.o comm.o param.o psif.o r8subs.o tdinlegw.o \
	wpwrng.o 
micxinim.o: comm.o param.o 
micxinit.o: comm.o diagwrng.o micgetr.o param.o tdtry.o zcunix.o 
micxiniz.o: comm.o diagwrng.o micgetr.o param.o psif.o tdxin13d.o zcunix.o 
mpilib.o: comm.o param.o r8subs.o 
netcdfrf.o: bcast.o comm.o pack21.o param.o
netcdfrw2.o: advnce.o bcast.o comm.o diagentr.o param.o prppr.o r8subs.o zcunix.o
ntdstore.o: comm.o param.o restvty.o 
ntloop.o: comm.o param.o 
pack21.o: bcast.o r8subs.o 
pltdf.o: comm.o lookup.o param.o r8subs.o 
pltdnz.o: aminmx.o comm.o param.o pltelec.o 
pltelec.o: aminmx.o comm.o param.o 
pltendn.o: aminmx.o comm.o param.o r8subs.o 
pltfvsv.o: aminmx.o bcast.o comm.o param.o pltdf.o r8subs.o 
pltinit.o: comm.o param.o 
pltlosc.o: comm.o param.o pltdf.o 
pltmain.o:  advnce.o aminmx.o bcast.o coefefad.o coefegad.o coeffpad.o \
	coefload.o coefmidt.o coefmidv.o coefrfad.o coefstup.o coefsyad.o \
	comm.o fle.o diagentr.o param.o pltdf.o pltdnz.o pltendn.o \
	pltfvsv.o pltlosc.o pltpower.o  \
	pltvec.o pltvflux.o r8subs.o
pltpower.o: aminmx.o comm.o param.o 
pltrun.o: aminmx.o comm.o param.o r8subs.o tdnflxs.o 
pltvec.o: advnce.o aminmx.o bcast.o coefefad.o coeffpad.o \
	coefmidt.o coefmidv.o coefrfad.o coefstup.o comm.o diagentr.o \
	param.o pltvectr.o prppr.o r8subs.o 
pltvectr.o: aminmx.o r8subs.o 
pltvflux.o: aminmx.o comm.o param.o r8subs.o 
profaxis.o: comm.o param.o 
profiles.o: cfpmodbe.o comm.o param.o profaxis.o tdinterp.o \
	tdxin13d.o 
prppr.o: bcast.o comm.o diagwrng.o param.o r8subs.o 
prpprctr.o: aminmx.o comm.o param.o pltdf.o pltmain.o 
psif.o: comm.o diagwrng.o param.o zcunix.o 
r8lsode.o: r8subs.o
r8subs.o: comm.o 
rdc_bplt.o: bcast.o comm.o param.o pltdf.o 
rdc_multi.o: bcast.o comm.o netcdfrf.o param.o r8subs.o rdc_bplt.o \
	tdnflxs.o zcunix.o 
restcon.o: comm.o param.o 
resthks.o: comm.o param.o 
restvty.o: comm.o param.o 
rf.o: comm.o param.o 
sigalloc.o: bcast.o comm.o param.o 
siggy.o: comm.o param.o 
sigmax.o: bcast.o comm.o param.o 
sigmaxwl.o: bcast.o comm.o param.o 
sigsetup.o: bcast.o comm.o param.o r8subs.o sigfn.o siggy.o 
sigv.o: bcast.o comm.o param.o sigalloc.o sigsetup.o sigv5d.o \
	tdnflxs.o 
sigv5d.o: bcast.o cfpleg.o comm.o param.o r8subs.o sigmax.o 
soucrit.o: cfpgamma.o comm.o param.o r8subs.o 
sounorm.o: bcast.o comm.o param.o soup.o 
soup.o: comm.o param.o psif.o r8subs.o soup0.o 
soup0.o: comm.o param.o 
sourc0.o: comm.o param.o 
sourcee.o: bcast.o comm.o param.o sounorm.o sourc0.o sourcef.o \
	sourceko.o sourcpwr.o 
sourcef.o: bcast.o comm.o param.o r8subs.o soup.o 
sourceko.o: bcast.o comm.o fle.o lookup.o param.o r8subs.o \
	soucrit.o tdoutput.o 
sourcpwr.o: bcast.o comm.o param.o 
synchrad.o: bcast.o comm.o param.o 
tdboothi.o: bcast.o comm.o param.o zcunix.o 
tdbootst.o: bcast.o comm.o param.o zcunix.o 
tdchief.o: achief1.o achiefn.o aclear.o aindfpa.o ainplt.o \
	ainpltpa.o ainsetpa.o ampfar.o cfpgamma.o coefmidt.o coefmidv.o \
	coefstup.o coefwti.o coefwtj.o comm.o diag.o diagentr.o diaggnde.o \
	diagimpd.o diagscal.o dsk_gr.o dskout.o eflditer.o esefld.o \
	netcdfrf.o ntdstore.o param.o pltendn.o pltinit.o pltmain.o \
	pltrun.o profiles.o r8subs.o restvty.o sigv.o tddiag.o tdinitl.o \
	tdnflxs.o tdnpadiag.o tdoutput.o tdplteq.o tdpltmne.o tdsxray.o \
	tdtloop.o tdtoaray.o tdtransp.o tdtrcon.o tdtrdfus.o tdtrfcop.o \
	tdtrrsou.o tdtrsavf.o urfchief.o urfwrite.o wparsou.o wpavg.o \
	wpelecf.o wpsavf.o wptrafx.o wptramu.o 
tddiag.o: comm.o param.o tdboothi.o tdbootst.o 
tdeqdsk.o: comm.o equilib.o firstdrv.o param.o r8subs.o zcunix.o 
tdfinterp.o: comm.o lookup.o param.o 
tdinitl.o: ainalloc.o aindflt1.o aindflt.o aindfpa.o aingeom.o \
	ainitial.o ainpla.o ainsetva.o ainspec.o ainvnorm.o ampfar.o \
	baviorbt.o comm.o diagwrng.o eqindflt.o eqinitl.o micxinil.o \
	micxinit.o micxiniz.o param.o profiles.o rdc_multi.o sigv.o \
	tddiag.o tdeqdsk.o tdnflxs.o tdnpadiag.o tdoutput.o tdplteq.o \
	tdpltmne.o tdreadf.o tdrmshst.o tdstin.o tdsxray.o tdtoarad.o \
	tdtoaray.o tdtraloc.o tdtrdfus.o tdtrflg.o tdtrmuy.o tdtrvint.o \
	tdtrwtl.o tdwrng.o tdxinitl.o urfindfl.o urfinitl.o urfsetup.o \
	wpalloc.o wpchgdy.o wpinitl.o wploweq.o wpmshchk.o wptrmuy.o 
tdinlegw.o: bcast.o comm.o param.o 
tdinterp.o: zcunix.o 
tdnflxs.o: comm.o diagwrng.o param.o 
tdnpa.o: bcast.o comm.o param.o tdfinterp.o tdnpalam.o 
tdnpa0.o: bcast.o comm.o eqfpsi.o param.o tdnpa.o tdnpalam.o \
	tdsetnpa.o tdsxrplt.o zcunix.o 
tdnpadiag.o: comm.o param.o tdnpa0.o 
tdnpalam.o: comm.o param.o tdnpabscs.o 
tdoutput.o: bcast.o cfpgamma.o comm.o param.o restcon.o resthks.o \
	tdnflxs.o tdtrflx.o 
tdplteq.o: comm.o frplteq.o param.o 
tdpltjop.o: aminmx.o comm.o param.o 
tdpltmne.o: aminmx.o comm.o param.o r8subs.o tdpltjop.o 
tdreadf.o: bcast.o comm.o param.o tdnflxs.o zcunix.o 
tdrmshst.o: comm.o eqorbit.o eqvolpsi.o param.o zcunix.o 
tdsetnpa.o: bcast.o comm.o param.o zcunix.o 
tdsetsxr.o: bcast.o comm.o param.o tddsig.o 
tdstin.o: comm.o param.o 
tdsxr.o: bcast.o cfpleg.o comm.o param.o r8subs.o tdnflxs.o \
	tdwrng.o zcunix.o 
tdsxr0.o: bcast.o comm.o eqfpsi.o param.o tdsetsxr.o tdsxr.o \
	tdsxrplt.o zcunix.o 
tdsxray.o: comm.o param.o tdsxr0.o 
tdsxrplt.o: aminmx.o comm.o param.o r8subs.o 
tdtloop.o: aindfpa.o comm.o param.o tdeqdsk.o tdtscout.o tdwritef.o 
tdtoarad.o: comm.o param.o 
tdtoaray.o: comm.o param.o 
tdtraloc.o: bcast.o comm.o param.o 
tdtransp.o: comm.o param.o r8subs.o tdtravct.o tdtrchk.o \
	tdtrrtov2.o tdtrrtov.o tdtrsym.o tdtrvtor2.o tdtrvtor.o tdtrwtl.o 
tdtranspn.o: comm.o param.o r8subs.o tdtravct.o tdtrvtor.o \
	tdtrwtl.o 
tdtravct.o: bcast.o comm.o param.o r8subs.o tdtrrtov2.o tdtrrtov.o \
	tdtrsym.o tdtrvtor2.o tdtrvtor.o tdtrwtl.o 
tdtrchk.o: bcast.o comm.o diagwrng.o param.o r8subs.o 
tdtrchkd.o: bcast.o comm.o param.o 
tdtrcon.o: comm.o param.o tdtrflx.o 
tdtrdfus.o: bcast.o comm.o param.o tdnflxs.o 
tdtrfcop.o: comm.o param.o r8subs.o 
tdtrflg.o: bcast.o comm.o param.o 
tdtrflx.o: bcast.o comm.o param.o r8subs.o 
tdtrmuy.o: comm.o micgetr.o param.o 
tdtrrsou.o: comm.o param.o r8subs.o tdtravct.o tdtrrtov2.o \
	tdtrvtor2.o 
tdtrrtov.o: comm.o param.o tdtrchkd.o 
tdtrrtov2.o: comm.o param.o r8subs.o 
tdtrsavf.o: comm.o param.o 
tdtrsym.o: comm.o param.o 
tdtrvint.o: bcast.o comm.o param.o 
tdtrvsou.o: advnce.o comm.o diagentr.o param.o tdtrvtor3.o 
tdtrvtor.o: comm.o param.o r8subs.o tdtrchkd.o 
tdtrvtor2.o: comm.o param.o r8subs.o 
tdtrvtor3.o: comm.o param.o r8subs.o 
tdtrwtl.o: comm.o param.o 
tdtry.o: comm.o param.o tdwrng.o 
tdtscinp.o: comm.o param.o tdeqdsk.o 
tdtscout.o: comm.o param.o tdinterp.o 
tdwritef.o: comm.o param.o 
tdwrng.o: param.o 
tdxin13d.o: param.o profaxis.o 
tdxin23d.o: param.o profaxis.o 
tdxin33d.o: param.o profaxis.o 
tdxinitl.o: bcast.o comm.o micgetr.o param.o profaxis.o tdinterp.o \
	tdpro.o tdtscinp.o tdwrng.o tdxin13d.o tdxin23d.o tdxin33d.o 
urfalloc.o: bcast.o comm.o param.o 
urfavg.o: comm.o param.o 
urfb0.o: bcast.o comm.o param.o tdnflxs.o 
urfbes.o: aminmx.o bcast.o comm.o param.o zcunix.o 
urfbplt.o: bcast.o comm.o param.o pltdf.o 
urfchief.o: comm.o param.o tdnflxs.o urfavg.o urfb0.o urfbes.o \
	urfbplt.o urfdamp0.o urffflx.o urfpack.o urfrays.o urfread.o \
	urfwrite.o 
urfdamp0.o: bcast.o comm.o param.o r8subs.o urfdamp1.o urfdamp2.o \
	urfdampa.o 
urfdamp1.o: comm.o param.o tdnflxs.o 
urfdamp2.o: bcast.o comm.o param.o tdnflxs.o urfmidv.o 
urfdampa.o: comm.o param.o 
urfdout.o: comm.o param.o 
urfedge.o: comm.o param.o 
urffflx.o: bcast.o comm.o param.o r8subs.o zcunix.o 
urfindfl.o: comm.o param.o 
urfinitl.o: comm.o param.o 
urfmidv.o: bcast.o comm.o param.o r8subs.o 
urfpack.o: comm.o param.o r8subs.o tdnflxs.o urfedge.o urfwrong.o 
urfrays.o: comm.o param.o 
urfread.o: comm.o netcdfrf.o param.o r8subs.o urfread_.o 
urfread_.o: comm.o param.o urfwrong.o 
urfsetup.o: bcast.o comm.o netcdfrf.o param.o urfalloc.o urfread_.o 
urfwrite.o: comm.o param.o urfwrite_.o 
urfwrite_.o: comm.o param.o urfwr0.o urfwr0c.o 
urfwrong.o: param.o 
vlf.o: bcast.o comm.o param.o r8subs.o urfedge.o urfwrong.o \
	vlfbplt.o vlfsetup.o zcunix.o 
vlfalloc.o: bcast.o comm.o param.o 
vlfbplt.o: bcast.o comm.o param.o pltdf.o r8subs.o urfwrong.o 
vlfsetup.o: bcast.o comm.o param.o vlfalloc.o zcunix.o 
vlh.o: bcast.o comm.o param.o vlhbplt.o vlhd.o 
vlhbplt.o: bcast.o comm.o param.o pltdf.o 
vlhd.o: comm.o param.o 
wpalloc.o: bcast.o comm.o param.o 
wparsou.o: bcast.o comm.o param.o 
wpavg.o: comm.o param.o 
wpbdry.o: comm.o param.o r8subs.o 
wpcheck.o: bcast.o comm.o param.o 
wpchgdy.o: comm.o param.o 
wpcthta.o: comm.o param.o 
wpelecf.o: advnce.o comm.o param.o r8subs.o 
wpinitl.o: bcast.o comm.o param.o r8subs.o wpwrng.o 
wploweq.o: comm.o param.o 
wpmshchk.o: comm.o param.o 
wpsavf.o: comm.o param.o r8subs.o 
wptrafx.o: bcast.o comm.o param.o r8subs.o wpbdry.o wpcheck.o \
	wpwrng.o znonsym.o 
wptramu.o: bcast.o comm.o param.o r8subs.o wpbdry.o wpcheck.o \
	znonsym.o 
wptrmuy.o: comm.o param.o 
wpvptb.o: comm.o param.o 
wpwrng.o: comm.o param.o 
