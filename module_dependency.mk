#------ Module generated dependencies, to be included in makefiles...
a_cqlp.o: abchief.o impavnc0.o 
abchief.o: cqlcomm.o param.o tdchief.o 
achief1.o: achiefn.o ainalloc.o aindflt1.o aindflt.o \
	aindfpa.o aingeom.o ainitial.o ainpla.o ainsetva.o ainspec.o \
	ainvnorm.o cfpgamma.o coefmidt.o coefmidv.o coefstup.o coefwti.o \
	coefwtj.o cqlcomm.o diagimpd.o diagscal.o eqindflt.o eqinitl.o \
	micxinit.o micxiniz.o ntloop.o param.o pltmain.o profiles.o \
	r8subs.o tdnflxs.o wploweq.o 
achiefn.o: cfpcoefc.o cfpgamma.o cqlcomm.o diag.o diaggnde.o dskout.o \
	efield.o exsweep.o finit.o impavnc0.o ntdstore.o ntloop.o param.o \
	pltmain.o pltrun.o restvty.o sourcee.o tdoutput.o 
aclear.o: bcast.o cqlcomm.o param.o 
advnce.o: bsl.o bsu.o cqlcomm.o r8subs.o 
ainalloc.o: bcast.o cqlcomm.o param.o 
aindflt.o: cqlcomm.o param.o 
aindflt1.o: cqlcomm.o param.o 
aindfpa.o: cqlcomm.o param.o 
aingeom.o: cqlcomm.o diagwrng.o eqalloc.o eqcoord.o eqelpse.o \
	equilib.o flxfn.o param.o zcunix.o 
ainitial.o: bavgmax.o coefefld.o cqlcomm.o diag.o diagentr.o \
	diaggnde.o efield.o finit.o losscone.o lossegy.o micxinim.o \
	ntdstore.o param.o pltmain.o restvty.o sourcee.o synchrad.o \
	tdoutput.o vlf.o vlh.o wpvptb.o 
ainpla.o: cfpmodbe.o cqlcomm.o param.o 
ainplt.o: cqlcomm.o param.o 
ainpltpa.o: cqlcomm.o param.o 
ainsetpa.o: cqlcomm.o param.o 
ainsetva.o: cqlcomm.o diagwrng.o param.o wpwrng.o 
ainspec.o: cqlcomm.o param.o 
ainvnorm.o: cqlcomm.o diagwrng.o param.o 
ampfar.o: bcast.o cqlcomm.o param.o r8subs.o 
bavdens.o: bcast.o cqlcomm.o exlin.o param.o psif.o r8subs.o 
bavgmax.o: bcast.o cqlcomm.o param.o psif.o 
baviorbt.o: bcast.o cqlcomm.o lookup.o micgnbnd.o param.o psif.o \
	r8subs.o zcunix.o 
bsl.o: bcast.o cqlcomm.o lookup.o param.o 
bsu.o: bcast.o cqlcomm.o lookup.o param.o 
cfpcoefc.o: cfpcoefn.o cfpcoefr.o cqlcomm.o param.o 
cfpcoefn.o: bavdens.o bcast.o cfpleg.o cfpmodbe.o cfpsymt.o cqlcomm.o \
	diagescl.o param.o r8subs.o 
cfpcoefr.o: bavdens.o bcast.o cfpleg.o cfpmodbe.o cfpsymt.o cqlcomm.o \
	diagescl.o param.o r8subs.o 
cfpgamma.o: cqlcomm.o param.o 
cfpleg.o: bcast.o cqlcomm.o param.o r8subs.o 
cfpsymt.o: cqlcomm.o param.o 
coefefad.o: cqlcomm.o param.o 
coefefld.o: bcast.o cqlcomm.o param.o 
coefegad.o: cqlcomm.o param.o 
coeffpad.o: cqlcomm.o param.o 
coefload.o: cqlcomm.o losstor.o param.o 
coefmidt.o: bcast.o cqlcomm.o param.o 
coefmidv.o: bcast.o cqlcomm.o param.o r8subs.o 
coefrfad.o: bcast.o cqlcomm.o param.o 
coefstup.o: bcast.o coefefad.o coefegad.o coeffpad.o coefload.o \
	coefrfad.o coefsyad.o cqlcomm.o param.o r8subs.o 
coefsyad.o: cqlcomm.o param.o 
coefwti.o: advnce.o bcast.o coefmidt.o cqlcomm.o param.o 
coefwtj.o: advnce.o bcast.o coefmidv.o cqlcomm.o param.o 
cqlcomm.o: param.o 
diag.o: cqlcomm.o param.o 
diagcfac.o: bcast.o zcunix.o 
diagdens.o: bcast.o cqlcomm.o param.o 
diagdenz.o: bcast.o cfpleg.o cqlcomm.o param.o r8subs.o 
diagentr.o: advnce.o bcast.o coefefad.o coefegad.o coefload.o \
	coefmidv.o coefrfad.o coefsyad.o cqlcomm.o diagwrng.o param.o \
	r8subs.o 
diagescl.o: bcast.o cqlcomm.o param.o r8subs.o 
diaggnde.o: bcast.o cqlcomm.o diagcfac.o diagdenz.o diaggnde2.o \
	diagwrng.o param.o r8subs.o soucrit.o 
diaggnde2.o: bcast.o cfpleg.o cqlcomm.o diagcfac.o diagdenz.o \
	diagwrng.o param.o r8subs.o 
diagimpd.o: advnce.o aminmx.o bcast.o cqlcomm.o diagdens.o diagentr.o \
	param.o r8subs.o 
diagscal.o: bcast.o cqlcomm.o param.o r8subs.o 
diagwrng.o: cqlcomm.o param.o 
diagxswt.o: cqlcomm.o diagdens.o diagentr.o param.o r8subs.o 
diagxswx.o: advnce.o bcast.o coefmidt.o coefmidv.o coefstup.o \
	cqlcomm.o diagdens.o diagentr.o param.o r8subs.o 
dsk_gr.o: cqlcomm.o param.o 
dskout.o: cqlcomm.o param.o 
efield.o: cfpgamma.o cqlcomm.o param.o restcon.o resthks.o 
eflditer.o: bcast.o coefefad.o cqlcomm.o param.o r8subs.o restcon.o \
	resthks.o soucrit.o 
eqalloc.o: bcast.o cqlcomm.o param.o 
eqcoord.o: aminmx.o cqlcomm.o eqfndpsi.o eqrhopsi.o param.o 
eqelpse.o: cqlcomm.o eqfn.o param.o 
eqflxavg.o: cqlcomm.o eqorbit.o param.o 
eqfn.o: cqlcomm.o eqwrng.o param.o 
eqfndpsi.o: cqlcomm.o eqflxavg.o eqfpsi.o eqonovrp.o eqorbit.o \
	eqvolpsi.o eqwrng.o param.o 
eqfninv.o: cqlcomm.o eqwrng.o param.o 
eqfpsi.o: cqlcomm.o eqwrng.o param.o zcunix.o 
eqindflt.o: cqlcomm.o param.o 
eqinitl.o: cqlcomm.o param.o 
eqjac.o: cqlcomm.o param.o zcunix.o 
eqonovrp.o: cqlcomm.o eqflxavg.o eqorbit.o param.o 
eqorbit.o: aminmx.o bcast.o cqlcomm.o eqfpsi.o eqjac.o eqrhs.o \
	eqwrng.o exlin.o param.o zcunix.o 
eqrhopsi.o: aminmx.o cqlcomm.o eqflxavg.o eqfpsi.o eqonovrp.o \
	eqorbit.o eqvolpsi.o eqwrng.o exlin.o param.o zcunix.o 
eqrhs.o: cqlcomm.o param.o zcunix.o 
eqtopeol.o: cqlcomm.o param.o r8subs.o zcunix.o 
equilib.o: cqlcomm.o eqtopeol.o eqwrng.o param.o r8subs.o zcunix.o 
eqvolpsi.o: cqlcomm.o eqorbit.o param.o 
eqwrng.o: param.o 
esefld.o: param.o 
exsweep.o: advnce.o coefmidt.o coefmidv.o coefstup.o coefwti.o \
	coefwtj.o cqlcomm.o diagxswt.o diagxswx.o exsweept.o exsweepx.o \
	param.o r8subs.o 
exsweept.o: advnce.o cqlcomm.o diagwrng.o param.o 
exsweepx.o: advnce.o cqlcomm.o diagentr.o diagwrng.o param.o 
finit.o: bcast.o cqlcomm.o lossorbm.o param.o tdreadf.o 
freyasou.o: bcast.o cqlcomm.o eqfpsi.o param.o r8subs.o sourcpwr.o tdnflxs.o \
	tdtoaray.o urfb0.o zcunix.o
frinitz.o: cqlcomm.o param.o
frnfreya.o: cqlcomm.o frplteq.o param.o
frsmooth.o: cqlcomm.o param.o
frstup.o: bcast.o cqlcomm.o param.o
frwrong.o: cqlcomm.o param.o
fle.o: bcast.o cqlcomm.o diagwrng.o micgetr.o param.o r8subs.o 
flxfn.o: cqlcomm.o diagwrng.o param.o 
frplteq.o: cqlcomm.o param.o r8subs.o tdnflxs.o 
hpalloc0.o: param.o 
ilut.o: r8subs.o 
impavnc.o: advnce.o bcast.o coefmidt.o coefmidv.o coefstup.o \
	coefwti.o coefwtj.o cqlcomm.o diagimpd.o esefld.o impchk.o impnorm.o \
	param.o r8subs.o tdtrvsou.o 
impavnc0.o: advnce.o bcast.o bsl.o bsu.o coefmidt.o coefmidv.o \
	coefstup.o coefwti.o coefwtj.o cqlcomm.o esefld.o ilut.o impchk.o \
	impnorm.o param.o r8subs.o tdtranspn.o tdtrvsou.o 
impchk.o: advnce.o bcast.o cqlcomm.o diagentr.o diagwrng.o param.o 
lookup.o: urfb0.o 
losscone.o: bcast.o cqlcomm.o lossorbm.o param.o r8subs.o zcunix.o 
lossegy.o: bcast.o cqlcomm.o param.o 
lossorbm.o: cqlcomm.o param.o r8subs.o 
losstor.o: cqlcomm.o param.o 
micgmbnd.o: cqlcomm.o param.o psif.o 
micgnbnd.o: cqlcomm.o micgmbnd.o param.o 
micxinil.o: bcast.o cqlcomm.o param.o psif.o r8subs.o tdinlegw.o \
	wpwrng.o 
micxinim.o: cqlcomm.o param.o 
micxinit.o: cqlcomm.o diagwrng.o micgetr.o param.o tdtry.o zcunix.o 
micxiniz.o: cqlcomm.o diagwrng.o micgetr.o param.o psif.o tdxin13d.o zcunix.o 
mpilib.o: cqlcomm.o param.o r8subs.o 
netcdfrf.o: bcast.o cqlcomm.o pack21.o param.o
netcdfrw2.o: advnce.o bcast.o cqlcomm.o coefefad.o coeffpad.o \
        coefmidt.o coefmidv.o coefrfad.o coefstup.o diagentr.o \
        param.o prppr.o r8subs.o tdfinterp.o zcunix.o
ntdstore.o: cqlcomm.o param.o restvty.o 
ntloop.o: cqlcomm.o param.o 
pack21.o: bcast.o r8subs.o 
pltdf.o: cqlcomm.o lookup.o param.o r8subs.o 
pltdnz.o: aminmx.o cqlcomm.o param.o pltelec.o 
pltelec.o: aminmx.o cqlcomm.o param.o 
pltendn.o: aminmx.o cqlcomm.o param.o r8subs.o 
pltfvsv.o: aminmx.o bcast.o cqlcomm.o param.o pltdf.o r8subs.o 
pltinit.o: cqlcomm.o param.o 
pltlosc.o: cqlcomm.o param.o pltdf.o 
pltmain.o:  advnce.o aminmx.o bcast.o coefefad.o coefegad.o coeffpad.o \
	coefload.o coefmidt.o coefmidv.o coefrfad.o coefstup.o coefsyad.o \
	cqlcomm.o fle.o diagentr.o param.o pltdf.o pltdnz.o pltendn.o \
	pltfvsv.o pltlosc.o pltpower.o  \
	pltvec.o pltvflux.o r8subs.o
pltpower.o: aminmx.o cqlcomm.o param.o 
pltrun.o: aminmx.o cqlcomm.o param.o r8subs.o tdnflxs.o 
pltvec.o: advnce.o aminmx.o bcast.o coefefad.o coeffpad.o \
	coefmidt.o coefmidv.o coefrfad.o coefstup.o cqlcomm.o diagentr.o \
	param.o pltvectr.o prppr.o r8subs.o 
pltvectr.o: aminmx.o r8subs.o 
pltvflux.o: aminmx.o cqlcomm.o param.o r8subs.o 
profaxis.o: cqlcomm.o param.o 
profiles.o: cfpmodbe.o cqlcomm.o param.o profaxis.o tdinterp.o \
	tdxin13d.o 
prppr.o: bcast.o cqlcomm.o diagwrng.o param.o r8subs.o 
prpprctr.o: aminmx.o cqlcomm.o param.o pltdf.o pltmain.o 
psif.o: cqlcomm.o diagwrng.o param.o zcunix.o 
r8lsode.o: r8subs.o
r8subs.o: cqlcomm.o 
rdc_bplt.o: bcast.o cqlcomm.o param.o pltdf.o 
rdc_multi.o: bcast.o cqlcomm.o netcdfrf.o param.o r8subs.o rdc_bplt.o \
	tdnflxs.o zcunix.o 
restcon.o: cqlcomm.o param.o 
resthks.o: cqlcomm.o param.o 
restvty.o: cqlcomm.o param.o 
rf.o: cqlcomm.o param.o 
sigalloc.o: bcast.o cqlcomm.o param.o 
siggy.o: cqlcomm.o param.o 
sigmax.o: bcast.o cqlcomm.o param.o 
sigmaxwl.o: bcast.o cqlcomm.o param.o 
sigsetup.o: bcast.o cqlcomm.o param.o r8subs.o sigfn.o siggy.o 
sigv.o: bcast.o cqlcomm.o param.o sigalloc.o sigsetup.o sigv5d.o \
	tdnflxs.o 
sigv5d.o: bcast.o cfpleg.o cqlcomm.o param.o r8subs.o sigmax.o 
soucrit.o: cfpgamma.o cqlcomm.o param.o r8subs.o 
sounorm.o: bcast.o cqlcomm.o param.o soup.o 
soup.o: cqlcomm.o param.o psif.o r8subs.o soup0.o 
soup0.o: cqlcomm.o param.o 
sourc0.o: cqlcomm.o param.o 
sourcee.o: bcast.o cqlcomm.o param.o sounorm.o sourc0.o sourcef.o \
	sourceko.o sourcpwr.o 
sourcef.o: bcast.o cqlcomm.o param.o r8subs.o soup.o 
sourceko.o: bcast.o cqlcomm.o fle.o lookup.o param.o r8subs.o \
	soucrit.o tdoutput.o 
sourcpwr.o: bcast.o cqlcomm.o param.o 
synchrad.o: bcast.o cqlcomm.o param.o 
tdboothi.o: bcast.o cqlcomm.o param.o zcunix.o 
tdbootst.o: bcast.o cqlcomm.o param.o zcunix.o 
tdchief.o: achief1.o achiefn.o aclear.o aindfpa.o ainplt.o \
	ainpltpa.o ainsetpa.o ampfar.o cfpgamma.o coefmidt.o coefmidv.o \
	coefstup.o coefwti.o coefwtj.o cqlcomm.o diag.o diagentr.o diaggnde.o \
	diagimpd.o diagscal.o dsk_gr.o dskout.o eflditer.o esefld.o \
	netcdfrf.o ntdstore.o param.o pltendn.o pltinit.o pltmain.o \
	pltrun.o profiles.o r8subs.o restvty.o sigv.o tddiag.o tdinitl.o \
	tdnflxs.o tdnpadiag.o tdoutput.o tdplteq.o tdpltmne.o tdsxray.o \
	tdtloop.o tdtoaray.o tdtransp.o tdtrcon.o tdtrdfus.o tdtrfcop.o \
	tdtrrsou.o tdtrsavf.o urfchief.o urfwrite.o wparsou.o wpavg.o \
	wpelecf.o wpsavf.o wptrafx.o wptramu.o 
tddiag.o: cqlcomm.o param.o tdboothi.o tdbootst.o 
tdeqdsk.o: cqlcomm.o equilib.o firstdrv.o param.o r8subs.o zcunix.o 
tdfinterp.o: cqlcomm.o lookup.o param.o 
tdinitl.o: ainalloc.o aindflt1.o aindflt.o aindfpa.o aingeom.o \
	ainitial.o ainpla.o ainsetva.o ainspec.o ainvnorm.o ampfar.o \
	baviorbt.o cqlcomm.o diagwrng.o eqindflt.o eqinitl.o micxinil.o \
	micxinit.o micxiniz.o param.o profiles.o rdc_multi.o sigv.o \
	tddiag.o tdeqdsk.o tdnflxs.o tdnpadiag.o tdoutput.o tdplteq.o \
	tdpltmne.o tdreadf.o tdrmshst.o tdstin.o tdsxray.o tdtoarad.o \
	tdtoaray.o tdtraloc.o tdtrdfus.o tdtrflg.o tdtrmuy.o tdtrvint.o \
	tdtrwtl.o tdwrng.o tdxinitl.o urfindfl.o urfinitl.o urfsetup.o \
	wpalloc.o wpchgdy.o wpinitl.o wploweq.o wpmshchk.o wptrmuy.o 
tdinlegw.o: bcast.o cqlcomm.o param.o 
tdinterp.o: zcunix.o 
tdnflxs.o: cqlcomm.o diagwrng.o param.o 
tdnpa.o: bcast.o cqlcomm.o param.o tdfinterp.o tdnpalam.o 
tdnpa0.o: bcast.o cqlcomm.o eqfpsi.o param.o tdnpa.o tdnpalam.o \
	tdsetnpa.o tdsxrplt.o zcunix.o 
tdnpadiag.o: cqlcomm.o param.o tdnpa0.o 
tdnpalam.o: cqlcomm.o param.o tdnpabscs.o 
tdoutput.o: bcast.o cfpgamma.o cqlcomm.o param.o restcon.o resthks.o \
	tdnflxs.o tdtrflx.o 
tdplteq.o: cqlcomm.o frplteq.o param.o 
tdpltjop.o: aminmx.o cqlcomm.o param.o 
tdpltmne.o: aminmx.o cqlcomm.o param.o r8subs.o tdpltjop.o 
tdreadf.o: bcast.o cqlcomm.o param.o tdnflxs.o zcunix.o 
tdrmshst.o: cqlcomm.o eqorbit.o eqvolpsi.o param.o zcunix.o 
tdsetnpa.o: bcast.o cqlcomm.o param.o zcunix.o 
tdsetsxr.o: bcast.o cqlcomm.o param.o tddsig.o 
tdstin.o: cqlcomm.o param.o 
tdsxr.o: bcast.o cfpleg.o cqlcomm.o param.o r8subs.o tdnflxs.o \
	tdwrng.o zcunix.o 
tdsxr0.o: bcast.o cqlcomm.o eqfpsi.o param.o tdsetsxr.o tdsxr.o \
	tdsxrplt.o zcunix.o 
tdsxray.o: cqlcomm.o param.o tdsxr0.o 
tdsxrplt.o: aminmx.o cqlcomm.o param.o r8subs.o 
tdtloop.o: aindfpa.o cqlcomm.o param.o tdeqdsk.o tdtscout.o tdwritef.o 
tdtoarad.o: cqlcomm.o param.o 
tdtoaray.o: cqlcomm.o param.o 
tdtraloc.o: bcast.o cqlcomm.o param.o 
tdtransp.o: cqlcomm.o param.o r8subs.o tdtravct.o tdtrchk.o \
	tdtrrtov2.o tdtrrtov.o tdtrsym.o tdtrvtor2.o tdtrvtor.o tdtrwtl.o 
tdtranspn.o: cqlcomm.o param.o r8subs.o tdtravct.o tdtrvtor.o \
	tdtrwtl.o 
tdtravct.o: bcast.o cqlcomm.o param.o r8subs.o tdtrrtov2.o tdtrrtov.o \
	tdtrsym.o tdtrvtor2.o tdtrvtor.o tdtrwtl.o 
tdtrchk.o: bcast.o cqlcomm.o diagwrng.o param.o r8subs.o 
tdtrchkd.o: bcast.o cqlcomm.o param.o 
tdtrcon.o: cqlcomm.o param.o tdtrflx.o 
tdtrdfus.o: bcast.o cqlcomm.o param.o tdnflxs.o 
tdtrfcop.o: cqlcomm.o param.o r8subs.o 
tdtrflg.o: bcast.o cqlcomm.o param.o 
tdtrflx.o: bcast.o cqlcomm.o param.o r8subs.o 
tdtrmuy.o: cqlcomm.o micgetr.o param.o 
tdtrrsou.o: cqlcomm.o param.o r8subs.o tdtravct.o tdtrrtov2.o \
	tdtrvtor2.o 
tdtrrtov.o: cqlcomm.o param.o tdtrchkd.o 
tdtrrtov2.o: cqlcomm.o param.o r8subs.o 
tdtrsavf.o: cqlcomm.o param.o 
tdtrsym.o: cqlcomm.o param.o 
tdtrvint.o: bcast.o cqlcomm.o param.o 
tdtrvsou.o: advnce.o cqlcomm.o diagentr.o param.o tdtrvtor3.o 
tdtrvtor.o: cqlcomm.o param.o r8subs.o tdtrchkd.o 
tdtrvtor2.o: cqlcomm.o param.o r8subs.o 
tdtrvtor3.o: cqlcomm.o param.o r8subs.o 
tdtrwtl.o: cqlcomm.o param.o 
tdtry.o: cqlcomm.o param.o tdwrng.o 
tdtscinp.o: cqlcomm.o param.o tdeqdsk.o 
tdtscout.o: cqlcomm.o param.o tdinterp.o 
tdwritef.o: cqlcomm.o param.o 
tdwrng.o: param.o 
tdxin13d.o: param.o profaxis.o 
tdxin23d.o: param.o profaxis.o 
tdxin33d.o: param.o profaxis.o 
tdxinitl.o: bcast.o cqlcomm.o micgetr.o param.o profaxis.o tdinterp.o \
	tdpro.o tdtscinp.o tdwrng.o tdxin13d.o tdxin23d.o tdxin33d.o 
urfalloc.o: bcast.o cqlcomm.o param.o 
urfavg.o: cqlcomm.o param.o 
urfb0.o: bcast.o cqlcomm.o param.o tdnflxs.o 
urfbes.o: aminmx.o bcast.o cqlcomm.o param.o zcunix.o 
urfbplt.o: bcast.o cqlcomm.o param.o pltdf.o 
urfchief.o: cqlcomm.o param.o tdnflxs.o urfavg.o urfb0.o urfbes.o \
	urfbplt.o urfdamp0.o urffflx.o urfpack.o urfrays.o urfread.o \
	urfwrite.o 
urfdamp0.o: bcast.o cqlcomm.o param.o r8subs.o urfdamp1.o urfdamp2.o \
	urfdampa.o 
urfdamp1.o: cqlcomm.o param.o tdnflxs.o 
urfdamp2.o: bcast.o cqlcomm.o param.o tdnflxs.o urfmidv.o 
urfdampa.o: cqlcomm.o param.o 
urfdout.o: cqlcomm.o param.o 
urfedge.o: cqlcomm.o param.o 
urffflx.o: bcast.o cqlcomm.o param.o r8subs.o zcunix.o 
urfindfl.o: cqlcomm.o param.o 
urfinitl.o: cqlcomm.o param.o 
urfmidv.o: bcast.o cqlcomm.o param.o r8subs.o 
urfpack.o: cqlcomm.o param.o r8subs.o tdnflxs.o urfedge.o urfwrong.o 
urfrays.o: cqlcomm.o param.o 
urfread.o: cqlcomm.o netcdfrf.o param.o r8subs.o urfread_.o 
urfread_.o: cqlcomm.o param.o urfwrong.o 
urfsetup.o: bcast.o cqlcomm.o netcdfrf.o param.o urfalloc.o urfread_.o 
urfwrite.o: cqlcomm.o param.o urfwrite_.o 
urfwrite_.o: cqlcomm.o param.o urfwr0.o urfwr0c.o 
urfwrong.o: param.o 
vlf.o: bcast.o cqlcomm.o param.o r8subs.o urfedge.o urfwrong.o \
	vlfbplt.o vlfsetup.o zcunix.o 
vlfalloc.o: bcast.o cqlcomm.o param.o 
vlfbplt.o: bcast.o cqlcomm.o param.o pltdf.o r8subs.o urfwrong.o 
vlfsetup.o: bcast.o cqlcomm.o param.o vlfalloc.o zcunix.o 
vlh.o: bcast.o cqlcomm.o param.o vlhbplt.o vlhd.o 
vlhbplt.o: bcast.o cqlcomm.o param.o pltdf.o 
vlhd.o: cqlcomm.o param.o 
wpalloc.o: bcast.o cqlcomm.o param.o 
wparsou.o: bcast.o cqlcomm.o param.o 
wpavg.o: cqlcomm.o param.o 
wpbdry.o: cqlcomm.o param.o r8subs.o 
wpcheck.o: bcast.o cqlcomm.o param.o 
wpchgdy.o: cqlcomm.o param.o 
wpcthta.o: cqlcomm.o param.o 
wpelecf.o: advnce.o cqlcomm.o param.o r8subs.o 
wpinitl.o: bcast.o cqlcomm.o param.o r8subs.o wpwrng.o 
wploweq.o: cqlcomm.o param.o 
wpmshchk.o: cqlcomm.o param.o 
wpsavf.o: cqlcomm.o param.o r8subs.o 
wptrafx.o: bcast.o cqlcomm.o param.o r8subs.o wpbdry.o wpcheck.o \
	wpwrng.o znonsym.o 
wptramu.o: bcast.o cqlcomm.o param.o r8subs.o wpbdry.o wpcheck.o \
	znonsym.o 
wptrmuy.o: cqlcomm.o param.o 
wpvptb.o: cqlcomm.o param.o 
wpwrng.o: cqlcomm.o param.o 
freya.o: aminmx.o bcast.o param.o
