#------ Module generated dependencies, to be included in makefiles...

$(BUILDDIR)/a_cqlp.$(oext): $(BUILDDIR)/abchief.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/impavnc0.$(oext)
$(BUILDDIR)/abchief.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdchief.$(oext)
$(BUILDDIR)/achief1.$(oext): $(BUILDDIR)/achiefn.$(oext) $(BUILDDIR)/ainalloc.$(oext) $(BUILDDIR)/aindflt1.$(oext) $(BUILDDIR)/aindfpa.$(oext)  \
	$(BUILDDIR)/aingeom.$(oext) $(BUILDDIR)/ainitial.$(oext) $(BUILDDIR)/ainpla.$(oext) $(BUILDDIR)/ainsetva.$(oext) $(BUILDDIR)/ainspec.$(oext) $(BUILDDIR)/ainvnorm.$(oext)  \
	$(BUILDDIR)/cfpgamma.$(oext) $(BUILDDIR)/coefmidt.$(oext) $(BUILDDIR)/coefmidv.$(oext) $(BUILDDIR)/coefstup.$(oext) $(BUILDDIR)/coefwti.$(oext) $(BUILDDIR)/coefwtj.$(oext)  \
	$(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/diagimpd.$(oext) $(BUILDDIR)/diagscal.$(oext) $(BUILDDIR)/eqinitl.$(oext)  \
	$(BUILDDIR)/micxinit.$(oext) $(BUILDDIR)/micxiniz.$(oext) $(BUILDDIR)/ntloop.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltmain.$(oext) $(BUILDDIR)/profiles.$(oext) $(BUILDDIR)/r8subs.$(oext)  \
	$(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/wploweq.$(oext)
$(BUILDDIR)/achiefn.$(oext): $(BUILDDIR)/cfpcoefc.$(oext) $(BUILDDIR)/cfpgamma.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/cqldiag.$(oext)  \
	$(BUILDDIR)/diaggnde.$(oext) $(BUILDDIR)/dskout.$(oext) $(BUILDDIR)/efield.$(oext) $(BUILDDIR)/exsweep.$(oext) $(BUILDDIR)/finit.$(oext) $(BUILDDIR)/impavnc0.$(oext) $(BUILDDIR)/ntdstore.$(oext)  \
	$(BUILDDIR)/ntloop.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltmain.$(oext) $(BUILDDIR)/pltrun.$(oext) $(BUILDDIR)/restvty.$(oext) $(BUILDDIR)/sourcee.$(oext) $(BUILDDIR)/tdoutput.$(oext)
$(BUILDDIR)/aclear.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/advnce.$(oext): $(BUILDDIR)/bsl.$(oext) $(BUILDDIR)/bsu.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/ainalloc.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/aindflt1.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/aindfpa.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/aingeom.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/eqalloc.$(oext) $(BUILDDIR)/eqcoord.$(oext) $(BUILDDIR)/eqelpse.$(oext)  \
	$(BUILDDIR)/equilib.$(oext) $(BUILDDIR)/flxfn.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/ainitial.$(oext): $(BUILDDIR)/bavgmax.$(oext) $(BUILDDIR)/coefefld.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqldiag.$(oext) $(BUILDDIR)/diagentr.$(oext)  \
	$(BUILDDIR)/diaggnde.$(oext) $(BUILDDIR)/efield.$(oext) $(BUILDDIR)/finit.$(oext) $(BUILDDIR)/losscone.$(oext) $(BUILDDIR)/lossegy.$(oext) $(BUILDDIR)/micxinim.$(oext)  \
	$(BUILDDIR)/ntdstore.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltmain.$(oext) $(BUILDDIR)/restvty.$(oext) $(BUILDDIR)/sourcee.$(oext) $(BUILDDIR)/synchrad.$(oext)  \
	$(BUILDDIR)/tdoutput.$(oext) $(BUILDDIR)/vlf.$(oext) $(BUILDDIR)/vlh.$(oext) $(BUILDDIR)/wpvptb.$(oext)
$(BUILDDIR)/ainpla.$(oext): $(BUILDDIR)/cfpmodbe.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/ainplt.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/ainpltpa.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/ainsetpa.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/ainsetva.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/wpwrng.$(oext)
$(BUILDDIR)/ainspec.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/ainvnorm.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/ampfar.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/bavdens.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/exlin.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/psif.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/bavgmax.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/psif.$(oext)
$(BUILDDIR)/baviorbt.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/lookup.$(oext) $(BUILDDIR)/micgnbnd.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/psif.$(oext)  \
	$(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/bsl.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/lookup.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/bsu.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/lookup.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/cfpcoefc.$(oext): $(BUILDDIR)/cfpcoefn.$(oext) $(BUILDDIR)/cfpcoefr.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/cfpcoefn.$(oext): $(BUILDDIR)/bavdens.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cfpleg.$(oext) $(BUILDDIR)/cfpmodbe.$(oext) $(BUILDDIR)/cfpsymt.$(oext) $(BUILDDIR)/cqlcomm.$(oext)  \
	$(BUILDDIR)/diagescl.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/cfpcoefr.$(oext): $(BUILDDIR)/bavdens.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cfpleg.$(oext) $(BUILDDIR)/cfpmodbe.$(oext) $(BUILDDIR)/cfpsymt.$(oext) $(BUILDDIR)/cqlcomm.$(oext)  \
	$(BUILDDIR)/diagescl.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/cfpgamma.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/cfpleg.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/cfpsymt.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coefefad.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coefefld.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coefegad.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coeffpad.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coefload.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/losstor.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coefmidt.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coefmidv.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/coefrfad.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coefstup.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/coefefad.$(oext) $(BUILDDIR)/coefegad.$(oext) $(BUILDDIR)/coeffpad.$(oext) $(BUILDDIR)/coefload.$(oext)  \
	$(BUILDDIR)/coefrfad.$(oext) $(BUILDDIR)/coefsyad.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/coefsyad.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coefwti.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/coefmidt.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/coefwtj.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/coefmidv.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/cqlcomm.$(oext): $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/cqlconf.$(oext): $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/cqldiag.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/diagcfac.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/diagdens.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/diagdenz.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cfpleg.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/diagentr.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/coefefad.$(oext) $(BUILDDIR)/coefegad.$(oext) $(BUILDDIR)/coefload.$(oext)  \
	$(BUILDDIR)/coefmidv.$(oext) $(BUILDDIR)/coefrfad.$(oext) $(BUILDDIR)/coefsyad.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext)  \
	$(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/diagescl.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/diaggnde.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/diagcfac.$(oext) $(BUILDDIR)/diagdenz.$(oext)  \
	$(BUILDDIR)/diaggnde2.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/soucrit.$(oext)
$(BUILDDIR)/diaggnde2.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cfpleg.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/diagcfac.$(oext)  \
	$(BUILDDIR)/diagdenz.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/diagimpd.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagdens.$(oext) $(BUILDDIR)/diagentr.$(oext)  \
	$(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/diagscal.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/diagwrng.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/diagxswt.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagdens.$(oext) $(BUILDDIR)/diagentr.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/diagxswx.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/coefmidt.$(oext) $(BUILDDIR)/coefmidv.$(oext) $(BUILDDIR)/coefstup.$(oext)  \
	$(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagdens.$(oext) $(BUILDDIR)/diagentr.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/dsk_gr.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/dskout.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/efield.$(oext): $(BUILDDIR)/cfpgamma.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/restcon.$(oext) $(BUILDDIR)/resthks.$(oext)
$(BUILDDIR)/eflditer.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/coefefad.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/restcon.$(oext)  \
	$(BUILDDIR)/resthks.$(oext) $(BUILDDIR)/soucrit.$(oext)
$(BUILDDIR)/eqalloc.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqcoord.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqfndpsi.$(oext) $(BUILDDIR)/eqrhopsi.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqelpse.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqfn.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqflxavg.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqorbit.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqfn.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqwrng.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqfndpsi.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqflxavg.$(oext) $(BUILDDIR)/eqfpsi.$(oext) $(BUILDDIR)/eqonovrp.$(oext) $(BUILDDIR)/eqorbit.$(oext)  \
	$(BUILDDIR)/eqvolpsi.$(oext) $(BUILDDIR)/eqwrng.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqfninv.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqwrng.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqfpsi.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqwrng.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/eqinitl.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqjac.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/eqonovrp.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqflxavg.$(oext) $(BUILDDIR)/eqorbit.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqorbit.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqfpsi.$(oext) $(BUILDDIR)/eqjac.$(oext) $(BUILDDIR)/eqrhs.$(oext)  \
	$(BUILDDIR)/eqwrng.$(oext) $(BUILDDIR)/exlin.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/eqrhopsi.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqflxavg.$(oext) $(BUILDDIR)/eqfpsi.$(oext) $(BUILDDIR)/eqonovrp.$(oext)  \
	$(BUILDDIR)/eqorbit.$(oext) $(BUILDDIR)/eqvolpsi.$(oext) $(BUILDDIR)/eqwrng.$(oext) $(BUILDDIR)/exlin.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/eqrhs.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/eqtopeol.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/equilib.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqtopeol.$(oext) $(BUILDDIR)/eqwrng.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/eqvolpsi.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqorbit.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/eqwrng.$(oext): $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/esefld.$(oext): $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/exsweep.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/coefmidt.$(oext) $(BUILDDIR)/coefmidv.$(oext) $(BUILDDIR)/coefstup.$(oext) $(BUILDDIR)/coefwti.$(oext)  \
	$(BUILDDIR)/coefwtj.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagxswt.$(oext) $(BUILDDIR)/diagxswx.$(oext) $(BUILDDIR)/exsweept.$(oext) $(BUILDDIR)/exsweepx.$(oext)  \
	$(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/exsweept.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/exsweepx.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagentr.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/finit.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/lossorbm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdreadf.$(oext)
$(BUILDDIR)/fle.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/micgetr.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/flxfn.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/freya.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zfreya.$(oext)
$(BUILDDIR)/freyasou.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/eqfpsi.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/sourcpwr.$(oext)  \
	$(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/tdtoaray.$(oext) $(BUILDDIR)/urfb0.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/frinitl.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/frinitz.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/frnbdep2.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/frnfreya.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/frplteq.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/frplteq.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdnflxs.$(oext)
$(BUILDDIR)/frset.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zfreya.$(oext)
$(BUILDDIR)/frsmooth.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/frsplft.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/frsuppor.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/frstup.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/frwrong.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/ilut.$(oext): $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/impavnc0.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/bsl.$(oext) $(BUILDDIR)/bsu.$(oext) $(BUILDDIR)/coefmidt.$(oext) $(BUILDDIR)/coefmidv.$(oext)  \
	$(BUILDDIR)/coefstup.$(oext) $(BUILDDIR)/coefwti.$(oext) $(BUILDDIR)/coefwtj.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/esefld.$(oext) $(BUILDDIR)/ilut.$(oext)  \
	$(BUILDDIR)/impchk.$(oext) $(BUILDDIR)/impnorm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdtranspn.$(oext) $(BUILDDIR)/tdtrvsou.$(oext)
$(BUILDDIR)/impchk.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagentr.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/lookup.$(oext): $(BUILDDIR)/urfb0.$(oext)
$(BUILDDIR)/losscone.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/lossorbm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/lossegy.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/lossorbm.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/losstor.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/micgetr.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/micgmbnd.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/psif.$(oext)
$(BUILDDIR)/micgnbnd.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/micgmbnd.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/micxinil.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/psif.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdinlegw.$(oext)  \
	$(BUILDDIR)/wpwrng.$(oext)
$(BUILDDIR)/micxinim.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/micxinit.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/micgetr.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdtry.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/micxiniz.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/micgetr.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/psif.$(oext) $(BUILDDIR)/tdxin13d.$(oext)  \
	$(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/cql3d_mpilib.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/netcdfrf.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/pack21.$(oext) $(BUILDDIR)/netcdfrw2.$(oext)  $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/netcdfrw2.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/coeffpad.$(oext) $(BUILDDIR)/coefmidt.$(oext) $(BUILDDIR)/coefstup.$(oext)  \
	$(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagentr.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/prppr.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdfinterp.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/ntdstore.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/restvty.$(oext)
$(BUILDDIR)/ntloop.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/pack21.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/pltdf.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/lookup.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/pltdnz.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltelec.$(oext)
$(BUILDDIR)/pltelec.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/pltendn.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/pltfvsv.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltdf.$(oext)  \
	$(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/pltinit.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/netcdfrw2.$(oext)
$(BUILDDIR)/pltlosc.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltdf.$(oext)
$(BUILDDIR)/pltmain.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/coefefad.$(oext) $(BUILDDIR)/coefegad.$(oext) $(BUILDDIR)/coeffpad.$(oext)  \
	$(BUILDDIR)/coefload.$(oext) $(BUILDDIR)/coefmidt.$(oext) $(BUILDDIR)/coefmidv.$(oext) $(BUILDDIR)/coefrfad.$(oext) $(BUILDDIR)/coefstup.$(oext) $(BUILDDIR)/coefsyad.$(oext)  \
	$(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/diagentr.$(oext) $(BUILDDIR)/fle.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltdf.$(oext) $(BUILDDIR)/pltdnz.$(oext)  \
	$(BUILDDIR)/pltendn.$(oext) $(BUILDDIR)/pltfvsv.$(oext) $(BUILDDIR)/pltlosc.$(oext) $(BUILDDIR)/pltpower.$(oext) $(BUILDDIR)/pltvec.$(oext) $(BUILDDIR)/pltvflux.$(oext)  \
	$(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/pltpower.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/pltrun.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdnflxs.$(oext)
$(BUILDDIR)/pltvec.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/coefefad.$(oext) $(BUILDDIR)/coeffpad.$(oext) $(BUILDDIR)/coefmidt.$(oext)  \
	$(BUILDDIR)/coefmidv.$(oext) $(BUILDDIR)/coefrfad.$(oext) $(BUILDDIR)/coefstup.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/diagentr.$(oext)  \
	$(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltvectr.$(oext) $(BUILDDIR)/prppr.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/pltvectr.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/pltvflux.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/profaxis.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/profiles.$(oext): $(BUILDDIR)/cfpmodbe.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/profaxis.$(oext) $(BUILDDIR)/tdinterp.$(oext)  \
	$(BUILDDIR)/tdxin13d.$(oext)
$(BUILDDIR)/prppr.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/prpprctr.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltdf.$(oext) $(BUILDDIR)/pltmain.$(oext)
$(BUILDDIR)/psif.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/r8lsode.$(oext): $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/r8subs.$(oext): $(BUILDDIR)/cqlcomm.$(oext)
$(BUILDDIR)/rdc_bplt.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltdf.$(oext)
$(BUILDDIR)/rdc_multi.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/netcdfrf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/rdc_bplt.$(oext)  \
	$(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/restcon.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/resthks.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/restvty.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/rf.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/sigalloc.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/siggy.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/sigmax.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/sigmaxwl.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/sigsetup.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/sigfn.$(oext) $(BUILDDIR)/siggy.$(oext)
$(BUILDDIR)/sigv.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/sigalloc.$(oext) $(BUILDDIR)/sigsetup.$(oext) $(BUILDDIR)/sigv5d.$(oext)  \
	$(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdnflxs.$(oext)
$(BUILDDIR)/sigv5d.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cfpleg.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/sigmax.$(oext)
$(BUILDDIR)/soucrit.$(oext): $(BUILDDIR)/cfpgamma.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/sounorm.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/soup.$(oext)
$(BUILDDIR)/soup.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/psif.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/soup0.$(oext)
$(BUILDDIR)/soup0.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/sourc0.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/sourcee.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/sounorm.$(oext) $(BUILDDIR)/sourc0.$(oext) $(BUILDDIR)/sourcef.$(oext)  \
	$(BUILDDIR)/sourceko.$(oext) $(BUILDDIR)/sourcpwr.$(oext)
$(BUILDDIR)/sourcef.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/soup.$(oext)
$(BUILDDIR)/sourceko.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/fle.$(oext) $(BUILDDIR)/lookup.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)  \
	$(BUILDDIR)/soucrit.$(oext) $(BUILDDIR)/tdoutput.$(oext)
$(BUILDDIR)/sourcpwr.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/synchrad.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdboothi.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdbootst.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdchief.$(oext): $(BUILDDIR)/achief1.$(oext) $(BUILDDIR)/achiefn.$(oext) $(BUILDDIR)/aclear.$(oext) $(BUILDDIR)/aindfpa.$(oext) $(BUILDDIR)/ainplt.$(oext) $(BUILDDIR)/ainpltpa.$(oext)  \
	$(BUILDDIR)/ainsetpa.$(oext) $(BUILDDIR)/ampfar.$(oext) $(BUILDDIR)/cfpgamma.$(oext) $(BUILDDIR)/coefmidt.$(oext) $(BUILDDIR)/coefmidv.$(oext) $(BUILDDIR)/coefstup.$(oext)  \
	$(BUILDDIR)/coefwti.$(oext) $(BUILDDIR)/coefwtj.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/cqldiag.$(oext) $(BUILDDIR)/diagentr.$(oext)  \
	$(BUILDDIR)/diaggnde.$(oext) $(BUILDDIR)/diagimpd.$(oext) $(BUILDDIR)/diagscal.$(oext) $(BUILDDIR)/dsk_gr.$(oext) $(BUILDDIR)/dskout.$(oext) $(BUILDDIR)/eflditer.$(oext)  \
	$(BUILDDIR)/esefld.$(oext) $(BUILDDIR)/netcdfrf.$(oext) $(BUILDDIR)/ntdstore.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltendn.$(oext) $(BUILDDIR)/pltinit.$(oext) $(BUILDDIR)/pltmain.$(oext)  \
	$(BUILDDIR)/pltrun.$(oext) $(BUILDDIR)/profiles.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/restvty.$(oext) $(BUILDDIR)/sigv.$(oext) $(BUILDDIR)/tddiag.$(oext) $(BUILDDIR)/tdinitl.$(oext)  \
	$(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/tdnpadiag.$(oext) $(BUILDDIR)/tdoutput.$(oext) $(BUILDDIR)/tdplteq.$(oext) $(BUILDDIR)/tdpltmne.$(oext) $(BUILDDIR)/tdsxray.$(oext)  \
	$(BUILDDIR)/tdtloop.$(oext) $(BUILDDIR)/tdtoaray.$(oext) $(BUILDDIR)/tdtransp.$(oext) $(BUILDDIR)/tdtrcon.$(oext) $(BUILDDIR)/tdtrdfus.$(oext) $(BUILDDIR)/tdtrfcop.$(oext)  \
	$(BUILDDIR)/tdtrrsou.$(oext) $(BUILDDIR)/tdtrsavf.$(oext) $(BUILDDIR)/urfchief.$(oext) $(BUILDDIR)/urfwrite.$(oext) $(BUILDDIR)/wparsou.$(oext) $(BUILDDIR)/wpavg.$(oext)  \
	$(BUILDDIR)/wpelecf.$(oext) $(BUILDDIR)/wpsavf.$(oext) $(BUILDDIR)/wptrafx.$(oext) $(BUILDDIR)/wptramu.$(oext)
$(BUILDDIR)/tddiag.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdboothi.$(oext) $(BUILDDIR)/tdbootst.$(oext)
$(BUILDDIR)/tdeqdsk.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/equilib.$(oext) $(BUILDDIR)/firstdrv.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdfinterp.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/lookup.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdinitl.$(oext): $(BUILDDIR)/ainalloc.$(oext) $(BUILDDIR)/aindflt1.$(oext) $(BUILDDIR)/aindfpa.$(oext) $(BUILDDIR)/aingeom.$(oext)  \
	$(BUILDDIR)/ainitial.$(oext) $(BUILDDIR)/ainpla.$(oext) $(BUILDDIR)/ainsetva.$(oext) $(BUILDDIR)/ainspec.$(oext) $(BUILDDIR)/ainvnorm.$(oext) $(BUILDDIR)/ampfar.$(oext)  \
	$(BUILDDIR)/baviorbt.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/eqinitl.$(oext)  \
	$(BUILDDIR)/micxinil.$(oext) $(BUILDDIR)/micxinit.$(oext) $(BUILDDIR)/micxiniz.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/profiles.$(oext) $(BUILDDIR)/rdc_multi.$(oext)  \
	$(BUILDDIR)/sigv.$(oext) $(BUILDDIR)/tddiag.$(oext) $(BUILDDIR)/tdeqdsk.$(oext) $(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/tdnpadiag.$(oext) $(BUILDDIR)/tdoutput.$(oext) $(BUILDDIR)/tdplteq.$(oext)  \
	$(BUILDDIR)/tdpltmne.$(oext) $(BUILDDIR)/tdreadf.$(oext) $(BUILDDIR)/tdrmshst.$(oext) $(BUILDDIR)/tdstin.$(oext) $(BUILDDIR)/tdsxray.$(oext) $(BUILDDIR)/tdtoarad.$(oext)  \
	$(BUILDDIR)/tdtoaray.$(oext) $(BUILDDIR)/tdtraloc.$(oext) $(BUILDDIR)/tdtrdfus.$(oext) $(BUILDDIR)/tdtrflg.$(oext) $(BUILDDIR)/tdtrmuy.$(oext) $(BUILDDIR)/tdtrvint.$(oext)  \
	$(BUILDDIR)/tdtrwtl.$(oext) $(BUILDDIR)/tdwrng.$(oext) $(BUILDDIR)/tdxinitl.$(oext) $(BUILDDIR)/urfinitl.$(oext) $(BUILDDIR)/urfsetup.$(oext)  \
	$(BUILDDIR)/wpalloc.$(oext) $(BUILDDIR)/wpchgdy.$(oext) $(BUILDDIR)/wpinitl.$(oext) $(BUILDDIR)/wploweq.$(oext) $(BUILDDIR)/wpmshchk.$(oext) $(BUILDDIR)/wptrmuy.$(oext)
$(BUILDDIR)/tdinlegw.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdinterp.$(oext): $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdnflxs.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdnpa.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdfinterp.$(oext) $(BUILDDIR)/tdnpalam.$(oext)
$(BUILDDIR)/tdnpa0.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/eqfpsi.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdnpa.$(oext)  \
	$(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdnpalam.$(oext) $(BUILDDIR)/tdsetnpa.$(oext) $(BUILDDIR)/tdsxrplt.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdnpadiag.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdnpa0.$(oext)
$(BUILDDIR)/tdnpalam.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdnpabscs.$(oext)
$(BUILDDIR)/tdoutput.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cfpgamma.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/netcdfrw2.$(oext) $(BUILDDIR)/param.$(oext) \
	$(BUILDDIR)/restcon.$(oext) $(BUILDDIR)/resthks.$(oext) $(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/tdtrflx.$(oext)
$(BUILDDIR)/tdplteq.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/frplteq.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdpltjop.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdpltmne.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdpltjop.$(oext)
$(BUILDDIR)/tdreadf.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/netcdfrw2.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdnflxs.$(oext) \
	$(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdrmshst.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/eqorbit.$(oext) $(BUILDDIR)/eqvolpsi.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdsetnpa.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdsetsxr.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tddsig.$(oext)
$(BUILDDIR)/tdstin.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdsxr.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cfpleg.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)  \
	$(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/tdwrng.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdsxr0.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/eqfpsi.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdsetsxr.$(oext)  \
	$(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdsxr.$(oext) $(BUILDDIR)/tdsxrplt.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/tdsxray.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdsxr0.$(oext)
$(BUILDDIR)/tdsxrplt.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/tdtloop.$(oext): $(BUILDDIR)/aindfpa.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdeqdsk.$(oext) $(BUILDDIR)/tdtscout.$(oext)  \
	$(BUILDDIR)/tdwritef.$(oext)
$(BUILDDIR)/tdtoarad.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtoaray.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtraloc.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtransp.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdtravct.$(oext) $(BUILDDIR)/tdtrchk.$(oext)  \
	$(BUILDDIR)/tdtrrtov.$(oext) $(BUILDDIR)/tdtrrtov2.$(oext) $(BUILDDIR)/tdtrsym.$(oext) $(BUILDDIR)/tdtrvtor.$(oext) $(BUILDDIR)/tdtrvtor2.$(oext) $(BUILDDIR)/tdtrwtl.$(oext)
$(BUILDDIR)/tdtranspn.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdtravct.$(oext) $(BUILDDIR)/tdtrvtor.$(oext)  \
	$(BUILDDIR)/tdtrwtl.$(oext)
$(BUILDDIR)/tdtravct.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdtrrtov.$(oext)  \
	$(BUILDDIR)/tdtrrtov2.$(oext) $(BUILDDIR)/tdtrsym.$(oext) $(BUILDDIR)/tdtrvtor.$(oext) $(BUILDDIR)/tdtrvtor2.$(oext) $(BUILDDIR)/tdtrwtl.$(oext)
$(BUILDDIR)/tdtrchk.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/diagwrng.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/tdtrchkd.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtrcon.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdtrflx.$(oext)
$(BUILDDIR)/tdtrdfus.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/netcdfrw2.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdnflxs.$(oext)
$(BUILDDIR)/tdtrfcop.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/tdtrflg.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtrflx.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/tdtrmuy.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/micgetr.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtrrsou.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdtravct.$(oext) $(BUILDDIR)/tdtrrtov2.$(oext)  \
	$(BUILDDIR)/tdtrvtor2.$(oext)
$(BUILDDIR)/tdtrrtov.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdtrchkd.$(oext)
$(BUILDDIR)/tdtrrtov2.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/tdtrsavf.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtrsym.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtrvint.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtrvsou.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/diagentr.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdtrvtor3.$(oext)
$(BUILDDIR)/tdtrvtor.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdtrchkd.$(oext)
$(BUILDDIR)/tdtrvtor2.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/tdtrvtor3.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/tdtrwtl.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdtry.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdwrng.$(oext)
$(BUILDDIR)/tdtscinp.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdeqdsk.$(oext)
$(BUILDDIR)/tdtscout.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdinterp.$(oext)
$(BUILDDIR)/tdwritef.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdwrng.$(oext): $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/tdxin13d.$(oext): $(BUILDDIR)/param.$(oext) $(BUILDDIR)/profaxis.$(oext)
$(BUILDDIR)/tdxin23d.$(oext): $(BUILDDIR)/param.$(oext) $(BUILDDIR)/profaxis.$(oext)
$(BUILDDIR)/tdxin33d.$(oext): $(BUILDDIR)/param.$(oext) $(BUILDDIR)/profaxis.$(oext)
$(BUILDDIR)/tdxinitl.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/micgetr.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/profaxis.$(oext)  \
	$(BUILDDIR)/tdinterp.$(oext) $(BUILDDIR)/tdpro.$(oext) $(BUILDDIR)/tdtscinp.$(oext) $(BUILDDIR)/tdwrng.$(oext) $(BUILDDIR)/tdxin13d.$(oext) $(BUILDDIR)/tdxin23d.$(oext)  \
	$(BUILDDIR)/tdxin33d.$(oext)
$(BUILDDIR)/urfalloc.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/urfavg.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/urfb0.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdnflxs.$(oext)
$(BUILDDIR)/urfbes.$(oext): $(BUILDDIR)/aminmx.$(oext) $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/urfbplt.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltdf.$(oext)
$(BUILDDIR)/urfchief.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/urfavg.$(oext) $(BUILDDIR)/urfb0.$(oext)  \
	$(BUILDDIR)/urfbes.$(oext) $(BUILDDIR)/urfbplt.$(oext) $(BUILDDIR)/urfdamp0.$(oext) $(BUILDDIR)/urffflx.$(oext) $(BUILDDIR)/urfpack.$(oext) $(BUILDDIR)/urfrays.$(oext)  \
	$(BUILDDIR)/urfread.$(oext) $(BUILDDIR)/urfwrite.$(oext)
$(BUILDDIR)/urfdamp0.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/urfdamp1.$(oext) $(BUILDDIR)/urfdamp2.$(oext)  \
	$(BUILDDIR)/urfdampa.$(oext)
$(BUILDDIR)/urfdamp1.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdnflxs.$(oext)
$(BUILDDIR)/urfdamp2.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/urfmidv.$(oext)
$(BUILDDIR)/urfdampa.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/urfdout.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/urfedge.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/urffflx.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/urfinitl.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/netcdfrw2.$(oext)
$(BUILDDIR)/urfmidv.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/urfpack.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/tdnflxs.$(oext) $(BUILDDIR)/urfedge.$(oext) $(BUILDDIR)/urfwrong.$(oext)
$(BUILDDIR)/urfrays.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/urfread.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/netcdfrf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/urfread_.$(oext)
$(BUILDDIR)/urfread_.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/urfwrong.$(oext)
$(BUILDDIR)/urfsetup.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/netcdfrf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/urfalloc.$(oext)  \
	$(BUILDDIR)/urfread_.$(oext)
$(BUILDDIR)/urfwrite.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/urfwrite_.$(oext)
$(BUILDDIR)/urfwrite_.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/urfwr0.$(oext) $(BUILDDIR)/urfwr0c.$(oext)
$(BUILDDIR)/urfwrong.$(oext): $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/vlf.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/urfedge.$(oext)  \
	$(BUILDDIR)/urfwrong.$(oext) $(BUILDDIR)/vlfbplt.$(oext) $(BUILDDIR)/vlfsetup.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/vlfalloc.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/vlfbplt.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltdf.$(oext) $(BUILDDIR)/r8subs.$(oext)  \
	$(BUILDDIR)/urfwrong.$(oext)
$(BUILDDIR)/vlfsetup.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/vlfalloc.$(oext) $(BUILDDIR)/zcunix.$(oext)
$(BUILDDIR)/vlh.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/vlhbplt.$(oext) $(BUILDDIR)/vlhd.$(oext)
$(BUILDDIR)/vlhbplt.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/pltdf.$(oext)
$(BUILDDIR)/vlhd.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpalloc.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wparsou.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpavg.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpbdry.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/wpcheck.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpchgdy.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpcthta.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpelecf.$(oext): $(BUILDDIR)/advnce.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/wpinitl.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/wpwrng.$(oext)
$(BUILDDIR)/wploweq.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpmshchk.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpsavf.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext)
$(BUILDDIR)/wptrafx.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/wpbdry.$(oext) $(BUILDDIR)/wpcheck.$(oext)  \
	$(BUILDDIR)/wpwrng.$(oext) $(BUILDDIR)/znonsym.$(oext)
$(BUILDDIR)/wptramu.$(oext): $(BUILDDIR)/bcast.$(oext) $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext) $(BUILDDIR)/r8subs.$(oext) $(BUILDDIR)/wpbdry.$(oext) $(BUILDDIR)/wpcheck.$(oext)  \
	$(BUILDDIR)/znonsym.$(oext)
$(BUILDDIR)/wptrmuy.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpvptb.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/wpwrng.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/param.$(oext)
$(BUILDDIR)/zcunix.$(oext): $(BUILDDIR)/cqlconf.$(oext)
$(BUILDDIR)/znonsym.$(oext): $(BUILDDIR)/cqlconf.$(oext)
$(BUILDDIR)/zfreya.$(oext): $(BUILDDIR)/cqlcomm.$(oext) $(BUILDDIR)/cqlconf.$(oext) $(BUILDDIR)/frsubs.$(oext) $(BUILDDIR)/param.$(oext)
