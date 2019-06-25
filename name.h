!     name.h

!..................................................................
!     NAMELIST (SETUP) DECLARATION FOR INPUT
!..................................................................

      namelist/setup/ &
	acoefne,acoefte, &
	ampfmod,nampfmax,nonampf,ampferr,bctimescal, &
	bnumb,btor,bth,bootst,bootcalc,bootupdt,bootsign,nonboot,jhirsh, &
	contrmin,constr,chang,colmodl, &
	deltabdb,denpar,droptol,dtr,dtr1, &
	eegy,eparc,eperc,simpbfac,epar,eper, &
	elecfld,elpar0,enorm,enorme,enormi,eleccomp, &
	elecin,elecin_t,elecscal,enein,enein_t,ennin_t, &
	enescal,enloss,epsthet, &
	enmin,enmax,ennb,ennin,ennl,ennscal,enmin_npa,enmax_npa, &
	eseswtch,xsink,esink,ephicc,esfac,eoved, &
	fds,fds_npa,fmass,f4d_out, &
	tavg,tavg1,tavg2, &
	gsla,gslb,gamaset,gamafac,gamegy, &
	iactst,ineg,idskf,idskrf,ichkpnt,implct, &
	iprone,iprote,iproti,iprozeff,iprovphi,iproelec,ipronn,iprocur, &
	tmdmeth,isigmas,isigtst,isigsgv1,isigsgv2, &
	pltflux1, &
	irzplt,izeff,ioutime, &
	iy, &
	jx, &
	kenorm,lfil,kfrsou,kpress,kfield,kspeci,fpld, &
	lmidpln,locquas,lbdry,lbdry0,lossfile,lossmode,lmidvel,laddbnd, &
	lz, &
	machine,meshy,manymat,netcdfnm,netcdfshort, &
	netcdfvecal,netcdfvecc,netcdfvece,netcdfvecrf,netcdfvecs, &
	nnspec,mpwr,megy,mtorloss,mmsv,msxr,mx, &
	nchgdy,ngauss,nlagran, &
	nlotp1,nlotp2,nlotp3,nlotp4, &
	nmax,ngen,nkconro,nplt3d,nrskip,nen,nv,nen_npa,nv_npa, &
	npaproc,npa_process, &
	nr_delta,nz_delta,nt_delta, &
	nr_f4d,nz_f4d,nv_f4d,nt_f4d, &
	npwr,negy,ntorloss,njene,njte,njti, &
	nstop,nondtr1,nplot,nsave,ncoef,nchec,ncont,nrstrt,nstps,nfpld, &
	noncntrl,nonel,noffel,nonvphi,noffvphi,nonavgf,nofavgf, &
	nonloss,noffloss,nummods,numixts, &
	numby,negyrg, &
	oldiag, &
	plt3d,pltvs,partner,paregy,peregy,pegy, &
	zeffin,zeffin_t,zeffscal,vphiplin,vphiplin_t,vphiscal, &
	pltdn,pltvecal,pltvecc,pltvecrf,pltvece, &
	pltstrm,pltflux,pltmag,pltsig,pltlim,pltlimm, &
	pltrst,plturfb,pltvflu,pltra,pltfvs,pltd,pltprpp,pltfofv,pltlos, &
	pltrdc,profpsi, &
	psimodel,pltpowe,pltend,pltinput,pltview, &
	qsineut,trapmod,trapredc,scatmod,scatfrac, &
	ryain,radmaj,radmin,rmirror,relativ, &
	reden,regy,rfacz,rzset,rd,roveram, &
	rovera,rya,radcoord, &
	sbdry,scheck,ndeltarho,softxry,npa_diag,symtrap,syncrad, &
	bremsrad,brfac,brfac1,brfacgm3,sigmamod,sigvcx,sigvi, &
	soln_method,tauegy,taunew,tein,tein_t,tescal,tiin,tiin_t,tiscal, &
	tauloss,temp,temppar, &
	tfac,tfacz,tbnd,tandem, &
	thetd,torloss,thet1,thet2,x_sxr,z_sxr, &
	rd_npa,thetd_npa,x_npa,z_npa,thet1_npa,thet2_npa,atten_npa, &
	updown, &
	veclnth,vnorm, &
	xfac,xpctlwr,xpctmdl,xlwr,xmdl,xsinkm, &
	xprpmax,ipxy,jpxy, &
	yreset,ylower,yupper, &
	mpwrzeff,npwrzeff,mpwrvphi,npwrvphi,mpwrxj,npwrxj, &
	mpwrelec,npwrelec, &
	redenc,redenb,temp_den,tempc,tempb,zeffc,zeffb,elecc,elecb, &
	vphic,vphib,xjc,xjb,xjin_t,totcrt,efswtch,efswtchn, &
	efiter,efflag,curr_edge,efrelax,efrelax1,currerr, &
	bctime,nbctime, &
	zmax, &
	fow,outorb,nmu,npfi,nsteps_orb,nptsorb,i_orb_width,iorb2, &
	j0_ini,j0_end,inc_j0, i0_ini,i0_end,inc_i0, &
	j2_ini,j2_end,inc_j2, i2_ini,i2_end,inc_i2

!..................................................................
!     Namelist (sousetup) for "sou" simple source routines.
!..................................................................

      namelist/sousetup/ &
	asorz,asor,flemodel, &
	nonso,noffso,nso,nsou, &
	pltso,mpwrsou,npwrsou, &
	scm2z,szm1z,scm2,sellm1,sellm2,seppm1, &
	sellm1z,sellm2z,seppm2,sem1,sem2, &
	seppm1z,sem1z,sem2z,sthm1z, &
	seppm2z,soucoord,knockon,komodel,nkorfn,nonko,noffko,soffvte, &
	soffpr,isoucof,faccof,jfl,xlfac,xlpctlwr,xlpctmdl,xllwr,xlmdl, &
	szm2z,sthm1,szm1,szm2
