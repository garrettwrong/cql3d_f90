c
c
      subroutine pltpower
      use param_mod
      use cqcomm_mod
      implicit integer (i-n), real*8 (a-h,o-z)
c
c     This routine plots the power due to the various physical mechanism
c     vs. time
c
      REAL RILIN
      REAL RPG1,RPG2
      REAL RNONCHA1(nonch),RNONCHA2(nonch)

c
      if (noplots.eq."enabled1") return

      do 100 k=1,ngen
        if (colmodl.eq.4 .and. k.eq.ngen) goto 100
c$$$        call gxglfr(0)
CPGPLT        CALL PGPAGE
        emin=ep90
        emax=-em90
        do 192 lu=-1,11
          call aminmx(pentr(1,k,lu,l_),1,nonch,1,hmin,hmax,idmin,idmax)
          if (hmin .lt. emin) emin=hmin
          if (hmax .gt. emax) emax=hmax
c     do 191 nc=1,nonch
c     h=pentr(nc,k,lu,l_)
c     if (h .lt. emin) emin=h
c     if (h .gt. emax) emax=h
c191  continue
 192    continue
        emax=emax+.05*(abs(emax))
        emin=emin-.05*(abs(emin))
CPGPLT        CALL PGSVP(.2,.8,.4,.95)
        DO I=1,NCH(L_)
           RNONCHA1(I)=ptime(i,l_)
        ENDDO
        RPG1=emin
        RPG2=emax
        IF ( RPG2-RPG1 .le. 1.e-16 ) THEN
           RPG2= RPG1+1.e-16
        ENDIF
CPGPLT        CALL PGSWIN(RNONCHA1(1),RNONCHA1(NCH(L_)),RPG1,RPG2)
CPGPLT        CALL PGBOX('BCNST',0.0,0,'BCNST',0.0,0)
CPGPLT        CALL PGLAB('Time (sec)','Power Den (W/cm\u3\d)',' ')
        ! Note: pentr(nch(l_),k,is,l_)=entr(k,is,l_)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,4,l_) !total - sum of all power sources
        ENDDO
CPGPLT        CALL PGSLS(1) ! 1-> solid
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)
        
        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,9,l_) !Power computed from df/dt. 
        ENDDO                          !Should equal entr(k,4,l_).
CPGPLT        CALL PGSLS(2) ! 2-> dashed
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,3,l_) !RF - radio frequency heating
        ENDDO
CPGPLT        CALL PGSLS(3) ! 3-> -.-.- 
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,2,l_) !ohmic - electric field term
        ENDDO
CPGPLT        CALL PGSLS(4) ! 4-> .....
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,1,l_) !collisional transfer from gens.
        ENDDO
CPGPLT        CALL PGSLS(1) ! 1-> solid
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,5,l_) !particle sources
        ENDDO
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,6,l_) !loss-lossmode
        ENDDO
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,7,l_) !losses-torloss
        ENDDO
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,-1,l_) !collisional transfer from Maxw. elec.
        ENDDO
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,0,l_) !collisional transfer from Maxw. ions
        ENDDO
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,8,l_) !losses due to runaway
        ENDDO
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)

        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,10,l_) !setting neg f to zero
        ENDDO
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)
        
        if (k .ne. kelecg .or. syncrad .eq. "disabled") go to 300
        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,11,l_) !synchrotron rad losses
        ENDDO
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)        
 300    continue
 
        if(tauegy(k,lr_).gt.0.)  then
        DO I=1,NCH(L_)
           RNONCHA2(I)=pentr(i,k,12,l_) !phenomenological energy losses
        ENDDO
CPGPLT        CALL PGLINE(nch(l_),RNONCHA1,RNONCHA2)        
        endif
        
c$$$ 200    format("COMPONENT BY COMPONENT POWER FLOW: SPECIES",i2,";",
c$$$     1    "ohmic--- electric field term",";",
c$$$     2    "col-mxe--- coulomb coll. with Maxwellian elec.",";",
c$$$     2    "col-mxi--- coulomb coll. with Maxwellian ions",";",
c$$$     3    "RF--- radio frequency heating",";",
c$$$     4    "total--- sum of all power sources ",";",
c$$$     5    "fusion--- fusion power output",";",
c$$$     6    "loss-o--- particle losses due to lossmode(k)",";",
c$$$     1    "loss-t--- losses due to torloss(k)",";",
c$$$     1    "runaway--- losses over high velocity terminator",";",
c$$$     1    "synchrd---synchrotron radiation losses",";",
c$$$     1    "los-egy---phenomenological energy loss","$")

        RILIN=4.
        write(t_,210) k ! Gen. species number
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,211) entr(k,4,l_),entr(k,9,l_) !sum over all components
                                        !and from df/dt (should be same)
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,212) entr(k,-1,l_) !collisional transfer from Maxwellian elec.
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,213) entr(k,0,l_)  !collisional transfer from Maxwellian ions
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,214) entr(k,1,l_)  !collisional transfer from gens.
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,215) entr(k,2,l_)  !ohmic drive
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,216) entr(k,3,l_)  !RF drive
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,217) entr(k,5,l_)  !particle sources
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,218) entr(k,6,l_),entr(k,7,l_) !loss-lossmode;losses-torloss
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,219) entr(k,8,l_)  !losses due to runaway
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,220) entr(k,10,l_) !setting neg f to zero
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_) 
        write(t_,221) entr(k,11,l_) !synchrotron rad losses
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)
        write(t_,222) entr(k,12,l_) !phenomenological energy losses
        RILIN=RILIN+1.
CPGPLT        CALL PGMTXT('B',RILIN,0.,0.,t_)

 210    format("Species k=",i2, "    Final powers in Watts/cc are:")
 211    format("sum over all comp=",1pe10.2, 3x,
     +         "From df/dt :",1pe10.2)
 212    format("collisional transfer from Maxwellian elec.=",1pe10.2)
 213    format("collisional transfer from Maxwellian ions=",1pe10.2)
 214    format("collisional transfer from gens.=",1pe10.2)
 215    format("ohmic drive=",1pe10.2)
 216    format("RF drive=",1pe10.2)
 217    format("particle sources=",1pe10.2)
 218    format("loss-lossmode(k)=", 1pe10.2, 3x,
     +         "losses-torloss(k)=",1pe10.2)
 219    format("losses due to runaway=",1pe10.2)
 220    format("setting neg f to zero=",1pe10.2)
 221    format("synchrotron rad losses=",1pe10.2)
 222    format("phenomenological energy losses=",1pe10.2)

 100  continue ! k=1,ngen
 
      return
      end
