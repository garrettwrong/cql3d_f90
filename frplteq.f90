!
module frplteq_mod

  !---BEGIN USE

  use r8subs_mod, only : dcopy
  use tdnflxs_mod, only : tdnflxs

  !---END USE
  character(len=8), public :: textt(200)

  save

contains

  subroutine frplteq(xpts,ypts,zpts,rpts,vx,vy,vz,ipts,curdep,nfrplt,frplt)
    use param_mod
    use comm_mod
    use r8subs_mod, only : dcopy
    implicit integer (i-n), real*8 (a-h,o-z)
    save
    !MPIINSERT_INCLUDE


    character(len=8) frplt

    dimension xpts(*),ypts(*),zpts(*),rpts(*),vx(*),vy(*),vz(*)

    REAL*4 RBOT,RTOP,ZBOT,ZTOP
    REAL*4 RTAB1(LFIELDA),RTAB2(LFIELDA)
    REAL*4 RPG1,RPG2, xyplotmax

    data nconskp /2/

    !..................................................................
    !     frplt="enabled":
    !     This routine plots out the contours (flux surfaces) on
    !     which the calculations are done and it plots out some of the
    !     birth points of the neutral beam source.
    !     frplt="plotwrit":
    !     Also outputs points to ascii file freya_points.txt.
    !     frplt="write_op":
    !     No plotting, just output to ascii file freya_points.txt.
    !..................................................................


    !MPIINSERT_IF_RANK_NE_0_RETURN
    ! make plots on mpirank.eq.0 only

    if (noplots.eq."enabled1") return

    if (frplt.eq."disabled") return

    if (frplt.eq."write_op") goto 200 !Skip plotting but save data into freya_points.txt

    call micfrplt

    !---- PLOTS in (R,Z) ---------------------------------------------------
    rmincon1=rmincon
    if(machine.eq."mirror") rmincon1=-rmaxcon

    if (zmaxcon.gt..5*(rmaxcon-rmincon1)) then
       delr=(rmaxcon-rmincon1)*.9/(2.*zmaxcon)
       zbot=.01 !.05
       ztop=.99 !.95
       rbot=.01+(.9-delr)/2.
       rtop=.99-(.9-delr)/2.
    else
       delz=(2.*zmaxcon)*.9/(rmaxcon-rmincon1)
       rbot=.01
       rtop=.99
       zbot=.01+(.9-delz)/2.
       ztop=.99-(.9-delz)/2.
    endif

          CALL PGPAGE
          CALL PGSVP(rbot,rtop,zbot,ztop)
    RBOT=rmincon
    RTOP=rmaxcon*1.2 ! give 20% more outside of last surface
    if(machine.eq."mirror") then
       RBOT=-RTOP ! to plot left and right sides of flux surfaces
    endif
    ZBOT=zmincon
    ZTOP=zmaxcon
          CALL PGSWIN(rbot,rtop,zbot,ztop)
          CALL PGWNAD(rbot,rtop,zbot,ztop)  ! limits
          CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
          if(machine.eq."mirror") then
          CALL PGLAB('X (cms)','Z (cms)', 'NBI Deposition')
          else
          CALL PGLAB('Major radius (cms)','Vert height (cms)', 'NBI Deposition')
          endif
    xyplotmax=0. ! to set limits in (X,Y) plots
    do 10 l=1,lrz,nconskp
       l1=l
       if (l1.gt.lrz-nconskp) l1=lrz
       call tdnflxs(l1)
       call dcopy(lfield,solr(1:lfield,lr_),1,tlorb1,1)
       call dcopy(lfield,solz(1:lfield,lr_),1,tlorb2,1)

       do 20 j=1,lorbit(lr_)
          solr(j,lr_)=tlorb1(lorbit(lr_)+1-j)
          solz(j,lr_)=tlorb2(lorbit(lr_)+1-j)
20     end do
       do j=1,lorbit(lr_)
          RTAB1(j)=solr(j,lr_)
          RTAB2(j)=solz(j,lr_)
       enddo
               CALL PGLINE(LORBIT(LR_),RTAB1,RTAB2)
       if(machine.eq."mirror") then
                  CALL PGLINE(LORBIT(LR_),-RTAB1,RTAB2) !mirror area to the left of Z-axis
       endif

       !       if eqsym.ne."none", still need to plot lower flux surface
       if (eqsym.ne."none") then
          do 30 j=1,lorbit(lr_)
             solz(j,lr_)=-solz(j,lr_)
30        end do
          DO J=1,LORBIT(LR_)
             RTAB2(J)=SOLZ(J,LR_)
          ENDDO
                     CALL PGLINE(LORBIT(LR_),RTAB1,RTAB2)
          if(machine.eq."mirror") then
                        CALL PGLINE(LORBIT(LR_),-RTAB1,RTAB2) !mirror area to the left of Z-axis
          endif
       endif

       call dcopy(lfield,tlorb1,1,solr(1:lfield,lr_),1)
       call dcopy(lfield,tlorb2,1,solz(1:lfield,lr_),1)
10  end do

    if(ncontr.gt.1) then
       ! YuP[2015/05/03] Add LCFS, if available
       ncontr_= min(ncontr,LFIELDA)
       r_surf= MAXVAL(rcontr)
       xyplotmax= max(xyplotmax,r_surf)
       do ilim=1,ncontr_
          RTAB1(ilim)=rcontr(ilim)
          RTAB2(ilim)=zcontr(ilim)
       enddo
               CALL PGSLS(2) ! 2-> dashed
               CALL PGLINE(ncontr_,RTAB1,RTAB2)
       if(machine.eq."mirror") then
                  CALL PGLINE(ncontr_,-RTAB1,RTAB2) !mirror area to the left of Z-axis
       endif
               CALL PGSLS(1) ! 1-> restore solid line
    endif
    if(nlimiter.gt.1) then
       ! YuP[2016] Add "last surface" (plasma border), if available
       nline= min(nlimiter,LFIELDA)
       r_surf= MAXVAL(rlimiter)
       xyplotmax= max(xyplotmax,r_surf)
       do ilim=1,nline
          RTAB1(ilim)=rlimiter(ilim)
          RTAB2(ilim)=zlimiter(ilim)
       enddo
               CALL PGSLW(lnwidth*2) ! bold
               CALL PGLINE(nline,RTAB1,RTAB2)
       if(machine.eq."mirror") then
                  CALL PGLINE(nline,-RTAB1,RTAB2) !mirror area to the left of Z-axis
       endif
               CALL PGSLW(lnwidth) ! restore
    endif

    !..................................................................
    !     Plot nfrplt birth points. These will all be projected onto
    !     one poloidal cross-section.
    !..................................................................

    iskip=1+ipts/nfrplt
    write(*,*)'frplteq: ipts,nfrplt,iskip',ipts,nfrplt,iskip
    if (iskip .eq. 0) then
       write(*,*) 'frplteq: iskip=0, beam missing plasma? '
    else
       do i=1,ipts,iskip
          RPG1=RPTS(I)
          if(machine.eq."mirror") then
             RPG1=XPTS(I)
             !area to the left of Z-axis is included,
             ! so the horizontal axis is X
          endif
          RPG2=ZPTS(I)
                      CALL PGPT1(RPG1,RPG2,17)
          xyplotmax= max(xyplotmax,RPG1) ! limits for plots in (X,Y)
       enddo
    endif



    !---- PLOTS in (X,Y) (top view) -------------------------------------------
    rbot=.1
    rtop=.9
          CALL PGPAGE
          CALL PGSVP(rbot,rtop,rbot,rtop)
    RBOT=-rmaxcon
    RTOP= rmaxcon
          CALL PGSWIN(rbot,rtop,rbot,rtop)
          CALL PGWNAD(-xyplotmax,xyplotmax,-xyplotmax,xyplotmax) ! limits
          CALL PGBOX('BCNST',0.,0,'BCNST',0.,0)
          CALL PGLAB('X (cms)','Y (cms)','NBI Deposition')
    ! Plot circles for the largest and smallest FP surfaces.
    nline=LFIELDA ! could be other number, but RTAB1 has LFIELDA size
    r_surf=rpcon(lrz)  ! R radius of largest FP surf, outboard
    do iline=1,nline
       tora= (iline-1)*twopi/(nline-1)
       RTAB1(iline)= r_surf*cos(tora)
       RTAB2(iline)= r_surf*sin(tora)
    enddo
          CALL PGLINE(nline,RTAB1,RTAB2)
    r_surf=rmcon(lrz)  ! R radius of largest FP surf, inboard
    do iline=1,nline
       tora= (iline-1)*twopi/(nline-1)
       RTAB1(iline)= r_surf*cos(tora)
       RTAB2(iline)= r_surf*sin(tora)
    enddo
          CALL PGLINE(nline,RTAB1,RTAB2)
    r_surf=rpcon(1)  ! R radius of smallest FP surf, outboard
    do iline=1,nline
       tora= (iline-1)*twopi/(nline-1)
       RTAB1(iline)= r_surf*cos(tora)
       RTAB2(iline)= r_surf*sin(tora)
    enddo
          CALL PGLINE(nline,RTAB1,RTAB2)
    r_surf=rmcon(1)  ! R radius of smallest FP surf, inboard
    do iline=1,nline
       tora= (iline-1)*twopi/(nline-1)
       RTAB1(iline)= r_surf*cos(tora)
       RTAB2(iline)= r_surf*sin(tora)
    enddo
          CALL PGLINE(nline,RTAB1,RTAB2)

    if(ncontr.gt.1) then
       ! YuP[2016] Add "last surface" (plasma border), if available
       nline=LFIELDA ! could be other number, but RTAB1 has LFIELDA size
       r_surf= MAXVAL(rcontr)
       do iline=1,nline
          tora= (iline-1)*twopi/(nline-1)
          RTAB1(iline)= r_surf*cos(tora)
          RTAB2(iline)= r_surf*sin(tora)
       enddo
               CALL PGSLS(2) ! 2-> dashed
               CALL PGLINE(nline,RTAB1,RTAB2)
               CALL PGSLS(1) ! 1 - solid
    endif
    if(nlimiter.gt.1) then
       ! YuP[2016] Add "last surface" (plasma border), if available
       nline=LFIELDA ! could be other number, but RTAB1 has LFIELDA size
       r_surf= MAXVAL(rlimiter)
       do iline=1,nline
          tora= (iline-1)*twopi/(nline-1)
          RTAB1(iline)= r_surf*cos(tora)
          RTAB2(iline)= r_surf*sin(tora)
       enddo
               CALL PGSLW(lnwidth*2) ! bold
               CALL PGLINE(nline,RTAB1,RTAB2)
               CALL PGSLW(lnwidth) ! restore
    endif

    !..................................................................
    !     Plot nfrplt birth points. These will all be projected onto
    !     one poloidal cross-section.
    !..................................................................

    iskip=1+ipts/nfrplt
    if (iskip .eq. 0) then
       write(*,*) 'frplteq: iskip=0, beam missing plasma? '
    else
       do i=1,ipts,iskip
          RPG1=XPTS(I)
          RPG2=YPTS(I)
                      CALL PGPT1(RPG1,RPG2,17)
       enddo
    endif



    !--------------------------------------------------------------------
    !sadness
200 continue ! to skip plots

    if (frplt.eq."plotwrit" .or. frplt.eq."write_op") then

       open(unit=19,file="freya_points.txt",status="replace")
       write(19,1000) 'Freya birth points: pnt number,x,y,Z,R,vx,vy,vz (cgs)'
       do i=1,ipts
          write(19,1001) i,xpts(i),ypts(i),zpts(i),rpts(i),vx(i),vy(i),vz(i)
       enddo
1000   format(a)
1001   format(i7,1x,7ES12.4E2)
       close(19)

    endif  ! On frplt


    !      STOP  ! TEMP
    return
  end subroutine frplteq

  subroutine micfrplt
    implicit integer (i-n), real*8 (a-h,o-z)
    save

    textt(1)="1$"
    textt(2)="2$"
    textt(3)="3$"
    textt(4)="4$"
    textt(5)="5$"
    textt(6)="6$"
    textt(7)="7$"
    textt(8)="8$"
    textt(9)="9$"
    textt(10)="10$"
    textt(11)="11$"
    textt(12)="12$"
    textt(13)="13$"
    textt(14)="14$"
    textt(15)="15$"
    textt(16)="16$"
    textt(17)="17$"
    textt(18)="18$"
    textt(19)="19$"
    textt(20)="20$"
    textt(21)="21$"
    textt(22)="22$"
    textt(23)="23$"
    textt(24)="24$"
    textt(25)="25$"
    textt(26)="26$"
    textt(27)="27$"
    textt(28)="28$"
    textt(29)="29$"
    textt(30)="30$"
    textt(31)="31$"
    textt(32)="32$"
    textt(33)="33$"
    textt(34)="34$"
    textt(35)="35$"
    textt(36)="36$"
    textt(37)="37$"
    textt(38)="38$"
    textt(39)="39$"
    textt(40)="40$"
    textt(41)="41$"
    textt(42)="42$"
    textt(43)="43$"
    textt(44)="44$"
    textt(45)="45$"
    textt(46)="46$"
    textt(47)="47$"
    textt(48)="48$"
    textt(49)="49$"
    textt(50)="50$"
    textt(51)="51$"
    textt(52)="52$"
    textt(53)="53$"
    textt(54)="54$"
    textt(55)="55$"
    textt(56)="56$"
    textt(57)="57$"
    textt(58)="58$"
    textt(59)="59$"
    textt(60)="60$"
    textt(61)="61$"
    textt(62)="62$"
    textt(63)="63$"
    textt(64)="64$"
    textt(65)="65$"
    textt(66)="66$"
    textt(67)="67$"
    textt(68)="68$"
    textt(69)="69$"
    textt(70)="70$"
    textt(71)="71$"
    textt(72)="72$"
    textt(73)="73$"
    textt(74)="74$"
    textt(75)="75$"
    textt(76)="76$"
    textt(77)="77$"
    textt(78)="78$"
    textt(79)="79$"
    textt(80)="80$"
    textt(81)="81$"
    textt(82)="82$"
    textt(83)="83$"
    textt(84)="84$"
    textt(85)="85$"
    textt(86)="86$"
    textt(87)="87$"
    textt(88)="88$"
    textt(89)="89$"
    textt(90)="90$"
    textt(91)="91$"
    textt(92)="92$"
    textt(93)="93$"
    textt(94)="94$"
    textt(95)="95$"
    textt(96)="96$"
    textt(97)="97$"
    textt(98)="98$"
    textt(99)="99$"
    textt(100)="100$"
    textt(101)="101$"
    textt(102)="102$"
    textt(103)="103$"
    textt(104)="104$"
    textt(105)="105$"
    textt(106)="106$"
    textt(107)="107$"
    textt(108)="108$"
    textt(109)="109$"
    textt(110)="110$"
    textt(111)="111$"
    textt(112)="112$"
    textt(113)="113$"
    textt(114)="114$"
    textt(115)="115$"
    textt(116)="116$"
    textt(117)="117$"
    textt(118)="118$"
    textt(119)="119$"
    textt(120)="120$"
    textt(121)="121$"
    textt(122)="122$"
    textt(123)="123$"
    textt(124)="124$"
    textt(125)="125$"
    textt(126)="126$"
    textt(127)="127$"
    textt(128)="128$"
    textt(129)="129$"
    textt(130)="130$"
    textt(131)="131$"
    textt(132)="132$"
    textt(133)="133$"
    textt(134)="134$"
    textt(135)="135$"
    textt(136)="136$"
    textt(137)="137$"
    textt(138)="138$"
    textt(139)="139$"
    textt(140)="140$"
    textt(141)="141$"
    textt(142)="142$"
    textt(143)="143$"
    textt(144)="144$"
    textt(145)="145$"
    textt(146)="146$"
    textt(147)="147$"
    textt(148)="148$"
    textt(149)="149$"
    textt(150)="150$"
    textt(151)="151$"
    textt(152)="152$"
    textt(153)="153$"
    textt(154)="154$"
    textt(155)="155$"
    textt(156)="156$"
    textt(157)="157$"
    textt(158)="158$"
    textt(159)="159$"
    textt(160)="160$"
    textt(161)="161$"
    textt(162)="162$"
    textt(163)="163$"
    textt(164)="164$"
    textt(165)="165$"
    textt(166)="166$"
    textt(167)="167$"
    textt(168)="168$"
    textt(169)="169$"
    textt(170)="170$"
    textt(171)="171$"
    textt(172)="172$"
    textt(173)="173$"
    textt(174)="174$"
    textt(175)="175$"
    textt(176)="176$"
    textt(177)="177$"
    textt(178)="178$"
    textt(179)="179$"
    textt(180)="180$"
    textt(181)="181$"
    textt(182)="182$"
    textt(183)="183$"
    textt(184)="184$"
    textt(185)="185$"
    textt(186)="186$"
    textt(187)="187$"
    textt(188)="188$"
    textt(189)="189$"
    textt(190)="190$"
    textt(191)="191$"
    textt(192)="192$"
    textt(193)="193$"
    textt(194)="194$"
    textt(195)="195$"
    textt(196)="196$"
    textt(197)="197$"
    textt(198)="198$"
    textt(199)="199$"
    textt(200)="200$"
    return
  end subroutine micfrplt

end module frplteq_mod
