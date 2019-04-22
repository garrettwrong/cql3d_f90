module diagentr_mod

!
!

contains

      subroutine diagentr(lefct,k)
      use param_mod
      use comm_mod
      use r8subs_mod
      use advnce_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!..................................................................
!     This routine computes the FSA energy transfer diagnostics
!     in Watts/cc.
!     Lefct defines the physical process under consideration.
!     Entr(k,lefct,l_) will contain the energy transfer to species "k"
!     from process lefct
!
!     lefct =
!
!     -1  means power to species "k" from background Maxwellian electrons
!     is computed.  This number (reflecting collisions) will in general
!     be negative, since the Maxwellians normally are at
!     a lower average energy than the general species.
!
!     0 means power to species "k" from all background Maxwellian ions
!     is computed.  This number (reflecting collisions) will usually
!     be negative, since the Maxwellians normally are at
!     a lower average energy than the heated general species.
!
!     1 reflects collisional power transfer to species "k" from all
!     general species (includes self-collisions - a measure of error).
!
!     2 reflects power to "k" from D.C. electric field.
!
!     3 reflects power to "k" from RF wave field.
!
!     4 entr(k,4,l_) will equal the sum of all PHYSICAL processes
!     contributing to energy transfer between and to species.
!
!     5 Ion beam source power transfer.
!
!     6 Loss orbit energy transfer (subroutine losscone).
!
!     7 Toroidal losses (subroutine losstor).
!
!     8 Power lost to system due to particles advecting over x=xmax
!     terminator.
!
!     9 Power computed from df/dt. Should equal entr(k,4,l_). In other
!     words 9 comes from the l.h.s. of the equation and 4 from the
!     r.h.s of the equation. [BH180807: Check if calc'd before or
!     after density rescaling.  Particles could be lost off the grid,
!     then rescaled back to, for example, achieve a SS.  See 13.]
!
!     10 Power gained by setting negative values of f(i,j,l_)
!     to zero after each time step.
!     This is not a physical process, but an error diagnostic.
!
!     11 Power due to synchrotron radiation (should be negative).
!
!     12 Power due to phenomenological energy losses (should
!     be negative).
!
!     13 Power due to rescaling density (BH180807: Need to add this calc).
!
!     14-
!
!BH180906:  Need to add power related to distr function scaling.
!
!


      ncentr=ncentr+1

!..................................................................
!     Zero out coefficient arrays.
!..................................................................

      call bcast(da,zero,iyjxp1)
      call bcast(db,zero,iyjxp1)
      call bcast(dc,zero,iyjxp1)

!..................................................................
!     Accordng to the value of lefct, execute a jump.
!..................................................................

      if (lefct.ne.4) entr(k,lefct,l_)=0.

      if (lefct.eq.-1) then
        do 98 j=1,jx
          do 97 i=1,iy
            da(i,j)=eal(i,j,k,1,l_)
            db(i,j)=ebl(i,j,k,1,l_)
            dc(i,j)=0.
 97       continue
 98     continue
        go to 1000  !l 227

      elseif (lefct.eq.0) then
        do 101 j=1,jx
          do 102 i=1,iy
            da(i,j)=eal(i,j,k,2,l_)
            db(i,j)=ebl(i,j,k,2,l_)
            dc(i,j)=0.
 102      continue
 101    continue
        go to 1000  !l 227

      elseif (lefct.eq.1) then
        if (colmodl .eq. 1) go to 3000
        do 103 j=1,jx
          do 1031 i=1,iy
!            da(i,j)=(1d0+0.5/(j+1)**2)*   !Kerbel Correction!
!     1              cal(i,j,k,l_)-(eal(i,j,k,1,l_)+eal(i,j,k,2,l_))
            da(i,j)=cal(i,j,k,l_)-(eal(i,j,k,1,l_)+eal(i,j,k,2,l_))
            db(i,j)=cbl(i,j,k,l_)-(ebl(i,j,k,1,l_)+ebl(i,j,k,2,l_))
            dc(i,j)=ccl(i,j,k,l_)
 1031     continue
 103    continue
        go to 1000  !l 227

      else if (lefct.eq.2) then
        if (abs(elecfld(lr_)).lt.1.e-9 .or. n.lt.nonel .or. n .gt. &
          noffel) go to 3000
        call coefefad(k)
        go to 1000

      elseif (lefct.eq.3) then
        xrf=0.
        if (vlhplse.ne."disabled" .and. (vlhmod.eq."enabled".or. &
          vlfmod.eq."enabled")) then
          if (vlhplse.eq."tauee") then
            timespn=tauee(lr_)*(vlhpon+vlhpoff)
            timeon=tauee(lr_)*vlhpon
          else
            timespn=vlhpon+vlhpoff
            timeon=vlhpon
          endif
          iperiod=timet/timespn
          tim1=iperiod*timespn
          timedf=timet-tim1
          if (timedf.le.timeon) then
            call coefrfad(k,xrf)
            go to 1000
          endif
        else if (n .ge. nonrf(k) .and. n .lt. noffrf(k)) then
          call coefrfad(k,xrf)
          go to 1000
        endif
        go to 3000

      elseif (lefct.eq.4) then
        go to 3000

      elseif (lefct.eq.5) then
        if (n .lt. nonso(k,1) .or. n .gt. noffso(k,1)) go to 3000
        call bcast(tam1,zero,jx)
        call dcopy(iyjx2,source(0:iyjx2,0,k,indxlr_),1,so,1)
        do 99 i=1,iy
          do 100 j=1,jx
            tam1(j)=tam1(j)+so(i,j)*cynt2(i,l_)*vptb(i,lr_)
 100      continue
 99     continue
        s=0.
        do 104 j=1,jx
          s=s+tam1(j)*tcsgm1(j)*cint2(j)
 104    continue
        entr(k,lefct,l_)=s*fions(k)*one_*1.602e-16*zmaxpsii(lr_)
        go to 2000

      elseif (lefct.eq.6 .or. lefct.eq.7) then
        if(lefct.eq.6) then
          call coefload(k)
          call dcopy(iyjx2,gon(0:iyjx2,0),1,temp1(0:iyjx2,0),1)
          go to 105
        endif
        call coefload(k)
        do 108 j=1,jx
          do 109 i=1,iy
            temp1(i,j)=1./taulos(i,j,indxlr_)
 109      continue
 108    continue
 105    continue
        call bcast(tam1,zero,jx)
        do 106 i=1,iy
          do 107 j=1,jx
            tam1(j)=tam1(j)-f(i,j,k,l_)*temp1(i,j) &
              *vptb(i,lr_)*cynt2(i,l_)
 107      continue
 106    continue
        s=0.
        do 110 j=1,jx
          s=s+tam1(j)*cint2(j)*tcsgm1(j)
 110    continue
        entr(k,lefct,l_)=s*one_*fions(k)*1.602e-16*zmaxpsii(lr_)
        go to 2000

      elseif (lefct.eq.8) then
        entr(k,8,l_)=sgain(4,k)/dtreff*tcsgm1(jx)*fions(k) &
          *1.602e-16*zmaxpsii(lr_)
        go to 2000

      elseif (lefct.eq.9) then
        call bcast(tam1,zero,jx)
        do 112 i=1,iy
          do 113 j=1,jx
            tam1(j)=tam1(j)+vptb(i,lr_)*(f(i,j,k,l_) &
              -f_(i,j,k,l_))*cynt2(i,l_)
 113      continue
 112    continue
        s=0.
        do 111 j=1,jx
          s=s+tam1(j)*tcsgm1(j)*cint2(j)
 111    continue
        entr(k,9,l_)=s*one_*fions(k)*1.602e-16*zmaxpsii(lr_)/dtreff
        go to 3000

      elseif (lefct.eq.10) then
        entr(k,10,l_)=engain(k)*1.602e-16*fions(k)/dtreff*zmaxpsii(lr_)
        go to 2000

      elseif (lefct.eq.11) then
        if (k.ne.kelecg .or. syncrad.eq."disabled") go to 3000
        call coefsyad(k)
        go to 1000

      elseif (lefct.eq.12) then
        if(tauegy(k,lr_).le.0.)  go to 3000
        call coefegad(k)

      elseif (lefct.ge.13) then
        call diagwrng(5)
      endif

!..................................................................
!     Compute the power by evaluating the r.h.s. of the pde.
!..................................................................

 1000 continue
      call bcast(tam1,zero,jx)
      call bcast(tam2,zero,jx)
!..................................................................
!     Redefine the coefficients at the j+1/2 mesh points
!..................................................................

      call coefmidv(da,1)
      call coefmidv(db,2)
      call coefmidv(dc,3)

!..................................................................
!     If splitting scheme was utilized..
!..................................................................

      if (implct .eq. "disabled") then
        call dcopy(iyjx2,f_(0:iyjx2,0,k,l_),1,temp1(0:iyjx2,0),1)
      endif

!..................................................................
!     Collisional powers
!..................................................................

      !YuP: this (i,j)-loop is the main "drain" for cpu time
      do 1001 i=1,iy
        !09-09-2015 checked: almost no change when skipping i=1 and i=iy
        s_cynt2=cynt2(i,l_)
        if (cqlpmod .ne. "enabled") then
          if (i .ge. itl .and. i .le. iyh) s_cynt2=2.*cynt2(i,l_)
          ! Factor of 2 above because of skipping half trapped region?
          if (i .gt. iyh .and. i .le. itu) go to 1001 !next i
        endif
        if (implct .eq. "disabled") then
          gfu_o=gfu(i,0,k)
          do 1002 j=1,jx-1
            gfu_n=gfu(i,j,k)
            tam1(j)=tam1(j)+(gfu_n-gfu_o)*s_cynt2 ! Difference of fluxes
            tam2(j)=tam2(j)+(gfu_n+gfu_o)*0.5*s_cynt2 ! Median flux
            gfu_o=gfu_n
 1002     continue
        else  !i.e., implct case
          gfi_o=gfi(i,0,k)
          do 1003 j=1,jx-1
            gfi_n=gfi(i,j,k)
            tam1(j)=tam1(j)+(gfi_n-gfi_o)*s_cynt2 ! Difference of fluxes
            tam2(j)=tam2(j)+(gfi_n+gfi_o)*0.5*s_cynt2 ! Median flux
            gfi_o=gfi_n !2x-faster way of calc. SUM[gfi(i,j)-gfi(i,j-1)]
 1003     continue
        endif
 1001 continue ! i
      tott=0.
      do 1004 j=1,jx-1
        tott=tott+tam1(j)*tcsgm1(j)  !tcsgm1=2*(clight/vnorm)**2*(gamma-1)
 1004 continue

!..................................................................
!     Above collisional calc uses power=integral du**3*energy*df/dt.
!     Would be interesting to compare with calc after u-integration
!     by parts, which sums up gfu*(u/gamma).  [BH180803].
!..................................................................



!..................................................................
!     Compute the local (in v) RF power deposition and cumulative sum
!..................................................................

      const= 1.602e-16*zmaxpsii(lr_)*fions(k)*one_
      const_m= 1.e-7*vnorm**2*fmass(k)*zmaxpsii(lr_)*one_
      if (lefct.eq.3) then
        !call bcast(pwrrf(1,k,l_),zero,jx)  !-YuP: not really needed
        !call bcast(pwrrfs(1,k,l_),zero,jx) !-YuP: not really needed
        do 1005 j=1,jx
          !pwrrf(j,k,l_)=tam1(j)*tcsgm1(j)*const*dxi(j)
          ![YuP 09-09-2015] The above expression is based on
          !difference of fluxes, like (F_j+0.5 - F_j-0.5)/dx .
          !Such form is prone to numerical errors because effectively
          !it involves the 2nd deriv. of distr.func.
          !Better - integrate analytically by parts,
          !then - the power is expressed through the flux itself
          !(not a derivative of flux). This is how it is done for powurf
          !(or powrf array), so why not the same approach here, too?
          ! Here is the alternative expression, based on flux
          !(instead of deriv. of flux):
          pwrrf(j,k,l_)= -tam2(j)*const_m*x(j)*gammi(j)
          !With this version, The profile of pwrrf is much better -
          !no negative part, and no "spikes" in plot versus x(j).
          !The total integrated powers - same as with original version.
          !BH180808: Actually, there is no particular problem with the
          ! negative part of pwrrf(j,,), or s"spikes".  After the
          ! integration by parts, a different understanding of pwrrf(j,,)
          ! may emerge.  I always thought the neg part pwrrf(j,,)
          ! represented the effect of digging out of the distribution
          ! function by the RF diffusion, and QL diffusing these
          ! particles to higher velocity.
 1005   continue ! j
        pwrrfs(1,k,l_)=dx(1)*pwrrf(1,k,l_)
        do j=2,jx
          pwrrfs(j,k,l_)=pwrrfs(j-1,k,l_)+dx(j)*pwrrf(j,k,l_)
        enddo
!BH100715        write(*,*)'pwrrfs(jx,k,l_),l_:',l_,pwrrfs(jx,k,l_)
      endif ! lefct.eq.3

      entr(k,lefct,l_)=tott*const
 2000 entr(k,4,l_)=entr(k,4,l_)+entr(k,lefct,l_)
 3000 continue
      return
      end
!
!
!=====================================================================
!
!
      subroutine diagentr_vol
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!     Integrate the FSA power densities over volume.

      if (cqlpmod.eq."enabled") return

!BH050616:  There is a code misfunction below, wherein
!           entrintr(k,lefct) appears not to be zeroed
!           out from time-step to time step in the
!           following.  However, when write(*,*) statements
!           are enabled, the problem disappears.  Have not
!           been able to track down problem.  Will try with
!           updated lf95 compiler at GA.
!           BH added the if(ll.eq.1)-statement below, and
!           code seems to be executing OK....  But needs more work!

!BH120807:  corrected previous use of do lefect=-1,15, since conflict w
!BH120807:  subrtine arg.  Maybe this was cause of BH050616 prblm below?
      do llefct=-1,15
         do k=1,ngen
            entrintr(k,llefct)=0.0
!            write(*,*)'diagentr: entrintr(k,llefct)',entrintr(k,llefct)
            do ll=1,lrz
!BH050616:  Added if statement since above zeroing of entrintr
!BH050616:  is not working?
               if(ll.eq.1) entrintr(k,llefct)=0.0
               lr_=lrindx(ll)
               entrintr(k,llefct)=entrintr(k,llefct)+ &
                    entr(k,llefct,ll)*dvol(lr_)
            enddo
!            write(*,*)'diagentr: entrintr(k,llefct)',entrintr(k,llefct)
         enddo
      enddo
!      write(*,*)'diagentr: entrintr',(entrintr(1,llefct),llefct=-1,15)
!      write(*,*)'diagentr: entr(1,5,*)',(entr(1,5,ll),ll=1,lrz)

      return
      end


!
!
!=====================================================================
!
!
      real*8 function gfi(i,j,k)
      use param_mod
      use comm_mod
      use advnce_mod, only : cl, fpj0, fpj, fpjp
!..................................................................
!     Express the velocity flux at (i,j+1/2)
!..................................................................
      implicit integer (i-n), real*8 (a-h,o-z)

      if((i.ne.itl .and. i.ne.itu) .or. symtrap.ne."enabled") then
       gfii= dc(i,j)*0.5*dyi(i,l_)*(fpjp(i,j,k)-fpj0(i-1,j,k))
       else
       gfii= &
        +cl(itl-1,j)*eyp5(itl-1,l_)*(fpjp(itl-1,j,k)-fpj(itl-1,j,k)) &
        +2.*cl(itl+1,j)*eyp5(itl,l_)*(fpj(itl+1,j,k)-fpj(itl,j,k)) &
        +cl(itu+1,j)*eyp5(itu,l_)*(fpj(itu+1,j,k)-fpj0(itu,j,k))
      endif

      gfi=da(i,j)*fpj(i,j,k) &
         +db(i,j)*exp5(j)*(f(i,j+1,k,l_)-f(i,j,k,l_)) +gfii

      end function gfi

!
!
!=====================================================================
!
!
      real*8 function gfu(i,j,k)
      use param_mod
      use comm_mod
      use advnce_mod, only : cdf, f1j, f2j
!..................................................................
!     Define the flux related quantities G ( H not needed for
!     the power calculation) used in the r.h.s. of the FP eqn.
!..................................................................
      implicit integer (i-n), real*8 (a-h,o-z)

      if(i .ne. itl  .and. i .ne. itu) then
        gfuu=dc(i,j)*(f1j(i+1,j)-f1j(i-1,j))*0.5*dyi(i,l_)
      else
        gfuu=cdf(j)
      endif


      gfu=da(i,j)*f2j(i,j) &
         +db(i,j)*(temp2(i,j+1)-temp2(i,j))*exp5(j) +gfuu

      end function gfu
!
!=====================================================================
!
end module diagentr_mod
