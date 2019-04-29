module eflditer_mod

  !---BEGIN USE

  use bcast_mod, only : bcast
  use coefefad_mod, only : coefefad
  use r8subs_mod, only : dscal
  use restcon_mod, only : restcon
  use resthks_mod, only : resthks
  use soucrit_mod, only : soucrit

  !---END USE

!
!

contains

      subroutine eflditer ! for k=kelecg (=1) only
      use param_mod
      use comm_mod
      use r8subs_mod, only : dscal
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save

!.........................................................................
!     This routine facilitates iteration to obtain a given current density.
!     Current is calculated.
!     If this current is not acceptably close to a target current,
!     the electric field is adjusted, and the distribution is reset
!     to the beginning of the time step, so that another trial can
!     be make.
!.........................................................................


!     Return, for cases not implemented:
      if (cqlpmod.eq."enabled" ) return
      if (ngen.gt.1 .or. kelecg.ne.1) stop 'check operation of eflditer'

!-YuP: moved this control outside of this subroutine:
!     Set (stop iteration) flag and return, if beyond iteration limit:
!-YuP      if (nefiter+1.gt.nefitera) then
!-YuP         nefiter=0
!-YuP         go to 999
!-YuP      endif

!     Set diagnostic variable elecn, on first iteration:
      if (nefiter.eq.1) then
!BH171230         elecfldn(nefiter,n,l_)=elecfld(lr_)
         elecn(l_,n,nefiter)=elecfld(lr_)
         do iter=2,nefitera
            elecn(l_,n,iter)=0.
         enddo
      endif

!-YuP: moved this control outside of this subroutine:
!     Set (stop iteration) flag and return, if n.lt.(control turn on step):
!-YuP      if (n.le.noncntrl) then
!-YuP         nefiter=0
!-YuP         go to 999
!-YuP      endif

!     Calculate current (following subroutine diaggnde).
!     This is parallel current at minimum B point on a flux surface.

      currtp(lr_)=0.0
      call bcast(tam3,zero,jx)
      do 20 i=1,iy
         do 10 j=1,jx
            tam3(j)=tam3(j)+f(i,j,1,l_)*cynt2(i,l_)*coss(i,l_)
 10      continue
 20   continue

      cn=0.
      do 40 j=1,jx
         currv(j,1,l_)=tam3(j)*x(j)*gammi(j)*cint2(j)*dxi(j)
         cn=cn+currv(j,1,l_)*dx(j)
 40   continue

!     Determine jxcrit (presently determined as
!     index for amin1(3.*clight,ucrit).
      call soucrit

      curra(1,l_)=0.0
      if (jxcrit(1,lr_).ge.jx) go to 50
      do j=jxcrit(1,lr_),jx
         curra(1,l_)=curra(1,l_)+currv(j,1,l_)*dx(j)
      enddo
 50   continue

      faccur=one_*vnorm*bnumb(1)*charge
      currm(1,l_)=faccur*cn
      psifct=psiovr(lr_)/onovrp(1,lr_)

!     psifct=<(B/B_0)/R>/<1/R>, where B_0 is min mag fld of a flux surf,
!     <...> is flux surface volume average.

!     curr, is obtainted by multiplying parallel current at the
!     midplane, currm, by psifct.  This gives the average j_parallel
!     over the fixed-toroidal-angle cross-section area between
!     neighboring flux surfaces.  This approximates the same average of
!     the toroidal current, in the usual tokamak case.
!     curr is not an accurate physically meaningful quantity,
!     and should be excised from the code. (BobH, 010301).

      curr(1,lr_)=currm(1,l_)*psifct
      currtp(lr_)=currtp(lr_)+curr(1,lr_)
      curra(1,l_)=faccur*psifct*curra(1,l_)
!     Following statement for consistency with diaggnde
      call dscal(jx,faccur*psifct,currv(1:jx,1,l_),1)
!     (End of following subroutine diaggnde.)

      if (efswtchn.eq."neo_hh") then
         call restcon
         call resthks(l_,lr_,lmdpln_, &
              zreshin(lr_),zreskim(lr_),zressau1,zressau2)
         currpar(lr_)=(curra(kelec,l_)+ &
              elecfld(lr_)/300./(zreshin(lr_)*sptzr(l_)))/3.e9
      else
         currpar(lr_)=currtp(lr_)/3.e9
      endif


!           Presently, the calculated parallel current
!           above is averaged over the area cross-section
!           (at cnst toroidal angle) between a flux surface.
!           This quantity is depracated, and should be replaced
!           in the future.  Here, for method5, we adjust currpar
!           back to parallel current at the minimum B position.
!           See comments in eflditer.f, above, on psifct.

            psifct1=1.
            if(efswtch.eq."method5") then
               psifct1=psiovr(lr_)/onovrp(1,lr_)
            endif


         currerrf=(currpar(lr_)/psifct1 -currxj0(lr_))/ &
                 ( sign(one,currxj0(lr_))* &
                 max(abs(currpar(lr_)/psifct1),abs(currxj0(lr_))) )

!     Compare fractional current error with allowable.
!     In necessary, set iteration flag, obtain new efield,
!        and adjust contribution of electric field coefficients.

      if (abs(currerrf).lt.currerr) then
         nefiter_(l_)=0  !-YuP: made it a function of flux surface
      else

         elecfld_old=elecfld(lr_) !-YuP:
!        Subtract old elecfld contribution to FP equation coefficients:
         elecfld(lr_)=-elecfld_old
         call coefefad(1)

!        Find new elecfld:
         corrf=efrelax1*currerrf
         !corrf=MIN(corrf,0.95d0)   !-YuP: limited by +0.95 from above
         !corrf=MAX(-0.95d0,corrf)  !-YuP: limited by -0.95 from below
         elecfld(lr_)=elecfld_old*(1.d0-corrf)
!-YuP         call coefefad(1) !-YuP: don't add el.field here;
                               ! It will be added during next iteration
                               ! through impavnc0->coefstup->coefefad

!        Save elecfld values for diagnostic purposes:
!BH171230         elecfldn(nefiter,n,l_)=elecfld(lr_)
         elecn(l_,n,nefiter)=elecfld(lr_)
      write(*,'(a,i5,4e12.4)')'elecfld,currpar/psifct,currxj0,currerrf', &
           lr_,elecfld(lr_),currpar(lr_)/psifct1,currxj0(lr_),currerrf

      endif



      return
      end
end module eflditer_mod
