module exsweep_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use coefmidt_mod, only : coefmidt
  use coefmidv_mod, only : coefmidv
  use coefstup_mod, only : coefstup
  use coefwti_mod, only : coefwti
  use coefwtj_mod, only : coefwtj
  use diagxswt_mod, only : diagxswt
  use diagxswx_mod, only : diagxswx
  use exsweept_mod, only : exsweept
  use exsweepx_mod, only : exsweepx
  use r8subs_mod, only : dcopy

  !---END USE

!
!

contains

      subroutine exsweep
      use param_mod
      use cqlcomm_mod
      use advnce_mod !here: in exsweep. To get alpx(), betx(), etc.
      use r8subs_mod, only : dcopy
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     Time advancement routine - uses splitting scheme
!     (called if implct.ne.enabled)
!..................................................................

      save

!.......................................................................

      navnc=navnc+1

!..................................................................
!     Copy current distribution into f_.
!..................................................................

      call dcopy(iyjx2*ngen,f(0:iy+1,0:jx+1,1:ngen,l_),1, &
                           f_(0:iy+1,0:jx+1,1:ngen,l_),1)

!..................................................................
!     loop over all time advanced species..
!..................................................................

      do 600 k=1,ngen

!..................................................................
!     normalized time step..
!..................................................................

        rbgn=1./dtr

!..................................................................
!     The routine coefstup(k) does the following:
!     Preset all coefficients to zero. da, db, and dc are the three
!     coefficients asscoiated with the velcoity (v) flux while
!     dd, de, df are associated with the theta flux. The contribution
!     from the Krook operator (cah(i,j)*f(i,j,l_)) is stored in cah.
!     Note cah does not contribute to da through df. These six
!     coefficients incorporate only wave-particle, collisional and
!     d.c. electric field effects. The ion source is represented
!     through so(i,j).
!     The routine coefstup(k) adds into the above arrays all the relevan
!     physics.
!..................................................................

        call coefstup(k)

!..................................................................
!     The coefficients of the equation are currently defined on the
!     same mesh as the distribution function f. The fluxes are best
!     defined (from the point of view of differencing and particle
!     conservation) on mid meshpoints. We undertake here to
!     interpolate the coefficients as needed to either (i,j+1/2)
!     (velocity flux) or to (i+1/2,j) (theta flux).
!     Finally to enforce boundary conditions (zero flux in general
!     except at the pass/trapped boundary) certain coefficients
!     are zeroed out or suitably averaged at specific mesh points.
!     The numbers 1,2,3 appearing in the calsetup0%ls below signify
!     which coefficient is being treated.
!..................................................................

!..................................................................
!     first the velocity flux..
!..................................................................

        call coefmidv(da,1)
        call coefmidv(db,2)
        call coefmidv(dc,3)

!..................................................................
!     the theta flux..
!..................................................................

        call coefmidt(dd,1)
        call coefmidt(de,2)
        call coefmidt(df,3)

!..................................................................
!     The differencing is done ala Chang and Cooper as generalized
!     by Karney to 2-D. This requires defining a couple of arrays
!     di(i,j,k,l_) and dj(i,j,k,l_) which are used to define f at mid mesh
!     points: to wit f(i,j+1/2,k,l_)=f(i,j+1,k,l_)*(1.-dj(i,j,k,l_))+
!     f(i,j,k,l_)*dj(i,j,k,l_). Similarly for theta.
!     The routine coefwtj and coefwti are used to
!     determine the dj and di arrays.
!..................................................................

        call coefwtj(k)
        call coefwti(k)

!..................................................................
!     copy the distribution function into temp1
!..................................................................

        call dcopy(iyjx2,f(0:iy+1,0:jx+1,k,l_),1,temp1(0:iy+1,0:jx+1),1)

!..................................................................
!     Initialize at x=0 for the forward recursion of the velocity
!     sweep. This boundary condition is a conservative 0 flux condition.
!     The quantities egg and fgg defined below in the forward recursion
!     will be used in the backward recursion to compute f (at n+1/2).
!     f(i,j,l_)=egg(i,j)*f(i,j+1,l_)+fgg(i,j)
!..................................................................

!dir$ ivdep
        do 20 i=1,iy

!..................................................................
!     da(i,0)=db(i,0)=dc(i,0)=0 forces G(i,0)=0. Thus the equation at
!     x=0 is effectively
!
!     df/dt = G(i,1)*integration coefficient
!
!     Recall the quantity G(i,j) is really defined at G(i,j+1/2).
!..................................................................

          fgg(i,1)=delx(i,1,k,l_)/betx(i,1,k,l_)
          egg(i,1)=alpx(i,1,k,l_)/betx(i,1,k,l_)
 20     continue

!..................................................................
!     Complete the forward recursion for egg and fgg
!..................................................................

        do 50 j=2,jxm1
!dir$ ivdep
          do 40 i=1,iy
            egg(i,j)=alpx(i,j,k,l_)/(betx(i,j,k,l_)-gamx(i,j,k,l_)*egg(i,j-1))
            fgg(i,j)=(delx(i,j,k,l_)+gamx(i,j,k,l_)*fgg(i,j-1))/ &
                     (betx(i,j,k,l_)-gamx(i,j,k,l_)*egg(i,j-1))
 40       continue
 50     continue

!..................................................................
!     A boundary equation is needed at x=xmax (j=jx). This is
!     effectively a free fall (hyperbolic) condition if characteristics
!     are leaving the system and zero flux with f(i,jx,l_) approximately
!     zero if characteristics are entering the system.
!     Note that we have arbitrarily set db(i,jx)=db(i,jx-1)=dc(i,jx)=
!     dc(i,jx-1)=0 for all i for purposes of time advancement. If
!     advection dominates (da dominates) setting diffusion equal to
!     zero does not pose a problem. If diffusion dominates then xmax
!     has to have been chosen high enough that f is close to zero. Thus
!     setting db and dc equal to zero here should not disturb the
!     solution.
!
!     Determine temc1 and temc2 such that
!
!     f(i,jx,l_)=temc1(i)*f(i,jx-1,l_)+temc2(i)
!..................................................................

!dir$ ivdep
        do 60 i=1,iy
          temc1(i)=gamx(i,jx,k,l_)/betx(i,jx,k,l_)
          temc2(i)=delx(i,jx,k,l_)/betx(i,jx,k,l_)
 60     continue

!..................................................................
!     Solve for f(i,jx,l_)=temp2 using the coefficients for the equation
!     developed in do loops 60 and 40.
!..................................................................

        do 70 i=1,iy
          temp2(i,jx)=(temc1(i)*fgg(i,jx-1)+temc2(i)) &
            /(1.-temc1(i)*egg(i,jx-1))
 70     continue

!..................................................................
!     Backsolve to determine solution at end of velocity split.
!..................................................................

        do 90 jj=1,jxm1
          j=jx-jj
          do 80 i=1,iy
            temp2(i,j)=egg(i,j)*temp2(i,j+1)+fgg(i,j)
 80       continue
 90     continue

!..................................................................
!     if desired call routine to check accuracy of differencing
!..................................................................

        if (iactst.ne. "disabled") call exsweepx(k)

!..................................................................
!     compute density diagnostics and if desired set negative values
!     of distribution function to 0
!..................................................................

        call diagxswx(k)
!..................................................................
!     store distribution function after velocity split in fxsp
!..................................................................

        call dcopy(iyjx2,temp1(0:iy+1,0:jx+1),1,fxsp(0:iy+1,0:jx+1,k,l_),1)

!..................................................................
!     velocity split done - commence with theta split
!
!
!     There will be three (forward) sweeps all of which tie into each
!     other at the pass/trapped boundary. Begin the forward sweeps at
!     y=0, pi/2 and pi. They terminate at  y(itl,l_), y(itl,l_) or y(itu
!     (the lower or upper p/t boundary. All sweeps utilize 0 flux
!     boundary conditions at the outside edges.
!..................................................................

!dir$ ivdep
        do 110 j=2,jx
          egg(1,j)=  alpy(1,  j,k,l_)/bety(1,  j,k,l_)
          egg(iyh,j)=gamy(iyh,j,k,l_)/bety(iyh,j,k,l_)
          egg(iy,j)= gamy(iy, j,k,l_)/bety(iy, j,k,l_)
          fgg(1,j)=  dely(1,  j,k,l_)/bety(1,  j,k,l_)
          fgg(iyh,j)=dely(iyh,j,k,l_)/bety(iyh,j,k,l_)
          fgg(iy,j)= dely(iy, j,k,l_)/bety(iy, j,k,l_)
 110    continue

!..................................................................
!     complete the first sweep over the three regions of the theta mesh
!..................................................................

        do 130 i=2,itl-1
!dir$ ivdep
          do 120 j=2,jx
            egg(i,j)=alpy(i,j,k,l_)/(bety(i,j,k,l_)-gamy(i,j,k,l_)*egg(i-1,j))
            fgg(i,j)=(dely(i,j,k,l_)+gamy(i,j,k,l_)*fgg(i-1,j))/(bety(i,j,k,l_) &
                    -gamy(i,j,k,l_)*egg(i-1,j))
 120      continue
 130    continue
        do 150 i=iy-1,iy+2-itl,-1
          do 140 j=2,jx
            egg(i,j)=gamy(i,j,k,l_)/(bety(i,j,k,l_)-alpy(i,j,k,l_)*egg(i+1,j))
            fgg(i,j)=(dely(i,j,k,l_)+alpy(i,j,k,l_)*fgg(i+1,j))/(bety(i,j,k,l_) &
                    -alpy(i,j,k,l_)*egg(i+1,j))
 140      continue
 150    continue
        do 170 i=iyh-1,itl+1,-1
          do 180 j=2,jx
            egg(i,j)=gamy(i,j,k,l_)/(bety(i,j,k,l_)-alpy(i,j,k,l_)*egg(i+1,j))
            fgg(i,j)=(dely(i,j,k,l_)+alpy(i,j,k,l_)*fgg(i+1,j))/(bety(i,j,k,l_) &
                    -alpy(i,j,k,l_)*egg(i+1,j))
 180      continue
 170    continue

!..................................................................
!     At the pass/trapped boundary conservation conditions hold. The
!     r.h.s. of the Fokker-Planck equation will be proportional to:
!
!     -H(itl-1/2) + 2.*H(itl+1/2) + H(itu+1/2).
!
!     This translates to an equation of the form:
!
!     -alpha*f(itl+1,l_) +beta*f(itl,l_) -gamma*f(itl+1,l_) -alpa*f(itu
!
!     = delta(itl)
!
!     and these coefficients are near neighbors to alpy, bety, etc
!     used at the other points of the domain. Compute alpha - delta
!     and store them in tam1-tam5
!..................................................................

        do 185 j=2,jx
          tam1(j)=alpy(itl,j,k,l_)
          tam2(j)=.5*(bety(itl,j,k,l_)+bety(itu,j,k,l_))
          tam3(j)=.5*gamy(itl,j,k,l_)
          tam4(j)=.5*alpy(itu,j,k,l_)
          tam5(j)=.5*(dely(itl,j,k,l_)+dely(itu,j,k,l_))
 185    continue

!..................................................................
!     The equation whose coefficients are defined in do loop 180
!     involves f at the points itl+1, itl-1, itu+1 as well as itl.
!     The recursive equations developed in do loops 120, 140, 160
!     can be employed to recast this equation so that it has only
!     one unknown: f(itl,j,l_). Solve for f(itl,j,l_)=f(itu,j,l_).
!..................................................................

        do 190 j=2,jx
          tam6(j)=tam2(j)-tam1(j)*egg(itl+1,j)-tam3(j)*egg(itl-1,j)- &
            tam4(j)*egg(iy+2-itl,j)
          tam7(j)=tam5(j)+tam1(j)*fgg(itl+1,j)+tam3(j)*fgg(itl-1,j)+ &
            tam4(j)*fgg(iy+2-itl,j)
          temp2(itl,j)=tam7(j)/tam6(j)
          temp2(itu,j)=temp2(itl,j)
 190    continue

!..................................................................
!     backsolve over three regions to obtain distribution
!..................................................................

        do 210 i=itl-1,1,-1
!dir$ ivdep
          do 200 j=2,jx
            temp2(i,j)=egg(i,j)*temp2(i+1,j)+fgg(i,j)
 200      continue
 210    continue
        do 230 i=itl+1,iyh
!dir$ ivdep
          do 220 j=2,jx
            temp2(i,j)=egg(i,j)*temp2(i-1,j)+fgg(i,j)
 220      continue
 230    continue
        do 250 i=iy+2-itl,iy
!dir$ ivdep
          do 240 j=2,jx
            temp2(i,j)=egg(i,j)*temp2(i-1,j)+fgg(i,j)
 240      continue
 250    continue

!..................................................................
!     Add the contribution of krook operator + source at v=0.
!     At x=0 f is clearly isotropic and will be averaged
!     over theta below.
!..................................................................

        do 260 i=1,iy
          temc1(i)=.5*so(i,1)+temp1(i,1)*rbgn
          temc2(i)=rbgn-.5*cah(i,1)
          temp2(i,1)=temc1(i)/temc2(i)
 260    continue

!..................................................................
!     update f in upper half of trapping region using symmetry.
!..................................................................

        do 270 j=1,jx
!dir$ ivdep
          do 271 i=itl+1,iyh
            temp2(iy+1-i,j)=temp2(i,j)
 271      continue
 270    continue

!..................................................................
!     difference solution - compare with right side.
!..................................................................

        if (iactst.ne."disabled") call exsweept(k)

!..................................................................
!     compute density diagnostics and set negative values of f to 0.
!..................................................................

        call diagxswt(k)

!..................................................................
!     redefine f at v=0 so it is unique
!..................................................................

        s=0.
        t=0.
        do 280 i=1,iy
          s=s+vptb(i,lr_)*cynt2(i,l_)
          t=t+vptb(i,lr_)*cynt2(i,l_)*temp2(i,1)
 280    continue
        do 512 i=1,iy
          temp1(i,1)=t/s
 512    continue

!..................................................................
!     store new distribution in f.
!..................................................................

        call dcopy(iyjx2,temp1(0:iy+1,0:jx+1),1,f(0:iy+1,0:jx+1,k,l_),1)
 600  continue
      return
      end subroutine exsweep


end module exsweep_mod
