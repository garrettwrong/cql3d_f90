module lossegy_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine lossegy
      use param_mod
      use cqlcomm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!
!mnt  this routine computes fp coefficients corresponding to the
!mnt  bounce average of the phenomonological energy loss term,
!mnt  or a Bremsstrahlung energy loss term:

!     Phenomonological loss:
!     df/dt = u**-2*d/du [1./tau(u,theta,lr_)*u**3*f/2.]
!     where   tau is energy loss time specified by
!     tau=rfact*tauegy*gamma**gamegy, u.le.vte
!     rfact*tauegy*gamma**gamegy/
!     (abs(u_par,lr_)/vte)**paregy*(u_per/vte)**pe
!     *(u/vte)**pegy),  u.gt.vte
!     rfact is a radially dependent modifier, as below.
!
!     Bremsstrahlung loss due to electron collisions with (only)
!     the Maxwellian ion species:
!     df/dt = u**-2*d/du [beta*u**3*f],
!     where beta=reden(kion,lr_)*sigmarad*clight**2/(u/gamma),
!           sigmarad=as given in Physics Vade Mecum, H.L.Anderson, Ed.,
!                    AIP, Sect. 16.07.F.
!     Sigmarad may be enhanced beyond energy brfacgm3*16.43/Z**(1/6),
!     using brfac,brfac1 (see cqlinput_help and coding below).
!     This is for purposes of stopping flow of runaway electrons
!     off the grid.



      data gam1/2.30097589089281/,   gam2/16.4312912649856/

      do 200 k=1,ngen

         call bcast(egylosa(0:(iy+1),0:(jx+1),k,indxlr_),zero,(iy+2)*(jx+2))
         !Note: allocate(egylosa(0:iy+1,0:jx+1,ngen,setup0%lrz)

!     Phenomenological:

         if(tauegy(k,lr_).le.0.)  go to 130
         coefnt=1./(2.*tauegy(k,lr_))
!     rfact modifier multiplies tau by (1.-rovera(lr_),lr_):
         if(regy(k).eq."enabled") coefnt=coefnt/(1.-rovera(lr_))
!
!..................................................................
!     Note: the do loop below uses vth(),
!     vth is the thermal velocity =sqrt(T/m) (at t=0 defined in ainpla).
!     But, T==temp(k,lr) can be changed in profiles.f,
!     in case of iprote (or iproti) equal to "prbola-t" or "spline-t"
!..................................................................
         do 100  j=1,jx
            vel=x(j)*vnorm
            coefnt1=coefnt*xsq(j)*x(j)/gamma(j)**gamegy(k)
            if(vel.le.vth(k,lr_))  then
               do 110  i=1,iy
                  egylosa(i,j,k,indxlr_)=vptb(i,lr_)*coefnt1
 110           continue
            else
               do 120  i=1,iy
                  egylosa(i,j,k,indxlr_)= &
                       vptb(i,lr_)*coefnt1*(vel/vth(k,lr_))** &
                       (paregy(k)+peregy(k)+pegy(k))*sincosba(i,k,lr_)
 120           continue
            endif
 100     continue
 130     continue

!     Bremsstrahlung:
         if (k.eq.kelecg .and. bremsrad.eq."enabled") then
            if(brfacgm3.lt.1.) stop 'Check brfacgm3, see cqlinput_help'
            do 150 kk=1,nionm
               kkk=kionm(kk)
               c1=fmass(kelecg)*clite2
               c2=(charge**2/c1)**2/137.*bnumb(kkk)**2
               c3=1./bnumb(kkk)**(1./3.)
               c5=sqrt(c3)
               c4=183.*c3
               c3=137.*c3
               gam3=c5*gam2
               gam4=brfacgm3*gam3
               sigmarad=0.
               do 160 j=1,jx
                  if (gamm1(j).lt.gam1) then
                     sigmarad=16./3.*c2
                  elseif (gamm1(j).gt.gam1.and. &
                         gamm1(j).le.gam3) then
!990131                     sigmarad=8.*c2*(alog(gamm1(j))-1./6.)
                     sigmarad=8.*c2*(log(gamm1(j))-1./6.)
!990131                     sigmarad=8.*c2*(alog(gamm1(j))-1./6.)
                     sigmarad=8.*c2*(log(gamm1(j))-1./6.)
                  elseif (gamm1(j).gt.gam3.and. &
                        gamm1(j).le.gam4) then
!990131                     sigmarad=4.*c2*(alog(c4)+1./18.)
                     sigmarad=4.*c2*(log(c4)+1./18.)
                  elseif (gamm1(j).gt.gam4) then
                     gam5=gamm1(j)-gam4
!990131                     sigmarad=4.*c2*(alog(c4)+1./18.)*(1.+brfac*gam5
                     sigmarad=4.*c2*(log(c4)+1./18.)*(1.+brfac*gam5 &
                          **brfac1)
                  endif
                  betau3=reden(kkk,lr_)*sigmarad*clite2*gamma(j)* &
                       xsq(j)/vnorm
                  do 170 i=1,iy
                     egylosa(i,j,k,indxlr_)=egylosa(i,j,k,indxlr_)+ &
                          vptb(i,lr_)*betau3
 170              continue
 160           continue
 150        continue
         endif

 200  continue

      return
      end subroutine lossegy


end module lossegy_mod
