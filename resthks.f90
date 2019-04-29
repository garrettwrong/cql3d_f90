module resthks_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine resthks(kl_,klr_,klmdpln,preshin,preskim,pressau1, &
        pressau2)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!     ------------------------------------------------------------------
!
!     Compute resistivity ratios (to Spitzer Z=1 value), xxxBH,010807
!     BH,010807: Spitzer with Zeff.ne.1 to be used with HH...,below.
!                Olivier evidently adjusted the other ratios also.
!
!     including collisionality and Zeff, according to
!     Hinton and Hazeltine (1976), Kim-Hirshman (1988) and Sauter (1993)
!
!     preshin=Hinton and Hazeltine ratio
!     preskim=Kim (Hirshman-Sigmar) ratio
!     pressau1=Sauter1
!     pressau1=Sauter2
!
!     Assume restcon called just before, i.e. xconn corresponds to kl_
!
      save
!.......................................................................
!     N(Zeff) from Spitzer such that N(1,2,4,16,inf)=
!                                               1.,0.85,0.74,0.63,0.58
      zzef=zeff(klr_)
!     gives 1.00045 for Z=1 and less than 1% error up to infinity
      znzeff=0.58+0.74/(0.76+zzef)
!
!l    1. Hinton-Hazeltine (Rev.Mod.Phys. 48, 239 (1976))
!        Also, see ONETWO manual, GA-A16178.
!
      za33=4./3.*(1.-1./zzef)*0.47 + 1./3.*(4./zzef-1.)*0.68
      zb33=4./3.*(1.-1./zzef)*0.20 + 1./3.*(4./zzef-1.)*0.32
      zc33=4./3.*(1.-1./zzef)*0.51 + 1./3.*(4./zzef-1.)*0.66
      zk033=1.46+(1.83-1.46)/zzef
      zfcol=zk033/(1.+za33*sqrt(starnue(klmdpln))+zb33*starnue(klmdpln)) &
        /(1.+zc33*starnue(klmdpln)*eps(klr_)**1.5)
      zftraph=zfcol*sqrt(eps(klr_)) - (zfcol-1.)*eps(klr_)

!%OS  preshin=znzeff/(1.-zftraph)
      preshin=1.0/(1.-zftraph)
!
!l    2. Kim from thesis, which is Hirshman-Sigmar formula
!
!     note: znzeffk is equivalent to: 0.60 + 0.59/(0.43+Zeff)
      zg=(1.-xconn)/xconn
      zftrapk=(zzef**2+1.41*zzef+(2.*zzef**2+2.66*zzef+0.75)*zg+ &
        (zzef**2+1.24*zzef+0.35)*zg**2)/ &
        (3.25*zzef**2+1.41*zzef+(3.25*zzef**2+1.39*zzef)*zg)
      preskim=zftrapk*(3.25*zzef**2+1.41*zzef)/(zzef**2+1.41*zzef)

!
!l    3. Sauter formula
!        O. Sauter, C. Angioni and Y.R. Lin-Liu, Phys. of Plasmas 6,
!        1834 (1999).
!
!%OS  pressau1=znzeff*(1.+1.24*zg)
      pressau1=1.0*(1.+1.24*zg)
      pressau2=(1.+1.24*zfcol/zk033*zg)

      return
      end
end module resthks_mod
