module sigfn_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      real(c_double) function sigfn(ec,kn)
      implicit integer (i-n), real(c_double) (a-h,o-z)
      save
!---------------------------------------------------------------------
!    This routine computes  reaction cross sections.
!    it is called from routine sigsetup.
!    A description of the meaning of kn=knumb is given there.
!    ec is energy in center of mass frame (keV), for fusion calc.
!    Reference for fusion reactions (kn:1==>4) is
!     H.-S. Bosch and G.M. Hale, Nucl. Fus. 32, 611 (1992),
!     which gives fusion cross-sections in millibarns.
!---------------------------------------------------------------------


      dimension a(5,4),b(4,4),erange(2,4),bg(4)
      data    a/6.927e4,7.454e8,2.050e6,5.2002e4,0.0, &
              5.7501e6,2.5226e3,4.5566e1,0.,0., &
              5.5576e4, 2.1054e2,-3.2638e-2,1.4987e-6,1.8181e-10, &
              5.3701e4,3.3027e2,-1.2706e-1,2.9327e-5,-2.5151e-9/
      data    b/6.38e1,-9.95e-1,6.981e-5,1.728e-4, &
              -3.1995e-3,-8.5530e-6,5.9014e-8,0., &
              0.,0.,0.,0., &
              0.,0.,0.,0./
      data    bg/34.3827,68.7508,31.3970,31.3970/
      data    erange/0.5,550., &
              0.3,900., &
              0.5,5000., &
              0.5,4900./


      if (kn.lt.1 .or. kn.gt.6) go to 200
      if (kn.eq.5) go to 260
      if (kn.eq.6) go to 270

      if (ec.gt.erange(1,kn) .and. ec.lt.erange(2,kn)) then

        s=a(1,kn)+ec*(a(2,kn)+ec*(a(3,kn)+ec*(a(4,kn)+ec*a(5,kn))))
        s=s/(1.+ec*(b(1,kn)+ec*(b(2,kn)+ec*(b(3,kn)+ec*b(4,kn)))))

        sigfn=s/(ec*exp(bg(kn)/sqrt(ec)))*1.e-27
        go to 200

      else

        sigfn=0.0
        go to 200

      endif

 260  continue
!     ion impact ionization
      a0=-4.203309e+1
      a1=3.557321
      a2=-1.045134
      a3=3.139238e-1
      a4=-7.454475e-2
      a5=8.459113e-3
      a6=-3.495444e-4
      eqs=ec
!990131      eqs=amax1(1.e-30,eqs)
!990131      elog=alog(eqs)
      small=1.e-30
      eqs=max(small,eqs)
      elog=log(eqs)
      siglog=(((((a6*elog+a5)*elog+a4)*elog+a3)*elog+a2)*elog+a1)*elog &
        +a0
      sigfn=exp(siglog)
      go to 200

 270  continue
!     charge exchange
      eqs=1.e+3*ec
!990131      eqs=amax1(1.,eqs)
      one=1.d0
      eqs=max(one,eqs)
!990131      enum=1.-6.731564469e-2*alog(eqs)
      enum=1.-6.731564469d-2*log(eqs)
!990131      enum=amax1(enum,0.)
      zero=0.d0
      enum=max(enum,zero)
      denom=1.+1.112e-15*eqs**3.3
      sigfn=6.937e-15*enum**2/denom

 200  return
      end
end module sigfn_mod
