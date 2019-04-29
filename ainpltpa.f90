module ainpltpa_mod

  !---BEGIN USE
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_double

  !---END USE

!
!

contains

      subroutine ainpltpa
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

      REAL RILIN

!..................................................................
!     This routine plots out the parameters set in the main module
!..................................................................

!MPIINSERT_INCLUDE


!.....................................................................
!     Plot out parameters (set in source code before compilation)
!.....................................................................


!MPIINSERT_IF_RANK_NE_0_RETURN
 ! make plots on mpirank.eq.0 only

      if (noplots.eq."enabled1") return

      CALL PGPAGE
      RILIN=0.
      CALL PGMTXT('T',-RILIN,0.,0.,"PARAMETER VALUES")

      write(t_,12998) version
12998 format("====> version =  ",a50)
      RILIN=RILIN+2.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,12999) precursr
12999 format("====> precursr = ",a50)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13000)
13000 format("ngena is the max. # of general (time advanced) species")
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13001) ngena
13001 format("====> ngena = ",i5)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13002)
13002 format("nmaxa is the max. # of background Maxwellian species")
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13003) nmaxa
13003 format("====> nmaxa = ",i5)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!      write(t_,13006)
13006 format("jx is the number of velocity mesh points")
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!      write(t_,13007) jx
13007 format("====> jx = ",i5)
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!      write(t_,13008)
13008 format("iy is the number of theta mesh points")
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!      write(t_,13009) iy
13009 format("====> iy = ",i5)
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13010)
13010 format("lza is the maximum number of z mesh points")
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13011) lza
13011 format("====> lza = ",i5)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!      write(t_,13012)
13012 format("mx is the highest order legendre polynomial employable")
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!      write(t_,13013) mx
13013 format("====> mx = ",i5)
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      if (lrzmax .gt. 1) then

         write(t_,13014)
13014    format("lrza is the maximum number of radial flux surfaces")
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,t_)

         write(t_,13015) lrza
13015    format("====> lrza = ",i5)
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,t_)

         write(t_,13016)
13016    format("analytic source routine parameters")
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,t_)

         write(t_,13017)
13017    format( "nsoa is the number of sources per species.")
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,t_)

         write(t_,13018) nsoa
13018    format("====> nsoa = ",i5)
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,t_)


      endif

!      write(t_,13019)
13019 format("noncha is the number of elements in arrays carrying time")
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!      write(t_,13020)
13020 format("plot information: last noncha time steps are remembered.")
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!      write(t_,13021) noncha
13021 format("====> noncha = ",i5)
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13022)
13022 format("nplota is max number of plot times for 2d and 3d plots.")
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13023) nplota
13023 format("====> nplota = ",i5)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13024)
13024 format("nbctimea is max number of points in arrays giving time")
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13025)
13025 format("dependent profile information.")
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13026) nbctimea
13026 format("====> nbctimea = ",i5)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13027)
13027 format("ndtr1a is maximum number of time step intervals dtr1().")
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13028) ndtr1a
13028 format("====> ndtr1a = ",i5)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13029)
13029 format("nefitera is the max number of iterations permitted for")
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13030)
13030 format("electric field per time step (to obtain target current).")
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      write(t_,13031) nefitera
13031 format("====> nefitera = ",i5)
      RILIN=RILIN+1.
      CALL PGMTXT('T',-RILIN,0.,0.,t_)


!      write(t_,13032)
13032 format("f(vpar,vperp) for reduced distn F(vpar)")
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

!      write(t_,13033) ipxy,jpxy
13033 format("====> ipxy = ",i5,"      jpxy = ",i5)
      RILIN=RILIN+1.
!      CALL PGMTXT('T',-RILIN,0.,0.,t_)

      if (frmodp.eq."enabled") then
!         write(t_,13034)
13034    format("npart is the number of ions launched by NFREYA")
         RILIN=RILIN+1.
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)

!         write(t_,13035) npart
13035    format("====> npart = ",i5)
         RILIN=RILIN+1.
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)
      endif

      if (urfmod.ne."disabled") then

!         write(t_,13036)
13036    format("nrayn is the number of rays.")
!         RILIN=RILIN+1.
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)

!         write(t_,13037) nrayn
13037    format("====> nrayn = ",i5)
!         RILIN=RILIN+1.
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)

!         write(t_,13038)
13038    format("nrayelts is the max number of ray elements per ray")
!         RILIN=RILIN+1.
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)

!         write(t_,13039) nrayelts
13039    format("====> nrayelts = ",i5)
!         RILIN=RILIN+1.
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)

         write(t_,13040)
13040    format("nmodsa is max number of wave modes or harmonics for")
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,t_)

         write(t_,13041)
13041    format("a single mode. CHECK code,  for values .ne. 3.")
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,t_)

         write(t_,13042) nmodsa
13042    format("====> nmodsa = ",i5)
         RILIN=RILIN+1.
         CALL PGMTXT('T',-RILIN,0.,0.,t_)

!         write(t_,13043)
13043    format("nharma is max harmonic for cyclotron interactions")
!         RILIN=RILIN+1.
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)

!         write(t_,13044) nharma
13044    format("====> nharma = ",i5)
!         RILIN=RILIN+1.
!         CALL PGMTXT('T',-RILIN,0.,0.,t_)

      endif

      return
      end
end module ainpltpa_mod
