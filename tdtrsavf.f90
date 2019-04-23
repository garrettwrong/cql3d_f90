module tdtrsavf_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine tdtrsavf
      use param_mod
      use comm_mod
      implicit integer (i-n), real*8 (a-h,o-z)

!..............................................................
!     This routine saves various versions of the distribution
!     function.
!     After transport the distribution is stored in frn.
!     After velocity advancement the distribution is stored
!     in fvn.
!     One time step back we have frn_1 and fvn_1 respectively.
!     We store in f either (frn_1+fvn_1+fvn*frn)*.25 or
!     (fvn+frn)*.5
!     f is used to compute collisional coefficients and to update
!     macroscopic quantities, such as temperature, current.
!     For purposes of the next velocity step we need an old version
!     of f, frn_2. We use either (frn+frn_1)*.5====>frn_2 or we use
!     frn=====>frn_2.
!     Next we update distributions: fvn===>fvn_1 and
!     frn====>frn_1
!     Then we copy fvn====>fvn_1
!..............................................................



      do 10 k=1,ngen
        do 20 j=1,jx
          do 30 l=1,lrors
            if (relaxtsp.eq."enabled" .and. nonadi.eq.2 .and. n.ge.2) &
              then
              do 41 i=1,iy_(l)
                f(i,j,k,l)=(frn(i,j,k,l)+fvn(i,j,k,l)+frn_1(i,j,k,l) &
                  +fvn_1(i,j,k,l))*.25
                frn_2(i,j,k,l)=(frn(i,j,k,l)+frn_1(i,j,k,l))*.5
!%OS  f(i,j,k,l)=frn_2(i,j,k,l)
 41           continue
              go to 42
            endif
!%OS
            if (relaxtsp.eq."enabled" .and. n.ge.2) then
              do 40 i=1,iy_(l)
!%OS  f(i,j,k,l)=(frn(i,j,k,l)+fvn(i,j,k,l)+frn_1(i,j,k,l)
!%OS  1            +fvn_1(i,j,k,l))*.25
                frn_2(i,j,k,l)=(frn(i,j,k,l)+frn_1(i,j,k,l))*.5
                f(i,j,k,l)=frn_2(i,j,k,l)
 40           continue
            else
              do 50 i=1,iy_(l)
!-out for now...f(i,j,k,l)=(frn(i,j,k,l)+fvn(i,j,k,l))*.5
                f(i,j,k,l)=frn(i,j,k,l)
                frn_2(i,j,k,l)=frn(i,j,k,l)
 50           continue
            endif
!%OS
 42         continue
!%OS
            do 60 i=1,iy_(l)
              frn_1(i,j,k,l)=frn(i,j,k,l)
              fvn_1(i,j,k,l)=fvn(i,j,k,l)
 60         continue
 30       continue
 20     continue
 10   continue
      return
      end
end module tdtrsavf_mod
