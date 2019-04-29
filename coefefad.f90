module coefefad_mod

  !---BEGIN USE

  !---END USE

!
!

contains

      subroutine coefefad(k)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)

!..................................................................
!     This routine adds in the contribution of the D.C. electric
!     field to the coefficients used for time advancement.
!..................................................................


!.......................................................................
!l    Include electrostatic field
!.......................................................................

      do 20 i=1,iy
        ztrda=-elparnw(l_)*charge/vnorm*coss(i,l_)
        ztrdd=elparnw(l_)*charge/vnorm*sinn(i,l_)**2
        do 21 j=1,jx
          da(i,j)=da(i,j)+bnumb(k)/fmass(k)*(elecfld(lr_)*cex(i,j,l_)+ &
            ztrda*xsq(j))
          dd(i,j)=dd(i,j)+bnumb(k)/fmass(k)*(elecfld(lr_)*cet(i,j,l_)+ &
            ztrdd*x(j))
 21     continue
 20   continue

      return
      end
end module coefefad_mod
