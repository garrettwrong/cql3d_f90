module sigmaxwl_mod

  !---BEGIN USE

  use bcast_mod, only : bcast

  !---END USE

!
!

contains

      subroutine sigmaxwl(k,kk)
      use param_mod
      use comm_mod
      implicit integer (i-n), real(c_double) (a-h,o-z)
!
!mnt  generate a maxwellian species "k" normalized to the general
!mnt  species "kk".


      fxllm2=2./3.*energy(k,lr_)/fions(k)
      fxppm2=fxllm2
      do 10 i=1,iy
      do 10 j=1,jx
      xll=x(j)*coss(i,l_)
      xpp=x(j)*sinn(i,l_)
      facx=-xll**2/fxllm2
      facy=-xpp**2/fxppm2
      temp3(i,j)=exp(facx+facy)
10    continue
      call bcast(tam1,zero,jx)
      do 20 i=1,iy
      do 21 j=1,jx
      tam1(j)=tam1(j)+temp3(i,j)*cynt2(i,l_)
21    continue
20    continue
      gn=0.
      do 22 j=1,jx
      gn=gn+tam1(j)*cint2(j)
22    continue
      sn=reden(k,lr_)/gn
      if (kk .eq. 0) go to 30
      sn=reden(k,lr_)/(reden(kk,lr_)*gn)
30    continue
      do 23 i=0,iyjx2-1
      temp3(i,0)=temp3(i,0)*sn
23    continue
      return
      end
end module sigmaxwl_mod
