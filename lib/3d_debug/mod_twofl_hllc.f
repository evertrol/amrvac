module mod_twofl_hllc
#include "amrvac.h"
  use mod_twofl_phys
  use mod_physics_hllc
  use mod_functions_bfield, only: mag

  implicit none
  private


  public :: twofl_hllc_init

  integer :: i_rho, i_e
  integer, allocatable :: i_mom(:)


contains

  subroutine twofl_hllc_init()
    use mod_physics_hllc
    use mod_global_parameters

    allocate(i_mom(1:ndir))

    phys_hllc_init_species => twofl_hllc_init_species

  end subroutine twofl_hllc_init

  subroutine twofl_hllc_init_species(ii, rho_, mom, e_)
    use mod_global_parameters
    integer, intent(in) :: ii
    integer, intent(out) :: rho_, mom(1:ndir),e_

    if(ii==1) then
      phys_diffuse_hllcd => twofl_diffuse_hllcd_c
      phys_get_lCD => twofl_get_lCD_c
      phys_get_wCD => twofl_get_wCD_c
      i_rho = rho_c_
      i_mom(1:ndir) = mom_c(1:ndir)
      i_e = e_c_
    else
      phys_diffuse_hllcd => twofl_diffuse_hllcd_n
      phys_get_lCD => twofl_get_lCD_n
      phys_get_wCD => twofl_get_wCD_n
      i_rho = rho_n_
      i_mom(1:ndir) = mom_n(1:ndir)
      i_e = e_n_
    end if

    rho_ = i_rho
    mom(:) = i_mom(:)
    e_ = i_e

  end subroutine twofl_hllc_init_species

  subroutine twofl_diffuse_hllcd_n(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,wLC,wRC,fLC,&
     fRC,patchf)

    ! when method is hllcd or hllcd1 then: 

    ! this subroutine is to enforce regions where we AVOID HLLC
    ! and use TVDLF instead: this is achieved by setting patchf to 4 in
    ! certain regions. An additional input parameter is nxdiffusehllc
    ! which sets the size of the fallback region.

    use mod_global_parameters

    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), intent(in)      :: wRC,wLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux),intent(in) :: fLC, fRC
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(inout)        :: patchf

    integer                                           :: ixOO1,ixOO2,ixOO3,&
       TxOOmin1,TxOOmin2,TxOOmin3,TxOOmax1,TxOOmax2,TxOOmax3

    ! In a user-controlled region around any point with flux sign change between
    ! left and right, ensure fallback to TVDLF
    do ixOO1= ixOmin1,ixOmax1
    do ixOO2= ixOmin2,ixOmax2
    do ixOO3= ixOmin3,ixOmax3
    
    TxOOmin1= max(ixOO1 - nxdiffusehllc*kr(idim,1), ixOmin1);
    TxOOmax1= min(ixOO1 + nxdiffusehllc*kr(idim,1), ixOmax1);
    
    
    TxOOmin2= max(ixOO2 - nxdiffusehllc*kr(idim,2), ixOmin2);
    TxOOmax2= min(ixOO2 + nxdiffusehllc*kr(idim,2), ixOmax2);
    
    
    TxOOmin3= max(ixOO3 - nxdiffusehllc*kr(idim,3), ixOmin3);
    TxOOmax3= min(ixOO3 + nxdiffusehllc*kr(idim,3), ixOmax3);
    
    if(abs(patchf(ixOO1,ixOO2,ixOO3)) == 1 .or. abs(patchf(ixOO1,ixOO2,&
       ixOO3)) == 4)Then
       if(any(fRC(ixOO1,ixOO2,ixOO3,1:nwflux)*fLC(ixOO1,ixOO2,ixOO3,&
          1:nwflux)<-smalldouble))Then
          where(Abs(patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2,&
             TxOOmin3:TxOOmax3))==1)
             patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2,TxOOmin3:TxOOmax3) = 4
          endwhere
       endif
    endif
    enddo
    enddo
    enddo

  end subroutine twofl_diffuse_hllcd_n

  subroutine twofl_get_lCD_n(wLC,wRC,fLC,fRC,cmin,cmax,idim,ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, whll,Fhll,lambdaCD,patchf)
    ! Calculate lambda at CD and set the patchf to know the orientation
    ! of the riemann fan and decide on the flux choice
    ! We also compute here the HLL flux and w value, for fallback strategy

    ! In this nul version, we simply compute nothing and ensure TVDLF fallback
    use mod_global_parameters

    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), intent(in)      :: wLC,wRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux), intent(in)  :: fLC,fRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(in)           :: cmax,cmin
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(inout)        :: patchf
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux), intent(out) :: Fhll,whll
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(out)          :: lambdaCD

    logical         , dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)     :: Cond_patchf
    double precision                       :: Epsilon
    integer                                :: iw,ix1,ix2,ix3

    !--------------------------------------------
    ! on entry, patch is preset to contain values from -2,1,2,4
    !      -2: take left flux, no computation here
    !      +2: take right flux, no computation here
    !      +4: take TVDLF flux, no computation here
    !       1: compute the characteristic speed for the CD

    Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))==1)
    lambdaCD=0.d0

    do iw=1,nwflux
      where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
        !============= compute HLL flux ==============!
        Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw)= (cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*fRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw) + cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw)-wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw)))/(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
        !======== compute intermediate HLL state =======!
        whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw) = (cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw)+fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw)-fRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw))/(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
      end where
    end do

    ! deduce the characteristic speed at the CD
    where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
       lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,i_mom(idim))/whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,i_rho)
    end where

    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
       if(Cond_patchf(ix1,ix2,ix3)) then
         ! double check whether obtained speed is in between min and max speeds given
         ! and identify in which part of the Riemann fan the time-axis is
         if(cmin(ix1,ix2,ix3) < zero .and. lambdaCD(ix1,ix2,&
            ix3)>zero.and.lambdaCD(ix1,ix2,ix3)<cmax(ix1,ix2,ix3)) then
            patchf(ix1,ix2,ix3) = -1
         else if(cmax(ix1,ix2,ix3) > zero .and. lambdaCD(ix1,ix2,&
            ix3) < zero.and.lambdaCD(ix1,ix2,ix3)>cmin(ix1,ix2,ix3)) then
            patchf(ix1,ix2,ix3) =  1
         else if(lambdaCD(ix1,ix2,ix3) >= cmax(ix1,ix2,ix3) .or. lambdaCD(ix1,&
            ix2,ix3) <= cmin(ix1,ix2,ix3)) then
            lambdaCD(ix1,ix2,ix3) = zero
            ! we will fall back to HLL flux case in this degeneracy
            patchf(ix1,ix2,ix3) =  3
         end if
       end if
    end do
    end do
    end do

    where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)== 3)
      Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=.false.
    end where

    ! handle the specific case where the time axis is exactly on the CD 
    if(any(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3).and.lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)==zero))then
      ! determine which sector (forward or backward) of the Riemann fan is smallest
      ! and select left or right flux accordingly
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        if(lambdaCD(ix1,ix2,ix3)==zero.and.Cond_patchf(ix1,ix2,ix3)) then
          if(-cmin(ix1,ix2,ix3)>=cmax(ix1,ix2,ix3)) then
            patchf(ix1,ix2,ix3) =  1
          else
            patchf(ix1,ix2,ix3) = -1
          end if
        end if
      end do
      end do
      end do
    end if

    ! eigenvalue lambda for contact is near zero: decrease noise by this trick
    if(flathllc)then
      Epsilon=1.0d-6
      where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3).and. dabs(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/max(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),Epsilon)< Epsilon  .and. &
         dabs(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/max(dabs(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)),Epsilon)< Epsilon)
        lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  zero
      end where
    end if

  end subroutine twofl_get_lCD_n

  subroutine twofl_get_wCD_n(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,&
     cmax,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
     ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,f)
    ! compute the intermediate state U*
    ! only needed where patchf=-1/1

    ! reference Li S., JCP, 203, 2005, 344-357
    ! reference T. Miyoski, Kusano JCP, 2008, 2005.

    use mod_global_parameters
    use mod_physics

    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), intent(in)      :: wRC,wLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux), intent(in):: whll, Fhll
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(in)           :: lambdaCD
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(in)           :: cmax,cmin
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux), intent(in):: fRC,fLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux),intent(out):: f
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)      :: wCD,wSub
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux)  :: fSub
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)           :: vSub,cspeed,pCD
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(in)           :: patchf

    double precision :: csmls
    integer                                      :: n, iw, ix1,ix2,ix3

    !-------------- auxiliary Speed and array-------------!
    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
       if(patchf(ix1,ix2,ix3)==1) then
         cspeed(ix1,ix2,ix3)=cmax(ix1,ix2,ix3)
         vSub(ix1,ix2,ix3)=wRC(ix1,ix2,ix3,i_mom(idim))/wRC(ix1,ix2,ix3,i_rho)
         wSub(ix1,ix2,ix3,:)=wRC(ix1,ix2,ix3,:)
         fSub(ix1,ix2,ix3,:)=fRC(ix1,ix2,ix3,:)
       else if(patchf(ix1,ix2,ix3)==-1) then
         cspeed(ix1,ix2,ix3)=cmin(ix1,ix2,ix3)
         vSub(ix1,ix2,ix3)=wLC(ix1,ix2,ix3,i_mom(idim))/wLC(ix1,ix2,ix3,i_rho)
         wSub(ix1,ix2,ix3,:)=wLC(ix1,ix2,ix3,:)
         fSub(ix1,ix2,ix3,:)=fLC(ix1,ix2,ix3,:)
       end if
    end do
    end do
    end do

    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
      if(abs(patchf(ix1,ix2,ix3))==1) then
        csmls=one/(cspeed(ix1,ix2,ix3)-lambdaCD(ix1,ix2,ix3))
        wCD(ix1,ix2,ix3,i_rho) = wSub(ix1,ix2,ix3,i_rho)*(cspeed(ix1,ix2,&
           ix3)-vSub(ix1,ix2,ix3))*csmls

        !------- Momentum ------!
        do iw=1, ndir
          if(iw /= idim)then
            ! eq. 21 22
            wCD(ix1,ix2,ix3,i_mom(iw))=(cspeed(ix1,ix2,ix3)*wSub(ix1,ix2,ix3,&
               i_mom(iw))-fSub(ix1,ix2,ix3,i_mom(iw)))*csmls
          else
            ! eq. 20
            wCD(ix1,ix2,ix3,i_mom(iw)) =  wCD(ix1,ix2,ix3,&
               i_rho) * lambdaCD(ix1,ix2,ix3)
          endif
        enddo
        if(phys_energy) then
          pCD(ix1,ix2,ix3) = wSub(ix1,ix2,ix3,i_rho)*(cspeed(ix1,ix2,&
             ix3)-vSub(ix1,ix2,ix3))*(lambdaCD(ix1,ix2,ix3)-vSub(ix1,ix2,&
             ix3))+fSub(ix1,ix2,ix3,i_mom(idim))-wSub(ix1,ix2,ix3,&
             i_mom(idim))*vSub(ix1,ix2,ix3)
          ! Eq 31
          wCD(ix1,ix2,ix3,i_e) = (cspeed(ix1,ix2,ix3)*wSub(ix1,ix2,ix3,&
             i_e) -fSub(ix1,ix2,ix3,i_e)+lambdaCD(ix1,ix2,ix3)*pCD(ix1,ix2,&
             ix3))*csmls
        end if
        do iw=1,nwflux
          ! f_i=fsub+lambda (wCD-wSub)
          f(ix1,ix2,ix3,iw)=fsub(ix1,ix2,ix3,iw)+cspeed(ix1,ix2,ix3)*(wCD(ix1,&
             ix2,ix3,iw)-wSub(ix1,ix2,ix3,iw))
        end do
      end if
    end do
    end do
    end do

  end subroutine twofl_get_wCD_n

  subroutine twofl_diffuse_hllcd_c(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,wLC,wRC,fLC,&
     fRC,patchf)
  ! when method is hllcd or hllcd1 then: 
  ! this subroutine is to enforce regions where we AVOID HLLC
  ! and use TVDLF instead: this is achieved by setting patchf to 4 in
  ! certain regions. An additional input parameter is nxdiffusehllc
  ! which sets the size of the fallback region.
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), intent(in)      :: wRC,wLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux),intent(in) :: fLC, fRC
    integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        intent(inout) :: patchf
    
    integer :: ixOO1,ixOO2,ixOO3,TxOOmin1,TxOOmin2,TxOOmin3,TxOOmax1,TxOOmax2,&
       TxOOmax3

    
    ! In a user-controlled region around any point with flux sign change between
    ! left and right, ensure fallback to TVDLF
    do ixOO1= ixOmin1,ixOmax1
    do ixOO2= ixOmin2,ixOmax2
    do ixOO3= ixOmin3,ixOmax3
      
      TxOOmin1= max(ixOO1 - nxdiffusehllc*kr(idim,1), ixOmin1);
      TxOOmax1= min(ixOO1 + nxdiffusehllc*kr(idim,1), ixOmax1);
      
      
      TxOOmin2= max(ixOO2 - nxdiffusehllc*kr(idim,2), ixOmin2);
      TxOOmax2= min(ixOO2 + nxdiffusehllc*kr(idim,2), ixOmax2);
      
      
      TxOOmin3= max(ixOO3 - nxdiffusehllc*kr(idim,3), ixOmin3);
      TxOOmax3= min(ixOO3 + nxdiffusehllc*kr(idim,3), ixOmax3);
      
      if(abs(patchf(ixOO1,ixOO2,ixOO3)) == 1 .or. abs(patchf(ixOO1,ixOO2,&
         ixOO3)) == 4)Then
         if(any(fRC(ixOO1,ixOO2,ixOO3,1:nwflux)*fLC(ixOO1,ixOO2,ixOO3,&
            1:nwflux)<-smalldouble))Then
           where(Abs(patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2,&
              TxOOmin3:TxOOmax3))==1)
             patchf(TxOOmin1:TxOOmax1,TxOOmin2:TxOOmax2,TxOOmin3:TxOOmax3) = 4
           endwhere
         endif
      endif
    enddo
    enddo
    enddo
  
  end subroutine twofl_diffuse_hllcd_c

  subroutine twofl_get_lCD_c(wLC,wRC,fLC,fRC,cmin,cmax,idim,ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, whll,Fhll,lambdaCD,patchf)
  
  ! Calculate lambda at CD and set the patchf to know the orientation
  ! of the riemann fan and decide on the flux choice
  ! We also compute here the HLL flux and w value, for fallback strategy
  
    use mod_global_parameters
    
    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), intent(in)      :: wLC,wRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux), intent(in):: fLC,fRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(in)           :: cmax,cmin
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(inout)        :: patchf
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux), intent(out) :: Fhll,whll
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(out)            :: lambdaCD
    
    logical         , dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)     :: Cond_patchf
    double precision                       :: Epsilon
    integer                                :: iw

    ! on entry, patch is preset to contain values from -2,1,2,4
    !      -2: take left flux, no computation here
    !      +2: take right flux, no computation here
    !      +4: take TVDLF flux, no computation here
    !       1: compute the characteristic speed for the CD
    
    Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))==1)
    
    do iw=1,nwflux
      where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      !============= compute HLL flux ==============!
      Fhll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw)= (cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*fRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw) + cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*(wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw)-wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw)))/(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      !======== compute intermediate HLL state =======!
      whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw) = (cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw)+fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw)-fRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw))/(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      endwhere
    enddo
    
    ! deduce the characteristic speed at the CD
    where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         i_mom(idim))/whll(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         i_rho)
    end where
    
    
    where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      ! double check whether obtained speed is in between min and max speeds given
      ! and identify in which part of the Riemann fan the time-axis is
      where(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)<zero.and.lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)>zero.and.lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)<cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
        patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = -1
      elsewhere(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)>zero.and.lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)<zero.and.lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)>cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
        patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  1
      elsewhere(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)>=cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3).or.lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) <= cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
        lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero
        ! we will fall back to HLL flux case in this degeneracy
        patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  3
      endwhere
    endwhere
    
    
    where(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)== 3)
      Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=.false.
    end where
    
    
    ! handle the specific case where the time axis is exactly on the CD 
    if(any(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)==zero.and.Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))then
      ! determine which sector (forward or backward) of the Riemann fan is smallest
      ! and select left or right flux accordingly
      where(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)==zero.and.Cond_patchf(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3))
        where(-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)>=cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
          patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  1
        elsewhere
          patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = -1
        endwhere
      endwhere
    endif
    
    ! eigenvalue lambda for contact is near zero: decrease noise by this trick
    if(flathllc)then
      Epsilon=1.0d-6
      where(Cond_patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3).and. dabs(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/max(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),Epsilon)< Epsilon  .and. &
         dabs(lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/max(dabs(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)),Epsilon)< Epsilon)
        lambdaCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) =  zero
      end where
    end if 

  end subroutine twofl_get_lCD_c

  subroutine twofl_get_wCD_c(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cmin,&
     cmax,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
     ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,f)
  ! compute the intermediate state U*
  ! only needed where patchf=-1/1
  
  ! reference Li S., JCP, 203, 2005, 344-357
  ! reference T. Miyoski, Kusano JCP, 2008, 2005.
    use mod_global_parameters
    use mod_physics
    integer, intent(in)                                      :: ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), intent(in)      :: wRC,wLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux), intent(in):: whll, Fhll
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(in)           :: lambdaCD
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(in)           :: cmax,cmin
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux), intent(in):: fRC,fLC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux),intent(out):: f
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)        :: wCD,wSub
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux)    :: fSub
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)             :: vSub,cspeed,pCD,VdotBCD
    integer         , dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), intent(in)           :: patchf

    integer                                        :: n, iw, idir,ix1,ix2,ix3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)           :: rho
    

    !-------------- auxiliary Speed and array-------------!
    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
       if(patchf(ix1,ix2,ix3)==1) then
         cspeed(ix1,ix2,ix3)=cmax(ix1,ix2,ix3)
         vSub(ix1,ix2,ix3)=wRC(ix1,ix2,ix3,i_mom(idim))/wRC(ix1,ix2,ix3,i_rho)
         wSub(ix1,ix2,ix3,:)=wRC(ix1,ix2,ix3,:)
         fSub(ix1,ix2,ix3,:)=fRC(ix1,ix2,ix3,:)
       else if(patchf(ix1,ix2,ix3)==-1) then
         cspeed(ix1,ix2,ix3)=cmin(ix1,ix2,ix3)
         vSub(ix1,ix2,ix3)=wLC(ix1,ix2,ix3,i_mom(idim))/wLC(ix1,ix2,ix3,i_rho)
         wSub(ix1,ix2,ix3,:)=wLC(ix1,ix2,ix3,:)
         fSub(ix1,ix2,ix3,:)=fLC(ix1,ix2,ix3,:)
       end if
    end do
    end do
    end do
    
    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
      if(abs(patchf(ix1,ix2,ix3))==1) then
        wCD(ix1,ix2,ix3,i_rho) = wSub(ix1,ix2,ix3,i_rho)*(cspeed(ix1,ix2,&
           ix3)-vSub(ix1,ix2,ix3))/(cspeed(ix1,ix2,ix3)-lambdaCD(ix1,ix2,ix3))
         !TODO tracer 
!        do n=1,twofl_n_tracer
!          iw = tracer(n)
!          wCD(ix^D,iw) = wSub(ix^D,iw)*(cspeed(ix^D)-vSub(ix^D))&
!             /(cspeed(ix^D)-lambdaCD(ix^D))
!        end do
        !==== Magnetic field ====!
        do idir=1,ndir
          ! case from eq 31
          wCD(ix1,ix2,ix3,mag(idir)) = whll(ix1,ix2,ix3,mag(idir))
        end do
        !------- Momentum ------!
        do iw=1, ndir
          if(iw /= idim)then
            ! eq. 21 22
            wCD(ix1,ix2,ix3,i_mom(iw))=(cspeed(ix1,ix2,ix3)*wSub(ix1,ix2,ix3,&
               i_mom(iw))-fSub(ix1,ix2,ix3,i_mom(iw))-wCD(ix1,ix2,ix3,&
               mag(idim))*wCD(ix1,ix2,ix3,mag(iw)))/(cspeed(ix1,ix2,&
               ix3)-lambdaCD(ix1,ix2,ix3))
          else
            ! eq. 20
            wCD(ix1,ix2,ix3,i_mom(iw)) =  wCD(ix1,ix2,ix3,&
               i_rho) * lambdaCD(ix1,ix2,ix3)
          endif
        enddo
        if(phys_energy) then
          VdotBCD(ix1,ix2,ix3) = sum(whll(ix1,ix2,ix3,i_mom(:))*whll(ix1,ix2,&
             ix3,mag(:)))/whll(ix1,ix2,ix3,i_rho)
          ! Eq 17
          pCD(ix1,ix2,ix3)  = wsub(ix1,ix2,ix3,i_rho)*(cspeed(ix1,ix2,&
             ix3)-vSub(ix1,ix2,ix3))*(lambdaCD(ix1,ix2,ix3)-vSub(ix1,ix2,&
             ix3))+fSub(ix1,ix2,ix3,i_mom(idim))-wsub(ix1,ix2,ix3,&
             i_mom(idim))*vSub(ix1,ix2,ix3)+ wCD(ix1,ix2,ix3,mag(idim))**2
          ! Eq 31
          wCD(ix1,ix2,ix3,i_e) = (cspeed(ix1,ix2,ix3)*wSub(ix1,ix2,ix3,&
             i_e) -fSub(ix1,ix2,ix3,i_e)+lambdaCD(ix1,ix2,ix3)*pCD(ix1,ix2,&
             ix3)-VdotBCD(ix1,ix2,ix3)*wCD(ix1,ix2,ix3,mag(idim)))/(cspeed(ix1,&
             ix2,ix3)-lambdaCD(ix1,ix2,ix3))
        end if
      end if
    end do
    end do
    end do

    do iw=1,nwflux
     if(iw == mag(idim)) then
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)=zero
     else if(twofl_glm .and. iw == psi_) then
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw)=zero
     else
       where(abs(patchf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))==1)
         ! f_i=fsub+lambda (wCD-wSub)
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            iw)=fsub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            iw)+cspeed(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*(wCD(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,iw)-wsub(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,iw))
       endwhere
     end if
    end do

  end subroutine twofl_get_wCD_c

end module mod_twofl_hllc
