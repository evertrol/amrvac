!> Module with Roe-type Riemann solver for hydrodynamics
module mod_hd_roe
  use mod_hd_phys
  use mod_physics_roe

  implicit none
  private

  integer :: soundRW_ = -1
  integer :: soundLW_ = -1
  integer :: entropW_ = -1
  integer :: shearW0_ = -1

  public :: hd_roe_init

contains

  subroutine hd_roe_init()
    use mod_global_parameters, only: entropycoef, nw

    integer :: iw

    if (hd_energy) then
       ! Characteristic waves
       soundRW_ = 1
       soundLW_ = 2
       entropW_ = 3
       shearW0_ = 3
       nworkroe = 3

       phys_average => hd_average
       phys_get_eigenjump => hd_get_eigenjump
       phys_rtimes => hd_rtimes
    else
       ! Characteristic waves
       soundRW_ = 1
       soundLW_ = 2
       shearW0_ = 2
       nworkroe = 1

       phys_average => hd_average_iso
       phys_get_eigenjump => hd_get_eigenjump_iso
       phys_rtimes => hd_rtimes_iso
    end if

    allocate(entropycoef(nw))

    do iw = 1, nw
       if (iw == soundRW_ .or. iw == soundLW_) then
          ! TODO: Jannis: what's this?
          entropycoef(iw) = 0.2d0
       else
          entropycoef(iw) = -1.0d0
       end if
    end do

  end subroutine hd_roe_init

  !> Calculate the Roe average of w, assignment of variables:
  !> rho -> rho, m -> v, e -> h
  subroutine hd_average(wL,wR,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,&
     wroe,workroe)
    use mod_global_parameters
    integer, intent(in)             :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3, idim
    double precision, intent(in)    :: wL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nw), wR(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3, nw)
    double precision, intent(inout) :: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nw)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nworkroe)
    double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, 1:3)
    integer                         :: idir

    ! call average2(wL,wR,x,ix^L,idim,wroe,workroe(ixG^T,1),workroe(ixG^T,2))
    workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        1) = sqrt(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
    workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        2) = sqrt(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))

    ! The averaged density is sqrt(rhoL*rhoR)
    wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       rho_)  = workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        1)*workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3, 2)

    ! Now the ratio sqrt(rhoL/rhoR) is put into workroe(ix^S, 1)
    workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        1) = workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        1)/workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3, 2)

    ! Roe-average velocities
    do idir = 1, ndir
       wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mom(idir)) = (wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mom(idir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_) * workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           1)+wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mom(idir))/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_))/(one+workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3, 1))
    end do

    ! Calculate enthalpyL, then enthalpyR, then Roe-average. Use tmp2 for pressure.
    call hd_get_pthermal(wL,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixmin1,&
       ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, 2))

    wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       e_)    = (workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        2)+wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,e_))/wL(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,rho_)

    call hd_get_pthermal(wR,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixmin1,&
       ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, 2))

    workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        2) = (workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        2)+wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,e_))/wR(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,rho_)
    wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       e_)    = (wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       e_)*workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        1) + workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        2))/(one+workroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3, 1))
  end subroutine hd_average

  subroutine average2(wL,wR,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,&
     wroe,tmp,tmp2)

    ! Calculate the Roe average of w, assignment of variables:
    ! rho -> rho, m -> v, e -> h

    use mod_global_parameters

    integer                                             :: ixmin1,ixmin2,&
       ixmin3,ixmax1,ixmax2,ixmax3,idim,idir
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nw)               :: wL,wR,wroe
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       ndim), intent(in) :: x
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)                  :: tmp,tmp2
    !-----------------------------------------------------------------------------


  end subroutine average2

  !> Calculate the il-th characteristic speed and the jump in the il-th 
  !> characteristic variable in the idim direction within ixL. 
  !> The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
  !> jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))
  subroutine hd_get_eigenjump(wL,wR,wroe,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
     ixmax3,il,idim,smalla,a,jump,workroe)
    use mod_global_parameters

    integer, intent(in) :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,il,idim
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nw):: wL,wR,wroe
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       ndim), intent(in) :: x
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)   :: smalla,a,jump
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nworkroe) :: workroe

    call geteigenjump2(wL,wR,wroe,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
       il,idim,smalla,a,jump, workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,2),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,3))

  end subroutine hd_get_eigenjump

  subroutine geteigenjump2(wL,wR,wroe,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
     ixmax3,il,idim,smalla,a,jump,csound,dpperc2,dvperc)

    ! Calculate the il-th characteristic speed and the jump in the il-th 
    ! characteristic variable in the idim direction within ixL. 
    ! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
    ! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

    use mod_global_parameters
    use mod_tvd

    integer                                             :: ixmin1,ixmin2,&
       ixmin3,ixmax1,ixmax2,ixmax3,il,idim,idir
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nw)               :: wL,wR,wroe
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       ndim), intent(in) :: x
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)                  :: smalla,a,jump,tmp,tmp2
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)                  :: csound,dpperc2,dvperc
    double precision                                    :: &
       kin_en(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)

    if(il==1)then
       !First calculate the square of the sound speed: c**2=(gamma-1)*(h-0.5*v**2)
       kin_en(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3) = 0.5d0 * sum(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3, mom(:))**2, dim=3+1)
       csound(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)=(hd_gamma-one)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,e_) - kin_en(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))
       ! Make sure that csound**2 is positive
       csound(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)=max(hd_gamma*smalldouble/wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_),csound(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))

       ! Calculate (pR-pL)/c**2
       ! To save memory we use tmp amnd tmp2 for pL and pR (hd_get_pthermal is OK)
       call hd_get_pthermal(wL,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp)
       call hd_get_pthermal(wR,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp2)
       dpperc2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(tmp2(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)-tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))/csound(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)

       !Now get the correct sound speed
       csound(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)=sqrt(csound(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))

       ! Calculate (vR_idim-vL_idim)/c
       dvperc(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,mom(idim))/wL(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_))/csound(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)

    endif

    if (il == soundRW_) then
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))+csound(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half*(dpperc2(&
          ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)+wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)*dvperc(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))
    else if (il == soundLW_) then
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))-csound(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half*(dpperc2(&
          ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)-wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)*dvperc(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))
    else if (il == entropW_) then
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-dpperc2(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)+wR(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_)
    else
       !Determine the direction of the shear wave
       idir=il-shearW0_; if(idir>=idim)idir=idir+1
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,mom(idir))/wR(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mom(idir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
    end if

    ! Calculate "smalla" or modify "a" based on the "typeentropy" switch
    ! Put left and right eigenvalues, if needed, into tmp and tmp2
    ! OK, since subroutines hd_get_pthermal and entropyfix do not use tmp and tmp2

    select case(typeentropy(il))
    case('yee')
       ! Based on Yee JCP 68,151 eq 3.23
       smalla(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=entropycoef(il)
    case('harten','powell')
       ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
       if (il == soundRW_) then
          call hd_get_pthermal(wL,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp)
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)+ sqrt(hd_gamma*tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
          call hd_get_pthermal(wR,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp2)
          tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)+ sqrt(hd_gamma*tmp2(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
       else if (il == soundLW_) then
          call hd_get_pthermal(wL,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp)
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)- sqrt(hd_gamma*tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
          call hd_get_pthermal(wR,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp2)
          tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)- sqrt(hd_gamma*tmp2(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
       else
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) =wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,rho_)
          tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,rho_)
       end if
    end select

    call entropyfix(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,il,tmp,tmp2,a,&
       smalla)

  end subroutine geteigenjump2

  !> Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe
  subroutine hd_rtimes(q,wroe,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,il,&
     idim,rq,workroe)
    use mod_global_parameters

    integer, intent(in)             :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3,iw,il,idim
    double precision, intent(in)    :: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw)
    double precision, intent(in)    :: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)
    double precision, intent(inout) :: rq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nworkroe)
    !-----------------------------------------------------------------------------
    call rtimes2(q,wroe,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,il,idim,&
       rq,workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1))
  end subroutine hd_rtimes

  subroutine rtimes2(q,wroe,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,il,&
     idim,rq,csound)

    ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

    use mod_global_parameters

    integer                            :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3,iw,il,idim,idir
    double precision                   :: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: q,rq,csound
    logical                            :: shearwave

    shearwave=il>shearW0_
    idir=idim
    if(shearwave)then
       ! Direction of shearwave increases with il plus idir==idim is jumped over
       idir=il-shearW0_; if(idir>=idim)idir=idir+1
    endif

    if (iw == rho_) then
       if(shearwave)then 
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
       else
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
       endif
    else if (iw == e_) then
       if (il == soundRW_) then
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,e_)+wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim))*csound(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
       else if (il == soundLW_) then
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,e_)-wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim))*csound(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
       else if (il == entropW_) then
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) * 0.5d0 * sum(wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3, mom(:))**2, dim=3+1)
       else
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idir))
       end if
    else
       if(iw==mom(idim))then
          if (il == soundRW_) then
             rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
                ixmin2:ixmax2,ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3,mom(idim))+csound(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3))
          else if (il == soundLW_) then
             rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
                ixmin2:ixmax2,ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3,mom(idim))-csound(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3))
          else if (il == entropW_) then
             rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
                ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3,mom(idim))
          else
             rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
          end if
       else
          if(shearwave)then
             if(iw==mom(idir))then
                rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
                   ixmin2:ixmax2,ixmin3:ixmax3)
             else
                rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
             endif
          else
             rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
                ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3,iw)
          endif
       endif
    end if

  end subroutine rtimes2

  subroutine hd_average_iso(wL,wR,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
     idim,wroe,workroe)

    ! Calculate the Roe average of w, assignment of variables:
    ! rho -> rho, m -> v
    use mod_global_parameters

    integer, intent(in)             :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3, idim
    double precision, intent(in)    :: wL(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nw), wR(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3, nw)
    double precision, intent(inout) :: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nw)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nworkroe)
    double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, 1:3)

    call average2_iso(wL,wR,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,&
       wroe,workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1))

  end subroutine hd_average_iso

  subroutine average2_iso(wL,wR,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
     idim,wroe,tmp)

    ! Calculate the Roe average of w, assignment of variables:
    ! rho -> rho, m -> v

    use mod_global_parameters

    integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,idir
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nw):: wL,wR,wroe
    double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:ndim)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3):: tmp
    !-----------------------------------------------------------------------------

    select case (typeaverage)
    case ('arithmetic')
       ! This is the simple arithmetic average
       wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
       do idir = 1, ndir
          wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
              mom(idir)) =  half * (wL(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,rho_) + wR(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idir))/wR(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,rho_))
       end do
    case ('roe','default')
       ! Calculate the Roe-average
       wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_)=sqrt(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_)*wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
       ! Roe-average velocities
       tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=sqrt(wL(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)/wR(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_))
       do idir=1,ndir
          wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idir))=(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)*tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)+wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idir))/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))/(one+tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
       end do
    end select

  end subroutine average2_iso

  subroutine hd_get_eigenjump_iso(wL,wR,wroe,x,ixmin1,ixmin2,ixmin3,ixmax1,&
     ixmax2,ixmax3,il,idim,smalla,a,jump,workroe)

    ! Calculate the il-th characteristic speed and the jump in the il-th 
    ! characteristic variable in the idim direction within ixL. 
    ! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
    ! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

    use mod_global_parameters

    integer, intent(in) :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,il,idim
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nw):: wL,wR,wroe
    double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:ndim)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)   :: smalla,a,jump
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nworkroe) :: workroe
    !-----------------------------------------------------------------------------
    call geteigenjump2_iso(wL,wR,wroe,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3,il,idim,smalla,a,jump,workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1))

  end subroutine hd_get_eigenjump_iso

  subroutine geteigenjump2_iso(wL,wR,wroe,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
     ixmax3,il,idim,smalla,a,jump,csound)

    ! Calculate the il-th characteristic speed and the jump in the il-th 
    ! characteristic variable in the idim direction within ixL. 
    ! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
    ! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

    use mod_global_parameters
    use mod_tvd

    integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,il,idim,idir
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nw):: wL,wR,wroe
    double precision, intent(in)    :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:ndim)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)   :: smalla,a,jump,tmp,tmp2
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)   :: csound
    DOUBLE PRECISION,PARAMETER:: qsmall=1.D-6
    !-----------------------------------------------------------------------------

    select case (typeaverage)
    case ('arithmetic')
       call hd_get_pthermal(wroe,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,csound)
       csound(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3) = sqrt(hd_gamma*csound(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)/wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
       ! This is the original simple Roe-solver
       if (il == soundRW_) then
          a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))+csound(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
          jump(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=half*((one-wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idim))/csound(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3))*(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))+(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim))-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim)))/csound(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
       else if (il == soundLW_) then
          a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))-csound(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
          jump(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=half*((one+wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idim))/csound(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3))*(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))-(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim))-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim)))/csound(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
       else
          ! Determine direction of shear wave
          idir=il-shearW0_; if(idir>=idim)idir=idir+1
          a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))
          jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idir))*(wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,rho_))+(wR(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idir))-wL(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idir)))
       end if
    case ('roe','default')
       call hd_get_pthermal(wroe,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,csound)
       call hd_get_pthermal(wL,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp)
       call hd_get_pthermal(wR,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp2)
       where(abs(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_)-wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_))<=qsmall*(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)))
          csound(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3) = sqrt(hd_gamma*csound(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)/wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
       elsewhere
          csound(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3) =  sqrt(hd_gamma*(tmp2(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)-tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3))/(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)))
       end where
       ! This is the Roe solver by Glaister
       ! based on P. Glaister JCP 93, 477-480 (1991)
       if (il == soundRW_) then
          a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))+csound(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
          jump(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=half*((wR(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))+wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)/csound(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim))/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim))/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)))
       else if (il == soundLW_) then
          a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))-csound(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
          jump(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=half*((wR(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))-wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)/csound(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim))/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idim))/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)))
       else
          ! Determine direction of shear wave
          idir=il-shearW0_; if(idir>=idim)idir=idir+1
          a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))
          jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,rho_)*(wR(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idir))/wR(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mom(idir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
       end if
    end select

    ! Calculate "smalla" or modify "a" based on the "typeentropy" switch
    ! Use tmp and tmp2 for the left and right eigenvalues if needed
    select case(typeentropy(il))
    case('yee')
       ! Based on Yee JCP 68,151 eq 3.23
       smalla(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=entropycoef(il)
    case('harten','powell')
       ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
       if (il == soundRW_) then
          call hd_get_pthermal(wL,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp)
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) =wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)+ sqrt(hd_gamma*tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
          call hd_get_pthermal(wR,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp2)
          tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)+ sqrt(hd_gamma*tmp2(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
       else if (il == soundLW_) then
          call hd_get_pthermal(wL,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp)
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) =wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)- sqrt(hd_gamma*tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
          call hd_get_pthermal(wR,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,tmp2)
          tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_)- sqrt(hd_gamma*tmp2(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
       else
          tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) =wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,rho_)
          tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,rho_)
       end if
    end select

    call entropyfix(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,il,tmp,tmp2,a,&
       smalla)

  end subroutine geteigenjump2_iso

  subroutine hd_rtimes_iso(q,wroe,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,&
     il,idim,rq,workroe)

    ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

    use mod_global_parameters

    integer, intent(in)             :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3,iw,il,idim
    double precision, intent(in)    :: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw)
    double precision, intent(in)    :: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)
    double precision, intent(inout) :: rq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nworkroe)
    !-----------------------------------------------------------------------------

    call rtimes2_iso(q,wroe,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,il,&
       idim,rq,workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1))

  end subroutine hd_rtimes_iso

  subroutine rtimes2_iso(q,wroe,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,&
     il,idim,rq,csound)

    ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

    use mod_global_parameters

    integer::          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,il,idim,&
       idir
    double precision:: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,nw)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3):: q,&
       rq,csound

    if(iw==rho_)then
       if (il == soundRW_ .or. il == soundLW_) then
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
       else
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
       end if
    else if(iw==mom(idim))then
       if (il == soundRW_) then
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idim))+csound(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3))
       else if (il == soundLW_) then
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mom(idim))-csound(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3))
       else
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
       end if
    else
       if (il == soundRW_ .or. il == soundLW_) then
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,iw)
       else
          !Determine direction of shear wave
          idir=il-shearW0_; if(idir>=idim)idir=idir+1
          if(iw==mom(idir)) then
             rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
                ixmin2:ixmax2,ixmin3:ixmax3)
          else
             rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
          endif
       end if
    endif

  end subroutine rtimes2_iso

end module mod_hd_roe
