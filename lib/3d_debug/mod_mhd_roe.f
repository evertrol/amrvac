!> Subroutines for Roe-type Riemann solver for MHD
module mod_mhd_roe
  use mod_mhd_phys
  use mod_functions_bfield, only: mag
  use mod_physics_roe

  implicit none
  private

  integer, parameter :: fastRW_ = 3,fastLW_=4,slowRW_=5,slowLW_=6 !Characteristic
  integer, parameter :: entroW_ = 8,diverW_=7,alfvRW_=1,alfvLW_=2 ! waves

  public :: mhd_roe_init

contains

  subroutine mhd_roe_init()
    use mod_global_parameters, only: entropycoef, nw
    integer :: il

    phys_average => mhd_average
    phys_get_eigenjump => mhd_get_eigenjump
    phys_rtimes => mhd_rtimes

    nworkroe=15
    allocate(entropycoef(nw))

    do il = 1, nw
       select case(il)
       case(fastRW_,fastLW_,slowRW_,slowLW_)
          entropycoef(il) = 0.2d0
       case(alfvRW_,alfvLW_)
          entropycoef(il) = 0.4d0
       case default
          entropycoef(il) = -1.0d0
       end select
    end do

  end subroutine mhd_roe_init

  ! Eight-wave MHD Riemann solver. See Powell, Notes on the eigensystem, Gombosi
  ! Calculate the wroe average of primitive variables in wL and wR, assignment:
  ! rho -> sqrho, m -> v, e -> p, B_idim -> B_idim, B_idir -> beta_idir
  ! Calculate also alpha_f,alpha_s,c_f,c_s,csound2,dp,rhodv
  !
  ! wL,wR,wroe are all interface centered quantities
  subroutine mhd_average(wL,wR,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
     idim,wroe,workroe)
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

    call average2(wL,wR,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,wroe,&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,2),&
        workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,3),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,4),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,5),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,6),&
        workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,7),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,8))

  end subroutine mhd_average

  ! Eight-wave MHD Riemann solver. See Powell, Notes on the eigensystem, Gombosi
  ! Calculate the wroe average of primitive variables in wL and wR, assignment:
  ! rho -> sqrho, m -> v, e -> p, B_idim -> B_idim, B_idir -> beta_idir
  ! Calculate also alpha_f,alpha_s,c_f,c_s,csound2,dp,rhodv
  !
  ! wL,wR,wroe are all interface centered quantities
  subroutine average2(wL,wR,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,&
     wroe,cfast,cslow,afast,aslow,csound2,dp, rhodv,tmp)
    use mod_global_parameters

    integer                               :: ixmin1,ixmin2,ixmin3,ixmax1,&
       ixmax2,ixmax3,idim,idir,jdir,iw
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nw) :: wL,wR,wroe
    double precision, intent(in)          :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:ndim)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)    :: cfast,cslow,afast,aslow,csound2,dp, rhodv,tmp

    if (ndir==1) call mpistop("MHD with d=11 is the same as HD")

    !Averaging primitive variables
    wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       rho_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
    do idir=1,ndir
      wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         mom(idir))=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         mom(idir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         mom(idir))/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
      wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         mag(idir))=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         mag(idir))+wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(idir)))
    end do
    ! Use afast and aslow for pressures pL and pR
    call mhd_get_pthermal(wL,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,afast)
    call mhd_get_pthermal(wR,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
       ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,aslow)

    if(mhd_energy) then
      wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         e_)=half*(afast(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)+aslow(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
      ! dp=pR-pL
      dp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aslow(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)-afast(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)
    else
      dp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aslow(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)-afast(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)
    end if

    !CONSERVATIVE rho*dv_idim=dm_idim-v_idim*drho
    rhodv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))-wL(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3,mom(idim))-wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       mom(idim))*(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))

    !Calculate csound2,cfast,cslow,alphafast and alphaslow

    ! get csound**2
    if(mhd_energy) then
      csound2(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)=mhd_gamma*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3,p_)/wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         rho_)
    else
      csound2(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3)=mhd_gamma*mhd_adiab*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
         ixmin3:ixmax3,rho_)**(mhd_gamma-one)
    end if

    ! aa=B**2/rho+a**2
    cfast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=sum(wroe(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,mag(:))**2,dim=ndim+1)/wroe(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,rho_)+csound2(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3)

    ! cs**2=0.5*(aa+dsqrt(aa**2-4*a**2*(b_i**2/rho)))
    cslow(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half*(cfast(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3)-dsqrt(cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3)**2-4d0*csound2(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       mag(idim))**2/wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)))

    ! cf**2=aa-cs**2
    cfast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=cfast(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3)-cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3)

    ! alpha_f**2=(a**2-cs**2)/(cf**2-cs**2)
    afast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(csound2(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3)-cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3))/(cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3)-cslow(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
    afast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=min(one,&
       max(afast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3),zero))

    ! alpha_s=dsqrt(1-alpha_f**2)
    aslow(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dsqrt(one-&
       afast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))

    ! alpha_f=dsqrt(alpha_f**2)
    afast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dsqrt(afast(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3))

    ! cf=dsqrt(cf**2)
    cfast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dsqrt(cfast(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3))

    ! cs=dsqrt(cs**2)
    cslow(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dsqrt(cslow(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3))

    !Replace the primitive variables with more useful quantities:
    ! rho -> dsqrt(rho)
    wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       rho_)=dsqrt(wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))

    ! Avoid sgn(b_idim)==0
    where(dabs(wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       mag(idim)))<smalldouble)wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       mag(idim))=smalldouble
    ! B_idir,jdir -> beta_idir,jdir
    idir=idim+1-ndir*(idim/ndir)
    if(ndir==2)then
       where(wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))>=zero)
          wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))=one
       elsewhere
          wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))=-one
       end where
    else
       !beta_j=B_j/dsqrt(B_i**2+B_j**2); beta_i=B_i/dsqrt(B_i**2+B_j**2)
       jdir=idir+1-ndir*(idir/ndir)
       tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=dsqrt(wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))**2+wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(jdir))**2)
       where(tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)>smalldouble)
          wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(idir))=wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(idir))/tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)
          wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(jdir))=wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(jdir))/tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)
       elsewhere
          wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(idir))=dsqrt(half)
          wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(jdir))=dsqrt(half)
       end where
    endif

  end subroutine average2

  ! Calculate the il-th characteristic speed and the jump in the il-th
  ! characteristic variable in the idim direction within ixL.
  ! The eigenvalues and the l=r**(-1) matrix is calculated from wroe.
  ! jump(il)=Sum_il l(il,iw)*(wR(iw)-wL(iw)), where w are the conservative
  ! variables. However part of the summation is done in advance and saved into
  ! bdv,bdb,dp and dv variables. "smalla" contains a lower limit for "a" to be
  ! used in the entropy fix.
  !
  ! All the variables are centered on the cell interface, thus the 
  ! "*C" notation is omitted for sake of brevity.
  subroutine mhd_get_eigenjump(wL,wR,wroe,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
     ixmax3,il,idim,smalla,a,jump,workroe)
    use mod_global_parameters

    integer, intent(in)                         :: ixmin1,ixmin2,ixmin3,ixmax1,&
       ixmax2,ixmax3,il,idim
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nw)       :: wL,wR,wroe
    double precision, intent(in)                :: x(ixGlo1:ixGhi1,&
       ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:ndim)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)          :: smalla,a,jump
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nworkroe) :: workroe

    call geteigenjump2(wL,wR,wroe,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
       il,idim,smalla,a,jump, workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1),workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,2),&
        workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,3),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,4),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,5),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,6),&
        workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,7),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,8),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,9),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,10),&
        workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,11),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,12),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,13))

  end subroutine mhd_get_eigenjump

  ! Calculate the il-th characteristic speed and the jump in the il-th
  ! characteristic variable in the idim direction within ixL.
  ! The eigenvalues and the l=r**(-1) matrix is calculated from wroe.
  ! jump(il)=Sum_il l(il,iw)*(wR(iw)-wL(iw)), where w are the conservative
  ! variables. However part of the summation is done in advance and saved into
  ! bdv,bdb,dp and dv variables. "smalla" contains a lower limit for "a" to be
  ! used in the entropy fix.
  !
  ! All the variables are centered on the cell interface, thus the 
  ! "*C" notation is omitted for sake of brevity.
  subroutine geteigenjump2(wL,wR,wroe,x,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
     ixmax3,il,idim,smalla,a,jump, cfast,cslow,afast,aslow,csound2,dp,rhodv,&
     bdv,bdb,cs2L,cs2R,cs2ca2L,cs2ca2R)
    use mod_global_parameters
    use mod_tvd

    integer                               :: ixmin1,ixmin2,ixmin3,ixmax1,&
       ixmax2,ixmax3,il,idim,idir,jdir
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       nw) :: wL,wR,wroe
    double precision, intent(in)          :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,1:ndim)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)    :: smalla,a,jump
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)    :: cfast,cslow,afast,aslow,csound2,dp,rhodv
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)    :: bdv,bdb
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)    :: aL,aR,cs2L,cs2R,cs2ca2L,cs2ca2R

    idir=idim+1-ndir*(idim/ndir)
    jdir=idir+1-ndir*(idir/ndir)

    if(il==fastRW_)then
       !Fast and slow waves use bdv=sqrho**2*sign(bx)*(betay*dvy+betaz*dvz)
       !                        bdb=sqrho*a*          (betay*dBy+betaz*dBz)
       bdv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))* (wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idir))/wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,mom(idir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_))
       if(ndir==3)bdv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)=bdv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)+wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mag(jdir))* (wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mom(jdir))/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mom(jdir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
       bdv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=bdv(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)*sign(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_)**2,wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,mag(idim)))

       bdb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))*(wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))-wL(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idir)))
       if(ndir==3)bdb(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)=bdb(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)+wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mag(jdir))*(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mag(jdir))-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mag(jdir)))
       bdb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=bdb(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)*dsqrt(csound2(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3))*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_)
    endif

    if(il==alfvRW_)then
       !Alfven waves use      bdv=0.5*sqrho**2*      (betaz*dvy-betay*dvz)
       !                      bdb=0.5*sqrho*sign(bx)*(betaz*dBy-betay*dBz)
       bdv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(jdir))* (wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idir))/wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,mom(idir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_)) -wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mag(idir))* (wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mom(jdir))/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
          mom(jdir))/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
       bdb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(jdir))*(wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))-wL(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))) -wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))*(wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(jdir))-wL(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(jdir)))
       bdv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=bdv(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)*half*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_)**2
       bdb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=bdb(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)*half*sign(wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_),wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,mag(idim)))
    endif

    select case(il)
    case(fastRW_)
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))+cfast(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half/csound2(&
          ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)*(afast(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)*(+cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*rhodv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)+dp(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))+aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*(-cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*bdv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)+bdb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)))
    case(fastLW_)
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))-cfast(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half/csound2(&
          ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)*(afast(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)*(-cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*rhodv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)+dp(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))+aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*(+cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*bdv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)+bdb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)))
    case(slowRW_)
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))+cslow(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half/csound2(&
          ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)*(aslow(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)*(+cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*rhodv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)+dp(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))+afast(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*(+cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*bdv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)-bdb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)))
    case(slowLW_)
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))-cslow(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half/csound2(&
          ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)*(aslow(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)*(-cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*rhodv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)+dp(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3))+afast(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*(-cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)*bdv(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)-bdb(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)))
    case(entroW_)
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3,rho_)-dp(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)/csound2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)
    case(diverW_)
       if(divbwave)then
          a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))
          jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mag(idim))-wL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3,mag(idim))
       else
          a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
          jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
       endif
    case(alfvRW_)
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))+dabs(wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idim)))/wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=+bdv(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)-bdb(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)
    case(alfvLW_)
       a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))-dabs(wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mag(idim)))/wroe(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)
       jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-bdv(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3)-bdb(ixmin1:ixmax1,ixmin2:ixmax2,&
          ixmin3:ixmax3)
    end select

    ! Calculate "smalla" or modify "a" based on the "typeentropy" switch

    select case(typeentropy(il))
    case('yee')
       ! Based on Yee JCP 68,151 eq 3.23
       smalla(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=entropycoef(il)
    case('harten','powell', 'ratio')
       ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
       ! Initialize left and right eigenvalues by velocities
       aL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= wL(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wL(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)
       aR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,mom(idim))/wR(ixmin1:ixmax1,&
          ixmin2:ixmax2,ixmin3:ixmax3,rho_)
       ! Calculate the final "aL" and "aR"
       select case(il)
       case(fastRW_)
          ! These quantities will be used for all the fast and slow waves
          ! Calculate soundspeed**2 and cs**2+ca**2.
          call mhd_get_csound2(wL,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,cs2L)
          call mhd_get_csound2(wR,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,cs2R)
          cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=cs2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)+sum(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(:))**2,dim=ndim+1)/wL(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,rho_)
          cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=cs2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)+sum(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(:))**2,dim=ndim+1)/wR(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,rho_)
          ! Save the discriminants into cs2L and cs2R
          cs2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=dsqrt(cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)**2-4d0*cs2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(idim))**2/wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))
          cs2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=dsqrt(cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)**2-4d0*cs2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(idim))**2/wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_))

          ! The left and right eigenvalues for the fast wave going to right
          aL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + dsqrt(half*(cs2ca2L(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + cs2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)))
          aR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + dsqrt(half*(cs2ca2R(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + cs2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)))
       case(fastLW_)
          aL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - dsqrt(half*(cs2ca2L(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + cs2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)))
          aR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - dsqrt(half*(cs2ca2R(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + cs2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)))
       case(slowRW_)
          aL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + dsqrt(half*(cs2ca2L(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - cs2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)))
          aR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + dsqrt(half*(cs2ca2R(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - cs2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)))
       case(slowLW_)
          aL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - dsqrt(half*(cs2ca2L(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - cs2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)))
          aR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - dsqrt(half*(cs2ca2R(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - cs2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)))
       case(entroW_,diverW_)
          ! These propagate by the velocity
       case(alfvRW_)
          ! Store the Alfven speeds into cs2ca2L and cs2ca2R
          cs2ca2L(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=dabs(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(idim)))/dsqrt(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))
          cs2ca2R(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)=dabs(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(idim)))/dsqrt(wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             rho_))

          aL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + cs2ca2L(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
          aR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) + cs2ca2R(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
       case(alfvLW_)
          aL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aL(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - cs2ca2L(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
          aR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=aR(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3) - cs2ca2R(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
       end select
    end select

    call entropyfix(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,il,aL,aR,a,&
       smalla)

  end subroutine geteigenjump2

  ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe
  subroutine mhd_rtimes(q,w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,il,&
     idim,rq,workroe)
    use mod_global_parameters

    integer, intent(in)             :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3, iw, il, idim
    double precision, intent(in)    :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nw), q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision, intent(inout) :: rq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nworkroe)

    call rtimes2(q,w,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,il,idim,rq,&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,1),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,2),&
        workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,3),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,4),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,5),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,6),&
        workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,7),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,14),&
       workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,15))

  end subroutine mhd_rtimes

  ! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe
  subroutine rtimes2(q,wroe,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,iw,il,&
     idim,rq, cfast,cslow,afast,aslow,csound2,dp,rhodv,bv,v2a2)
    use mod_global_parameters

    integer                            :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3,iw,il,idim,idir,jdir
    double precision                   :: wroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3,nw)
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: q,rq
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: cfast,cslow,afast,aslow,csound2,dp,rhodv
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3) :: bv,v2a2

    idir=idim+1-ndir*(idim/ndir)
    jdir=idir+1-ndir*(idir/ndir)

    if(iw == rho_) then
      select case(il)
      case(fastRW_,fastLW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*afast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)
      case(slowRW_,slowLW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)
      case(entroW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)
      case(diverW_,alfvRW_,alfvLW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
      end select
    else if(iw == e_) then
      if(il==fastRW_)then
        ! Store 0.5*v**2+(2-gamma)/(gamma-1)*a**2
        v2a2(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)=half*sum(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3,mom(:))**2,dim=ndim+1)+ &
           (two-mhd_gamma)/(mhd_gamma-one)*csound2(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)
        ! Store sgn(bx)*(betay*vy+betaz*vz) in bv
        bv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wroe(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))*wroe(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3,mom(idir))
        if(ndir==3)bv(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)=bv(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           mag(jdir))*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           mom(jdir))
        bv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=bv(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*sign(one,wroe(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3,mag(idim)))
      else if(il==alfvRW_)then
        !Store betaz*vy-betay*vz in bv
        bv(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(wroe(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3,mag(jdir))*wroe(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3,mom(idir))-wroe(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3,mag(idir))*wroe(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3,mom(jdir)))
      endif

      select case(il)
      case(fastRW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*(-aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*bv(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+afast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*(v2a2(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*(cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           mom(idim)))))
      case(fastLW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*(+aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*bv(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+afast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*(v2a2(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*(cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)-wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           mom(idim)))))
      case(slowRW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*(+afast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*bv(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*(v2a2(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*(cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           mom(idim)))))
      case(slowLW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*(-afast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*bv(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*(v2a2(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)+cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)*(cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)-wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
           mom(idim)))))
      case(entroW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*half*sum(wroe(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3,mom(:))**2,dim=ndim+1)
      case(diverW_)
        if(divbwave)then
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mag(idim))
        else
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= zero
        endif
      case(alfvRW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=+q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*bv(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)
      case(alfvLW_)
        rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-q(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*bv(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3)
      end select
    else if(any(mom(:)==iw)) then
      if(iw==mom(idim))then
        select case(il)
        case(fastRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*afast(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)+cfast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
        case(fastLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*afast(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)-cfast(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
        case(slowRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)+cslow(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
        case(slowLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*(wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)-cslow(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
        case(entroW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,iw)
        case(diverW_,alfvLW_,alfvRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
        end select
      else
        select case(il)
        case(fastRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*(afast(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)-aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(1)-mom(1)+iw)*sign(one,wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mag(idim))))
        case(fastLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*(afast(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)+aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*cslow(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(1)-mom(1)+iw)*sign(one,wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mag(idim))))
        case(slowRW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*(aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)+afast(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(1)-mom(1)+iw)*sign(one,wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mag(idim))))
        case(slowLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*(aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)-afast(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*cfast(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             mag(1)-mom(1)+iw)*sign(one,wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,mag(idim))))
        case(entroW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3,iw)
        case(diverW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
        case(alfvRW_)
          if(iw==mom(idir))then
            rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=+q(ixmin1:ixmax1,&
               ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,mag(jdir))
          else
            rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-q(ixmin1:ixmax1,&
               ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,mag(idir))
          endif
        case(alfvLW_)
          if(iw==mom(idir))then
            rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-q(ixmin1:ixmax1,&
               ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,mag(jdir))
          else
            rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=+q(ixmin1:ixmax1,&
               ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,mag(idir))
          endif
        end select
      end if ! iw=m_idir,m_jdir
    else if(any(mag(:)==iw)) then
      if(iw==mag(idim))then
        if(il==diverW_ .and. divbwave)then
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)
        else
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
        endif
      else
        select case(il)
        case(fastRW_,fastLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=+q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*aslow(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*dsqrt(csound2(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3))*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)/wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)
        case(slowRW_,slowLW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-q(ixmin1:ixmax1,&
             ixmin2:ixmax2,ixmin3:ixmax3)*afast(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3)*dsqrt(csound2(ixmin1:ixmax1,ixmin2:ixmax2,&
             ixmin3:ixmax3))*wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
             iw)/wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)
        case(entroW_,diverW_)
          rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=zero
        case(alfvRW_,alfvLW_)
          if(iw==mag(idir))then
            rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=-q(ixmin1:ixmax1,&
               ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,mag(jdir))/sign(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,rho_),wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,mag(idim)))
          else
            rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=+q(ixmin1:ixmax1,&
               ixmin2:ixmax2,ixmin3:ixmax3)*wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,mag(idir))/sign(wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,rho_),wroe(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,mag(idim)))
          end if
        end select
      end if ! iw=b_idir,b_jdir
    end if

  end subroutine rtimes2

end module mod_mhd_roe
