!> Module containing Roe solver for scalar advection
module mod_rho_roe
  use mod_physics_roe
  use mod_rho_phys

  implicit none
  private

  public :: rho_roe_init

contains

  subroutine rho_roe_init()
    use mod_physics_roe

    nworkroe = 1

    phys_average         => rho_average
    phys_get_eigenjump   => rho_get_eigenjump
    phys_rtimes          => rho_rtimes
  end subroutine rho_roe_init

  subroutine rho_average(wL, wR, x, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,&
      idim, wroe, workroe)
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

    wroe(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        rho_)=half*(wL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
        rho_)+wR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3, rho_))
  end subroutine rho_average

  subroutine rho_get_eigenjump(wL, wR, wC, x, ixmin1,ixmin2,ixmin3,ixmax1,&
     ixmax2,ixmax3, il, idim, smalla, a, jump, workroe)

    ! Calculate the characteristic speed a and the jump in the
    ! characteristic variable in the idim direction within ixL.
    ! For a scalar equation the characteristic and conservative variables coincide
    ! The characteristic speed is just the velocity, but it should be averaged
    ! for the cell interfaces

    use mod_global_parameters

    integer, intent(in)                          :: ixmin1,ixmin2,ixmin3,&
       ixmax1,ixmax2,ixmax3, il, idim
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
        nw)       :: wL, wR, wC
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)           :: smalla, a, jump, v
    double precision, dimension(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
        nworkroe) :: workroe
    double precision, intent(in)                 :: x(ixGlo1:ixGhi1,&
       ixGlo2:ixGhi2,ixGlo3:ixGhi3, 1:3)
    integer                                      :: jxmin1,jxmin2,jxmin3,&
       jxmax1,jxmax2,jxmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3

    jxmin1=ixmin1+kr(idim,1);jxmin2=ixmin2+kr(idim,2)
    jxmin3=ixmin3+kr(idim,3);jxmax1=ixmax1+kr(idim,1)
    jxmax2=ixmax2+kr(idim,2);jxmax3=ixmax3+kr(idim,3);
    ixCmin1=ixmin1;ixCmin2=ixmin2;ixCmin3=ixmin3; ixCmax1=jxmax1
    ixCmax2=jxmax2;ixCmax3=jxmax3;

    ! No entropy fix
    smalla(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)= -one
    ! The velocity is independent of w in the transport equation,
    ! but it may depend on the location
    call rho_get_v(wL, x, ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idim, v, .true.)

    a(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(v(jxmin1:jxmax1,&
       jxmin2:jxmax2,jxmin3:jxmax3)+v(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3))/2

    jump(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=wR(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3, rho_)-wL(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3, rho_)

  end subroutine rho_get_eigenjump

  subroutine rho_rtimes(q, w, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, iw,&
      il, idim, rq, workroe)

    ! Multiply q by R(il, iw), where R is the right eigenvalue matrix at wC.
    ! For a scalar equation the R matrix is unity

    use mod_global_parameters
    integer, intent(in)             :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
       ixmax3, iw, il, idim
    double precision, intent(in)    :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nw), q(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
    double precision, intent(inout) :: rq(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3)
    double precision, intent(inout) :: workroe(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3, nworkroe)

    rq(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=q(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3)

  end subroutine rho_rtimes

end module mod_rho_roe
