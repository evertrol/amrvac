!> To get a RBSL magnetic flux rope in 3D (Titov 2018 ApJL 852, L21)
module mod_rbsl
  implicit none

contains


  subroutine rbsl(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,np,a,F_flx,positive_helicity,x,&
     x_axis,Atotal,Bfr)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    ! resolution of flux rope axis
    integer, intent(in)             :: np
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! coordinates of flux rope axis (integral path)
    double precision, intent(in)    :: x_axis(np,1:ndim)
    ! cross-sectional radius of flux rope
    double precision, intent(in)    :: a
    ! net magnetic flux along flux rope axis
    double precision, intent(in)    :: F_flx
    ! is positive helicity
    logical, intent(in) :: positive_helicity
    ! vector potential of flux rope
    double precision, intent(out)   :: Atotal(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! magnetic field of flux rope
    double precision, optional, intent(out)   :: Bfr(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)

    ! net current along flux rope axis for azimuthal magnetic field
    double precision   :: I_cur
    ! vector potential for azimuthal magnetic field
    double precision   :: AIx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    ! vector potential for axial magnetic field
    double precision   :: AFx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    integer :: ix1,ix2,ix3, ixp, idirmin
    double precision :: r_mag, KIr, KFr, dl, re_pi, sqrt1r, f52r,fsqrt6
    double precision :: Rpl(1:ndim), r_vec(1:ndim), Rcr(1:ndim)

    if(positive_helicity) then
      I_cur = 5.d0*sqrt(2.d0)*F_flx/(3.d0*4.d0*dpi*a)
    else
      I_cur =-5.d0*sqrt(2.d0)*F_flx/(3.d0*4.d0*dpi*a)
    end if

    re_pi=1.d0/dpi
    AIx = 0.d0
    AFx = 0.d0
    fsqrt6=1.d0/sqrt(6.d0)
    do ix3=ixOmin3,ixOmax3
do ix2=ixOmin2,ixOmax2
do ix1=ixOmin1,ixOmax1
      do ixp=1,np
        ! position vector from source point to field point
        r_vec(:) = (x(ix1,ix2,ix3,:) - x_axis(ixp,:))/a
        r_mag = sqrt(sum(r_vec(:)**2))
        ! calculate tangential vector of the axis which is a circle
        if (ixp == 1) then
          Rpl(:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(np,:))
        else if (ixp==np) then
          Rpl(:) = 0.5d0*(x_axis(1,:)-x_axis(ixp-1,:))
        else
          Rpl(:) = 0.5d0*(x_axis(ixp+1,:)-x_axis(ixp-1,:))
        end if
        ! Rpl X r_vec
        Rcr(1) = Rpl(2)*r_vec(3) - Rpl(3)*r_vec(2)
        Rcr(2) = Rpl(3)*r_vec(1) - Rpl(1)*r_vec(3)
        Rcr(3) = Rpl(1)*r_vec(2) - Rpl(2)*r_vec(1)
        if (r_mag < 1.d0) then
          sqrt1r=sqrt(1.d0-r_mag**2)
          f52r=5.d0-2.d0*r_mag**2
          KIr = 2.d0*re_pi*(asin(r_mag)/r_mag + f52r*third*sqrt1r)
          KFr = 2.d0*re_pi/r_mag**2*(asin(r_mag)/r_mag-sqrt1r) + &
             2.d0*re_pi*sqrt1r + f52r*0.5d0*fsqrt6*(1.d0 - &
             2.d0*re_pi*asin((1.d0+2.d0*r_mag**2)/f52r))
        else
          KIr = 1.d0/r_mag
          KFr = KIr**3
        endif
        AIx(ix1,ix2,ix3,:) = AIx(ix1,ix2,ix3,:) + KIr*Rpl(:)
        AFx(ix1,ix2,ix3,:) = AFx(ix1,ix2,ix3,:) + KFr*Rcr(:)
      end do
      AIx(ix1,ix2,ix3,:) = AIx(ix1,ix2,ix3,:)*I_cur/a
      AFx(ix1,ix2,ix3,:) = AFx(ix1,ix2,ix3,:)*F_flx*0.25d0*re_pi/a**2
    end do
end do
end do
    Atotal=AIx+AFx
    if(present(Bfr)) call curlvector(Atotal,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Bfr,&
       idirmin,1,ndir)

  end subroutine rbsl


end module mod_rbsl
