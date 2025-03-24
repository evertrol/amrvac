module mod_venk
  ! Venkatakrishnan limiter
  !
  ! 2019.10.11 coded up by nanami;
  !
  ! see Venkatakrishnan 1993 for this limiter;


  implicit none
  private

  public :: venklimiter

contains

  subroutine venklimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,dxdim,w,wLC,wRC)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
    double precision, intent(in)    :: dxdim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) 
    !> local
    integer                         :: iMmin1,iMmin2,iMmin3,iMmax1,iMmax2,&
       iMmax3, iMmmin1,iMmmin2,iMmmin3,iMmmax1,iMmmax2,iMmmax3, iMpmin1,&
       iMpmin2,iMpmin3,iMpmax1,iMpmax2,iMpmax3
    integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
       iLmmax3, iLpmin1,iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,&
       iLppmin2,iLppmin3,iLppmax1,iLppmax2,iLppmax3
    double precision                :: wmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),wmin(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: westp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),westm(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: phi1(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),phi2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: phi3(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),phi4(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: phim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),phip(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),phi(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: deltap(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),deltam(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: eps2
    double precision, parameter     :: venk_omega = 1.d-12
    double precision, parameter     :: venk_k = 0.3d0

    iMmin1=iLmin1;iMmin2=iLmin2;iMmin3=iLmin3;iMmax1=iLmax1;iMmax2=iLmax2
    iMmax3=iLmax3;
    iMmax1=iLmax1+kr(idims,1);iMmax2=iLmax2+kr(idims,2)
    iMmax3=iLmax3+kr(idims,3);
    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
    iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
    iMmmin1=iMmin1-kr(idims,1);iMmmin2=iMmin2-kr(idims,2)
    iMmmin3=iMmin3-kr(idims,3);iMmmax1=iMmax1-kr(idims,1)
    iMmmax2=iMmax2-kr(idims,2);iMmmax3=iMmax3-kr(idims,3);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
    iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
    iMpmin1=iMmin1+kr(idims,1);iMpmin2=iMmin2+kr(idims,2)
    iMpmin3=iMmin3+kr(idims,3);iMpmax1=iMmax1+kr(idims,1)
    iMpmax2=iMmax2+kr(idims,2);iMpmax3=iMmax3+kr(idims,3);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
    iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);

    eps2 = (venk_k * dxdim) ** 3

    wmax(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) = max(w(iMmmin1:iMmmax1,iMmmin2:iMmmax2,iMmmin3:iMmmax3,&
       1:nwflux), w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux),&
        w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,iMpmin3:iMpmax3,1:nwflux))
    wmin(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) = min(w(iMmmin1:iMmmax1,iMmmin2:iMmmax2,iMmmin3:iMmmax3,&
       1:nwflux), w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux),&
        w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,iMpmin3:iMpmax3,1:nwflux))
    !> use central difference approximation as (eq.5) and take phi = 1
    westp(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) = w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) + (w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,iMpmin3:iMpmax3,&
       1:nwflux) - w(iMmmin1:iMmmax1,iMmmin2:iMmmax2,iMmmin3:iMmmax3,&
       1:nwflux)) * 0.25d0
    westm(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) = w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) - (w(iMpmin1:iMpmax1,iMpmin2:iMpmax2,iMpmin3:iMpmax3,&
       1:nwflux) - w(iMmmin1:iMmmax1,iMmmin2:iMmmax2,iMmmin3:iMmmax3,&
       1:nwflux)) * 0.25d0

    !> (eq.30) & (eq.31)
    deltap = 0
    deltam = 0
    phi1 = 0
    phi2 = 0
    phip = 1
    where(westp(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) .gt. w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = wmax(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) - w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux)
      deltam(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = westp(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) - w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux)
      deltam(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = sign(dabs(deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)) + venk_omega, deltam(iMmin1:iMmax1,&
         iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      phi1(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = (deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux)**2 + eps2) + 2.d0 * deltap(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)
      phi2(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux)**2 + 2.d0 * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)**2 + deltap(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) + eps2
      phip(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = phi1(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) / phi2(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3, 1:nwflux)
    elsewhere(westp(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) .lt. w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = wmin(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) - w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux)
      deltam(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = westp(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) - w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux)
      deltam(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = sign(dabs(deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)) + venk_omega, deltam(iMmin1:iMmax1,&
         iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      phi1(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = (deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux)**2 + eps2) + 2.d0 * deltap(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)
      phi2(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux)**2 + 2.d0 * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)**2 + deltap(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) + eps2
      phip(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = phi1(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) / phi2(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3, 1:nwflux)
    elsewhere
      phip(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux) = one
   endwhere

    deltap = 0
    deltam = 0
    phi3 = 0
    phi4 = 0
    phim = 0
    where(westm(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) .lt. w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = - (wmax(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) - w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      deltam(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = westm(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) - w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux)
      deltam(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = sign(dabs(deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)) + venk_omega, deltam(iMmin1:iMmax1,&
         iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      phi3(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = (deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux)**2 + eps2) + 2.d0 * deltap(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)
      phi4(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux)**2 + 2.d0 * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)**2 + deltap(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) + eps2
      phim(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = phi3(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) / phi4(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3, 1:nwflux)
    elsewhere(westm(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) .gt. w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = - (wmin(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) - w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      deltam(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = westm(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) - w(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux)
      deltam(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = sign(dabs(deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)) + venk_omega, deltam(iMmin1:iMmax1,&
         iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
      phi3(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = (deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux)**2 + eps2) + 2.d0 * deltap(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)
      phi4(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = deltap(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux)**2 + 2.d0 * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux)**2 + deltap(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) * deltam(iMmin1:iMmax1,iMmin2:iMmax2,&
         iMmin3:iMmax3,1:nwflux) + eps2
      phim(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) = phi3(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
         1:nwflux) / phi4(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3, 1:nwflux)
    elsewhere
      phim(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux) = one
    endwhere
    !> (eq.3)
    phi(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux) = min(phim(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,&
       1:nwflux),phip(iMmin1:iMmax1,iMmin2:iMmax2,iMmin3:iMmax3,1:nwflux))
    !> (eq.5)
    wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) + 0.25d0 * (w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux)-w(iLmmin1:iLmmax1,&
       iLmmin2:iLmmax2,iLmmin3:iLmmax3,1:nwflux)) * phi(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
    wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux) - 0.25d0 * (w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3,1:nwflux)-w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)) * phi(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux)

  end subroutine venklimiter

end module mod_venk
