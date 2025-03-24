!> Module containing the MP5 (fifth order) flux scheme
module mod_mp5

  implicit none
  private

  public :: MP5limiter
  public :: MP5limitervar
  public :: MP5limiterL
  public :: MP5limiterR

contains

  !> MP5 limiter from Suresh & Huynh 1997 Following the convention of Mignone et
  !> al. 2010. Needs at least three ghost cells
  subroutine MP5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
     iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC,wRC)
    use mod_global_parameters
    use mod_weno, only: fix_limiter1

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) 
    ! .. local ..
    integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
       iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
       iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
       iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,iLpppmin3,iLpppmax1,&
       iLpppmax2,iLpppmax3
    integer                         :: idmin1,idmin2,idmin3,idmax1,idmax2,&
       idmax3, idpmin1,idpmin2,idpmin3,idpmax1,idpmax2,idpmax3, idppmin1,&
       idppmin2,idppmin3,idppmax1,idppmax2,idppmax3, idmmin1,idmmin2,idmmin3,&
       idmmax1,idmmax2,idmmax3, iemin1,iemin2,iemin3,iemax1,iemax2,iemax3,&
        iemmin1,iemmin2,iemmin3,iemmax1,iemmax2,iemmax3, iepmin1,iepmin2,&
       iepmin3,iepmax1,iepmax2,iepmax3, ieppmin1,ieppmin2,ieppmin3,ieppmax1,&
       ieppmax2,ieppmax3
    integer                         :: iw
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc,&
        flim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)  :: wRCtmp, wLCtmp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.d0, alpha=4.0d0
    !double precision                :: alpha
    !----------------------------------------------------------------------------

    ! Variable alpha:
    !alpha = float(nstep)/courantpar - one

    ! Left side:
    ! range to process:
    !iLmin^D=ixmin^D-kr(idims,^D);iLmax^D=ixmax^D;

    ! HALL
    ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
    ! also add one ghost zone!
    !   {iL^L=iL^L^LADD1;}

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
    iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
    iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
    iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
    iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);

    f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLmmmin1:iLmmmax1,&
       iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3,&
       1:nwflux) - 13.0d0* w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
       1:nwflux) + 47.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux) - 3.0d0*  w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
          iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw)-w(iLmin1:iLmax1,iLmin2:iLmax2,&
          iLmin3:iLmax3,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)-w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
          iLmmin3:iLmmax3,iw))
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
       fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
          iLmin3:iLmax3)
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
          iLmin3:iLmax3)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1;idmax2=iLmax2;idmax3=iLmax3; idmin1=iLmin1-kr(idims,1)
    idmin2=iLmin2-kr(idims,2);idmin3=iLmin3-kr(idims,3);
    idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
    idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
    idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
    idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
    idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
    idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

    iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2)
    iemax3=idmax3+kr(idims,3); iemin1=idmin1;iemin2=idmin2;iemin3=idmin3;
    iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
    iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
    iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
    iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
    iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
    iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);

    d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3,1:nwflux) = w(iepmin1:iepmax1,&
       iepmin2:iepmax2,iepmin3:iepmax3,1:nwflux)-2.0d0*w(iemin1:iemax1,&
       iemin2:iemax2,iemin3:iemax3,1:nwflux)+w(iemmin1:iemmax1,iemmin2:iemmax2,&
       iemmin3:iemmax3,1:nwflux)

    do iw=1,nwflux
       a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3,iw)-d(idpmin1:idpmax1,idpmin2:idpmax2,&
          idpmin3:idpmax3,iw)
       b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idpmin1:idpmax1,&
          idpmin2:idpmax2,idpmin3:idpmax3,iw)-d(idmin1:idmax1,idmin2:idmax2,&
          idmin3:idmax3,iw)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,a,b,tmp)
       a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3,iw)
       b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idpmin1:idpmax1,&
          idpmin2:idpmax2,idpmin3:idpmax3,iw)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3,iw) = tmp3(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3,1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3,1:nwflux)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = half*(3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
       1:nwflux)) + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       iLmmin3:iLmmax3,1:nwflux)

    fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = max(min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
       fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
       min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

    fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = min(max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
       fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
       max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

    do iw=1,nwflux
       a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmin(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       c(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmax(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,a,b,c,tmp)
       flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = tmp(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)-w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)) .le. eps)
       wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
          1:nwflux) = f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
    elsewhere
       wLCtmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
          1:nwflux) = flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
    end where

    ! Right side:
    ! the interpolation from the right is obtained when the left-hand formula is applied to
    ! data mirrored about the interface.  
    ! thus substitute: 
    ! i-2 -> i+3
    ! i-1 -> i+2
    ! i   -> i+1
    ! i+1 -> i
    ! i+2 -> i-1

    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
    iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

    f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3,&
       1:nwflux) - 13.0d0* w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3,1:nwflux) + 47.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3,1:nwflux) + 27.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3,1:nwflux) - 3.0d0*  w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       iLmmin3:iLmmax3,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
          iLpmin3:iLpmax3,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLpmin1:iLpmax1,&
          iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw)-w(iLppmin1:iLppmax1,&
          iLppmin2:iLppmax2,iLppmin3:iLppmax3,iw))
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
       fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLpmin1:iLpmax1,&
          iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw) + tmp(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLpmin1:iLpmax1,&
          iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
          iLmin3:iLmax3)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2)
    idmax3=iLmax3+kr(idims,3); idmin1=iLmin1;idmin2=iLmin2;idmin3=iLmin3;
    idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
    idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
    idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
    idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
    idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
    idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

    iemax1=idmax1;iemax2=idmax2;iemax3=idmax3; iemin1=idmin1-kr(idims,1)
    iemin2=idmin2-kr(idims,2);iemin3=idmin3-kr(idims,3);
    iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
    iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
    iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
    iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
    iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
    iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);
    ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
    ieppmin3=iepmin3+kr(idims,3);ieppmax1=iepmax1+kr(idims,1)
    ieppmax2=iepmax2+kr(idims,2);ieppmax3=iepmax3+kr(idims,3);

    d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3,1:nwflux) = w(iemin1:iemax1,&
       iemin2:iemax2,iemin3:iemax3,1:nwflux)-2.0d0*w(iepmin1:iepmax1,&
       iepmin2:iepmax2,iepmin3:iepmax3,1:nwflux)+w(ieppmin1:ieppmax1,&
       ieppmin2:ieppmax2,ieppmin3:ieppmax3,1:nwflux)

    do iw=1,nwflux
       a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3,iw)-d(idmmin1:idmmax1,idmmin2:idmmax2,&
          idmmin3:idmmax3,iw)
       b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmmin1:idmmax1,&
          idmmin2:idmmax2,idmmin3:idmmax3,iw)-d(idmin1:idmax1,idmin2:idmax2,&
          idmin3:idmax3,iw)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,a,b,tmp)
       a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3,iw)
       b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmmin1:idmmax1,&
          idmmin2:idmmax2,idmmin3:idmmax3,iw)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3,iw) = tmp3(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3,1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3,1:nwflux)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = half*(3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3,1:nwflux) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3,1:nwflux)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux)

    fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = max(min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux),w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
       min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

    fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = min(max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux),w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
       max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

    do iw=1,nwflux
       a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmin(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       c(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmax(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,a,b,c,tmp)
       flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = tmp(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux))  .le. eps)
       wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
          1:nwflux) = f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
    elsewhere
       wRCtmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
          1:nwflux) = flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
    end where

    call fix_limiter1(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
       iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,wLCtmp,wRCtmp,wLC,wRC)

  end subroutine MP5limiter

  subroutine MP5limiterL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC)
    use mod_global_parameters
    !use mod_weno, only: fix_onelimiter1

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) 
    ! .. local ..
    integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
       iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
       iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
       iLppmax1,iLppmax2,iLppmax3
    integer                         :: idmin1,idmin2,idmin3,idmax1,idmax2,&
       idmax3, idpmin1,idpmin2,idpmin3,idpmax1,idpmax2,idpmax3, idppmin1,&
       idppmin2,idppmin3,idppmax1,idppmax2,idppmax3, idmmin1,idmmin2,idmmin3,&
       idmmax1,idmmax2,idmmax3, iemin1,iemin2,iemin3,iemax1,iemax2,iemax3,&
        iemmin1,iemmin2,iemmin3,iemmax1,iemmax2,iemmax3, iepmin1,iepmin2,&
       iepmin3,iepmax1,iepmax2,iepmax3, ieppmin1,ieppmin2,ieppmin3,ieppmax1,&
       ieppmax2,ieppmax3
    integer                         :: iw
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc,&
        flim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.d0, alpha=4.0d0
    !double precision                :: alpha
    !double precision, dimension(ixI^S,1:nw)  :: wLCtmp

    ! Variable alpha:
    !alpha = float(nstep)/courantpar - one

    ! Left side:
    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
    iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
    iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
    iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
    iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);

    f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLmmmin1:iLmmmax1,&
       iLmmmin2:iLmmmax2,iLmmmin3:iLmmmax3,&
       1:nwflux) - 13.0d0* w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
       1:nwflux) + 47.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) + 27.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux) - 3.0d0*  w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
          iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw)-w(iLmin1:iLmax1,iLmin2:iLmax2,&
          iLmin3:iLmax3,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)-w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
          iLmmin3:iLmmax3,iw))
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
       fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
          iLmin3:iLmax3)
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
          iLmin3:iLmax3)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1;idmax2=iLmax2;idmax3=iLmax3; idmin1=iLmin1-kr(idims,1)
    idmin2=iLmin2-kr(idims,2);idmin3=iLmin3-kr(idims,3);
    idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
    idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
    idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
    idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
    idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
    idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

    iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2)
    iemax3=idmax3+kr(idims,3); iemin1=idmin1;iemin2=idmin2;iemin3=idmin3;
    iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
    iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
    iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
    iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
    iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
    iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);

    d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3,1:nwflux) = w(iepmin1:iepmax1,&
       iepmin2:iepmax2,iepmin3:iepmax3,1:nwflux)-2.0d0*w(iemin1:iemax1,&
       iemin2:iemax2,iemin3:iemax3,1:nwflux)+w(iemmin1:iemmax1,iemmin2:iemmax2,&
       iemmin3:iemmax3,1:nwflux)

    do iw=1,nwflux
       a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3,iw)-d(idpmin1:idpmax1,idpmin2:idpmax2,&
          idpmin3:idpmax3,iw)
       b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idpmin1:idpmax1,&
          idpmin2:idpmax2,idpmin3:idpmax3,iw)-d(idmin1:idmax1,idmin2:idmax2,&
          idmin3:idmax3,iw)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,a,b,tmp)
       a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3,iw)
       b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idpmin1:idpmax1,&
          idpmin2:idpmax2,idpmin3:idpmax3,iw)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3,iw) = tmp3(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3,1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3,1:nwflux)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = half*(3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,iLmmin3:iLmmax3,&
       1:nwflux)) + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       iLmmin3:iLmmax3,1:nwflux)

    fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = max(min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
       fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
       min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

    fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = min(max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
       fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
       max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

    do iw=1,nwflux
       a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmin(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       c(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmax(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,a,b,c,tmp)
       flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = tmp(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)-w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)) .le. eps)
       wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
          1:nwflux) = f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
    elsewhere
       wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
          1:nwflux) = flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
    end where

    !call fix_onelimiter1(ixI^L,iL^L,wLCtmp,wLC)

  end subroutine MP5limiterL

  subroutine MP5limiterR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wRC)
    use mod_global_parameters
    !use mod_weno, only: fix_onelimiter1

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    ! .. local ..
    integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
       iLmmax3, iLpmin1,iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,&
       iLppmin2,iLppmin3,iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,&
       iLpppmin3,iLpppmax1,iLpppmax2,iLpppmax3
    integer                         :: idmin1,idmin2,idmin3,idmax1,idmax2,&
       idmax3, idpmin1,idpmin2,idpmin3,idpmax1,idpmax2,idpmax3, idppmin1,&
       idppmin2,idppmin3,idppmax1,idppmax2,idppmax3, idmmin1,idmmin2,idmmin3,&
       idmmax1,idmmax2,idmmax3, iemin1,iemin2,iemin3,iemax1,iemax2,iemax3,&
        iemmin1,iemmin2,iemmin3,iemmax1,iemmax2,iemmax3, iepmin1,iepmin2,&
       iepmin3,iepmax1,iepmax2,iepmax3, ieppmin1,ieppmin2,ieppmin3,ieppmax1,&
       ieppmax2,ieppmax3
    integer                         :: iw
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc,&
        flim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.d0, alpha=4.0d0
    !double precision                :: alpha
    !double precision, dimension(ixI^S,1:nw)  :: wRCtmp

    ! Right side:
    ! the interpolation from the right is obtained when the left-hand formula is applied to
    ! data mirrored about the interface.  
    ! thus substitute: 
    ! i-2 -> i+3
    ! i-1 -> i+2
    ! i   -> i+1
    ! i+1 -> i
    ! i+2 -> i-1

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
    iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
    iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
    iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);
    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
    iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

    f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = 1.0d0/60.0d0 * (2.0d0* w(iLpppmin1:iLpppmax1,&
       iLpppmin2:iLpppmax2,iLpppmin3:iLpppmax3,&
       1:nwflux) - 13.0d0* w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3,1:nwflux) + 47.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3,1:nwflux) + 27.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3,1:nwflux) - 3.0d0*  w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       iLmmin3:iLmmax3,1:nwflux))

    ! get fmp and ful:
    do iw=1,nwflux
       a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
          iLpmin3:iLpmax3,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLpmin1:iLpmax1,&
          iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw)-w(iLppmin1:iLppmax1,&
          iLppmin2:iLppmax2,iLppmin3:iLppmax3,iw))
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
       fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLpmin1:iLpmax1,&
          iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw) + tmp(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = w(iLpmin1:iLpmax1,&
          iLpmin2:iLpmax2,iLpmin3:iLpmax3,iw) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
          iLmin3:iLmax3)
    end do ! iw loop

    ! get dm4:
    idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2)
    idmax3=iLmax3+kr(idims,3); idmin1=iLmin1;idmin2=iLmin2;idmin3=iLmin3;
    idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
    idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
    idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
    idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
    idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
    idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

    iemax1=idmax1;iemax2=idmax2;iemax3=idmax3; iemin1=idmin1-kr(idims,1)
    iemin2=idmin2-kr(idims,2);iemin3=idmin3-kr(idims,3);
    iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
    iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
    iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
    iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
    iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
    iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);
    ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
    ieppmin3=iepmin3+kr(idims,3);ieppmax1=iepmax1+kr(idims,1)
    ieppmax2=iepmax2+kr(idims,2);ieppmax3=iepmax3+kr(idims,3);

    d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3,1:nwflux) = w(iemin1:iemax1,&
       iemin2:iemax2,iemin3:iemax3,1:nwflux)-2.0d0*w(iepmin1:iepmax1,&
       iepmin2:iepmax2,iepmin3:iepmax3,1:nwflux)+w(ieppmin1:ieppmax1,&
       ieppmin2:ieppmax2,ieppmin3:ieppmax3,1:nwflux)

    do iw=1,nwflux
       a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3,iw)-d(idmmin1:idmmax1,idmmin2:idmmax2,&
          idmmin3:idmmax3,iw)
       b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmmin1:idmmax1,&
          idmmin2:idmmax2,idmmin3:idmmax3,iw)-d(idmin1:idmax1,idmin2:idmax2,&
          idmin3:idmax3,iw)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,a,b,tmp)
       a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3,iw)
       b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmmin1:idmmax1,&
          idmmin2:idmmax2,idmmin3:idmmax3,iw)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
       call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,&
          idmin2,idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
       dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3,iw) = tmp3(idmin1:idmax1,&
          idmin2:idmax2,idmin3:idmax3)
    end do

    ! get fmd:
    fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3,1:nwflux))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3,1:nwflux)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = half*(3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3,1:nwflux) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3,1:nwflux)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux)

    fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = max(min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux),w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
       min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

    fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux) = min(max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux),w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)),&
       max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,1:nwflux),&
       ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)))

    do iw=1,nwflux
       a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmin(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       c(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = fmax(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3,iw)
       call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,a,b,c,tmp)
       flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw) = tmp(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
    end do

    ! check case
    where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3,&
       1:nwflux))*(f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
       1:nwflux))  .le. eps)
       wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
          1:nwflux) = f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
    elsewhere
       wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,&
          1:nwflux) = flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,1:nwflux)
    end where
  
    !call fix_onelimiter1(ixI^L,iL^L,wRCtmp,wRC)

  end subroutine MP5limiterR

  subroutine minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,minm)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(out):: minm(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    minm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = (sign(one,&
       a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))+sign(one,&
       b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))/2.0d0 * min(abs(a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)),abs(b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))

  end subroutine minmod

  subroutine median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,c,med)

    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: a(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        c(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(out):: med(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision             :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = c(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) - a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    med(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) + (sign(one,tmp1(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3))+sign(one,tmp2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)))/2.0d0 * min(abs(tmp1(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)),abs(tmp2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)))

  end subroutine median

  !> MP5 limiter from Suresh & Huynh 1997
  !> Following the convention of Mignone et al. 2010.
  !> Needs at least three ghost cells.  Set nghostcells=3.
  ! for one variable only: no fixing applied
  subroutine MP5limitervar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC,wRC)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(inout) :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),wLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) 

    integer                         :: iLmmin1,iLmmin2,iLmmin3,iLmmax1,iLmmax2,&
       iLmmax3, iLmmmin1,iLmmmin2,iLmmmin3,iLmmmax1,iLmmmax2,iLmmmax3, iLpmin1,&
       iLpmin2,iLpmin3,iLpmax1,iLpmax2,iLpmax3, iLppmin1,iLppmin2,iLppmin3,&
       iLppmax1,iLppmax2,iLppmax3, iLpppmin1,iLpppmin2,iLpppmin3,iLpppmax1,&
       iLpppmax2,iLpppmax3
    integer                         :: idmin1,idmin2,idmin3,idmax1,idmax2,&
       idmax3, idpmin1,idpmin2,idpmin3,idpmax1,idpmax2,idpmax3, idppmin1,&
       idppmin2,idppmin3,idppmax1,idppmax2,idppmax3, idmmin1,idmmin2,idmmin3,&
       idmmax1,idmmax2,idmmax3, iemin1,iemin2,iemin3,iemax1,iemax2,iemax3,&
        iemmin1,iemmin2,iemmin3,iemmax1,iemmax2,iemmax3, iepmin1,iepmin2,&
       iepmin3,iepmax1,iepmax2,iepmax3, ieppmin1,ieppmin2,ieppmin3,ieppmax1,&
       ieppmax2,ieppmax3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)  :: f, fmp, fmin, fmax, ful, dm4, d, fmd, flc, flim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)  :: wRCtmp, wLCtmp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: tmp, tmp2, tmp3, a, b, c
    double precision, parameter     :: eps=0.0d0, alpha=4.0d0
    !double precision                :: alpha

    ! Variable alpha:
    !alpha = float(nstep)/courantpar - one

    ! Left side:
    ! range to process:
    !iLmin^D=ixmin^D-kr(idims,^D);iLmax^D=ixmax^D;

    ! HALL
    ! For Hall, we need one more reconstructed layer since currents are computed in getflux:
    ! also add one ghost zone!
    !   {iL^L=iL^L^LADD1;}

    ! iL^L holds the indices of interfaces to reconstruct to.  Convention is that a center index holds the _right-side_ interface.  

    iLmmin1=iLmin1-kr(idims,1);iLmmin2=iLmin2-kr(idims,2)
    iLmmin3=iLmin3-kr(idims,3);iLmmax1=iLmax1-kr(idims,1)
    iLmmax2=iLmax2-kr(idims,2);iLmmax3=iLmax3-kr(idims,3);
    iLmmmin1=iLmmin1-kr(idims,1);iLmmmin2=iLmmin2-kr(idims,2)
    iLmmmin3=iLmmin3-kr(idims,3);iLmmmax1=iLmmax1-kr(idims,1)
    iLmmmax2=iLmmax2-kr(idims,2);iLmmmax3=iLmmax3-kr(idims,3);
    iLpmin1=iLmin1+kr(idims,1);iLpmin2=iLmin2+kr(idims,2)
    iLpmin3=iLmin3+kr(idims,3);iLpmax1=iLmax1+kr(idims,1)
    iLpmax2=iLmax2+kr(idims,2);iLpmax3=iLmax3+kr(idims,3);
    iLppmin1=iLpmin1+kr(idims,1);iLppmin2=iLpmin2+kr(idims,2)
    iLppmin3=iLpmin3+kr(idims,3);iLppmax1=iLpmax1+kr(idims,1)
    iLppmax2=iLpmax2+kr(idims,2);iLppmax3=iLpmax3+kr(idims,3);

    f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 1.0d0/60.0d0 * (2.0d0* &
       w(iLmmmin1:iLmmmax1,iLmmmin2:iLmmmax2,&
       iLmmmin3:iLmmmax3) - 13.0d0* w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       iLmmin3:iLmmax3) + 47.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3) + 27.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3) - 3.0d0*  w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3))

    ! get fmp and ful:
    a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,iLpmin3:iLpmax3)-w(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)
    b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3)-w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       iLmmin3:iLmmax3))
    call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
       iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
    fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)
    ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)

    ! get dm4:
    idmax1=iLmax1;idmax2=iLmax2;idmax3=iLmax3; idmin1=iLmin1-kr(idims,1)
    idmin2=iLmin2-kr(idims,2);idmin3=iLmin3-kr(idims,3);
    idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
    idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
    idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
    idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
    idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
    idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

    iemax1=idmax1+kr(idims,1);iemax2=idmax2+kr(idims,2)
    iemax3=idmax3+kr(idims,3); iemin1=idmin1;iemin2=idmin2;iemin3=idmin3;
    iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
    iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
    iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
    iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
    iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
    iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);

    d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3) = w(iepmin1:iepmax1,&
       iepmin2:iepmax2,iepmin3:iepmax3)-2.0d0*w(iemin1:iemax1,iemin2:iemax2,&
       iemin3:iemax3)+w(iemmin1:iemmax1,iemmin2:iemmax2,iemmin3:iemmax3)

    a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
       idmin2:idmax2,idmin3:idmax3)-d(idpmin1:idpmax1,idpmin2:idpmax2,&
       idpmin3:idpmax3)
    b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idpmin1:idpmax1,&
       idpmin2:idpmax2,idpmin3:idpmax3)-d(idmin1:idmax1,idmin2:idmax2,&
       idmin3:idmax3)
    call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
       idmin3,idmax1,idmax2,idmax3,a,b,tmp)
    a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
       idmin2:idmax2,idmin3:idmax3)
    b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idpmin1:idpmax1,&
       idpmin2:idpmax2,idpmin3:idpmax3)
    call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
       idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
    call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
       idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
    dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = tmp3(idmin1:idmax1,&
       idmin2:idmax2,idmin3:idmax3)

    ! get fmd:
    fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = &
       half*(3.0d0*w(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3) - w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       iLmmin3:iLmmax3)) + 4.0d0/3.0d0*dm4(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       iLmmin3:iLmmax3)

    fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = max(min(w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3),fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)),&
       min(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3),ful(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3),flc(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)))

    fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = min(max(w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3),w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3),fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)),&
       max(w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3),ful(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3),flc(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)))

    call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
       iLmin3,iLmax1,iLmax2,iLmax3,fmin,f,fmax,tmp)
    flim(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = tmp(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3)

    ! check case
    where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)-w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3))*(f(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)) .le. eps)
       wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
    elsewhere
       wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flim(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
    end where

    ! Right side:
    ! the interpolation from the right is obtained when the left-hand formula is applied to
    ! data mirrored about the interface.  
    ! thus substitute: 
    ! i-2 -> i+3
    ! i-1 -> i+2
    ! i   -> i+1
    ! i+1 -> i
    ! i+2 -> i-1

    iLpppmin1=iLppmin1+kr(idims,1);iLpppmin2=iLppmin2+kr(idims,2)
    iLpppmin3=iLppmin3+kr(idims,3);iLpppmax1=iLppmax1+kr(idims,1)
    iLpppmax2=iLppmax2+kr(idims,2);iLpppmax3=iLppmax3+kr(idims,3);

    f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = 1.0d0/60.0d0 * (2.0d0* &
       w(iLpppmin1:iLpppmax1,iLpppmin2:iLpppmax2,&
       iLpppmin3:iLpppmax3) - 13.0d0* w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3) + 47.0d0* w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3) + 27.0d0* w(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3) - 3.0d0*  w(iLmmin1:iLmmax1,iLmmin2:iLmmax2,&
       iLmmin3:iLmmax3))

    ! get fmp and ful:
    a(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3)-w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3)
    b(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = alpha*(w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,iLpmin3:iLpmax3)-w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3))
    call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
       iLmin3,iLmax1,iLmax2,iLmax3,a,b,tmp)
    fmp(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,iLpmin3:iLpmax3) + tmp(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)
    ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,iLpmin3:iLpmax3) + b(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)

    ! get dm4:
    idmax1=iLmax1+kr(idims,1);idmax2=iLmax2+kr(idims,2)
    idmax3=iLmax3+kr(idims,3); idmin1=iLmin1;idmin2=iLmin2;idmin3=iLmin3;
    idmmin1=idmin1-kr(idims,1);idmmin2=idmin2-kr(idims,2)
    idmmin3=idmin3-kr(idims,3);idmmax1=idmax1-kr(idims,1)
    idmmax2=idmax2-kr(idims,2);idmmax3=idmax3-kr(idims,3);
    idpmin1=idmin1+kr(idims,1);idpmin2=idmin2+kr(idims,2)
    idpmin3=idmin3+kr(idims,3);idpmax1=idmax1+kr(idims,1)
    idpmax2=idmax2+kr(idims,2);idpmax3=idmax3+kr(idims,3);

    iemax1=idmax1;iemax2=idmax2;iemax3=idmax3; iemin1=idmin1-kr(idims,1)
    iemin2=idmin2-kr(idims,2);iemin3=idmin3-kr(idims,3);
    iemmin1=iemin1-kr(idims,1);iemmin2=iemin2-kr(idims,2)
    iemmin3=iemin3-kr(idims,3);iemmax1=iemax1-kr(idims,1)
    iemmax2=iemax2-kr(idims,2);iemmax3=iemax3-kr(idims,3);
    iepmin1=iemin1+kr(idims,1);iepmin2=iemin2+kr(idims,2)
    iepmin3=iemin3+kr(idims,3);iepmax1=iemax1+kr(idims,1)
    iepmax2=iemax2+kr(idims,2);iepmax3=iemax3+kr(idims,3);
    ieppmin1=iepmin1+kr(idims,1);ieppmin2=iepmin2+kr(idims,2)
    ieppmin3=iepmin3+kr(idims,3);ieppmax1=iepmax1+kr(idims,1)
    ieppmax2=iepmax2+kr(idims,2);ieppmax3=iepmax3+kr(idims,3);

    d(iemin1:iemax1,iemin2:iemax2,iemin3:iemax3) = w(iemin1:iemax1,&
       iemin2:iemax2,iemin3:iemax3)-2.0d0*w(iepmin1:iepmax1,iepmin2:iepmax2,&
       iepmin3:iepmax3)+w(ieppmin1:ieppmax1,ieppmin2:ieppmax2,&
       ieppmin3:ieppmax3)

    a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmin1:idmax1,&
       idmin2:idmax2,idmin3:idmax3)-d(idmmin1:idmmax1,idmmin2:idmmax2,&
       idmmin3:idmmax3)
    b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = 4.0d0*d(idmmin1:idmmax1,&
       idmmin2:idmmax2,idmmin3:idmmax3)-d(idmin1:idmax1,idmin2:idmax2,&
       idmin3:idmax3)
    call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
       idmin3,idmax1,idmax2,idmax3,a,b,tmp)
    a(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmin1:idmax1,&
       idmin2:idmax2,idmin3:idmax3)
    b(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = d(idmmin1:idmmax1,&
       idmmin2:idmmax2,idmmin3:idmmax3)
    call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
       idmin3,idmax1,idmax2,idmax3,a,b,tmp2)
    call minmod(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,idmin1,idmin2,&
       idmin3,idmax1,idmax2,idmax3,tmp,tmp2,tmp3)
    dm4(idmin1:idmax1,idmin2:idmax2,idmin3:idmax3) = tmp3(idmin1:idmax1,&
       idmin2:idmax2,idmin3:idmax3)

    ! get fmd:
    fmd(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = (w(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3)+w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3))/2.0d0-dm4(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)/2.0d0

    !get flc: 
    flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = &
       half*(3.0d0*w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3) - w(iLppmin1:iLppmax1,iLppmin2:iLppmax2,&
       iLppmin3:iLppmax3)) + 4.0d0/3.0d0*dm4(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3)

    fmin(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = &
       max(min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3),&
       w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3),fmd(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3)),min(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3),ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)))

    fmax(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = &
       min(max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,iLpmin3:iLpmax3),&
       w(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3),fmd(iLmin1:iLmax1,&
       iLmin2:iLmax2,iLmin3:iLmax3)),max(w(iLpmin1:iLpmax1,iLpmin2:iLpmax2,&
       iLpmin3:iLpmax3),ful(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3),&
       flc(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)))

    call median(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,iLmin2,&
       iLmin3,iLmax1,iLmax2,iLmax3,fmin,f,fmax,flim)

    ! check case
    where ((f(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3)-w(iLpmin1:iLpmax1,&
       iLpmin2:iLpmax2,iLpmin3:iLpmax3))*(f(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3)-fmp(iLmin1:iLmax1,iLmin2:iLmax2,&
       iLmin3:iLmax3))  .le. eps)
       wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = f(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
    elsewhere
       wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3) = flim(iLmin1:iLmax1,&
          iLmin2:iLmax2,iLmin3:iLmax3)
    end where

  end subroutine MP5limitervar

end module mod_mp5
