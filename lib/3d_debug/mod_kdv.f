!> Module for including kdv source term in simulations
!> adds \f$-\delta^2*\sum_i \partial_{iii} \rho \f$ over dimensions i
module mod_kdv
  implicit none

  !> source split or not
  logical :: kdv_split= .false.
  !> forefactor \f$ \delta^2\f$  of \f$ \partial_{iii} \f$ term
  double precision :: kdv_delta = 1.0d0
  !> switch for second order [1] or fourth order [2] central FD for \f$ \partial_{iii}\f$
  !> Note: fourth order needs 3 nghostcells, all assume equidistant grid
  integer :: kdv_order = 1

contains
  !> Read this module's parameters from a file
  subroutine kdv_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /kdv_list/ kdv_split, kdv_delta, kdv_order

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, kdv_list, end=111)
111    close(unitpar)
    end do

  end subroutine kdv_params_read

  !> Initialize the module
  subroutine kdv_init()
    use mod_global_parameters

    call kdv_params_read(par_files)

  end subroutine kdv_init

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine kdv_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
     qsourcesplit,active)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active

    integer          :: idir, lxmin1,lxmin2,lxmin3,lxmax1,lxmax2,lxmax3,kxmin1,&
       kxmin2,kxmin3,kxmax1,kxmax2,kxmax3,jxmin1,jxmin2,jxmin3,jxmax1,jxmax2,&
       jxmax3,hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3,gxmin1,gxmin2,gxmin3,&
       gxmax1,gxmax2,gxmax3,fxmin1,fxmin2,fxmin3,fxmax1,fxmax2,fxmax3
    double precision :: skdv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    if(qsourcesplit .eqv. kdv_split) then
      active = .true.
      skdv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
      select case(kdv_order)
        case(1)
          do idir=1,ndim
             ! The source is based on the time centered wCT
             kxmin1=ixOmin1+2*kr(idir,1);kxmin2=ixOmin2+2*kr(idir,2)
             kxmin3=ixOmin3+2*kr(idir,3);kxmax1=ixOmax1+2*kr(idir,1)
             kxmax2=ixOmax2+2*kr(idir,2);kxmax3=ixOmax3+2*kr(idir,3);
             jxmin1=ixOmin1+kr(idir,1);jxmin2=ixOmin2+kr(idir,2)
             jxmin3=ixOmin3+kr(idir,3);jxmax1=ixOmax1+kr(idir,1)
             jxmax2=ixOmax2+kr(idir,2);jxmax3=ixOmax3+kr(idir,3);
             hxmin1=ixOmin1-kr(idir,1);hxmin2=ixOmin2-kr(idir,2)
             hxmin3=ixOmin3-kr(idir,3);hxmax1=ixOmax1-kr(idir,1)
             hxmax2=ixOmax2-kr(idir,2);hxmax3=ixOmax3-kr(idir,3);
             gxmin1=ixOmin1-2*kr(idir,1);gxmin2=ixOmin2-2*kr(idir,2)
             gxmin3=ixOmin3-2*kr(idir,3);gxmax1=ixOmax1-2*kr(idir,1)
             gxmax2=ixOmax2-2*kr(idir,2);gxmax3=ixOmax3-2*kr(idir,3);
             ! 2nd order centered difference for -\partial_xxx \rho 
             ! warning: needs 2 ghostcells, equidistant grid
             skdv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=skdv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)+(wCT(kxmin1:kxmax1,kxmin2:kxmax2,&
                kxmin3:kxmax3,iw_rho)-2.0d0*wCT(jxmin1:jxmax1,jxmin2:jxmax2,&
                jxmin3:jxmax3,iw_rho) +2.0d0*wCT(hxmin1:hxmax1,hxmin2:hxmax2,&
                hxmin3:hxmax3,iw_rho)-wCT(gxmin1:gxmax1,gxmin2:gxmax2,&
                gxmin3:gxmax3,iw_rho)) /(2.0d0 *dxlevel(idir)**3)
          enddo
        case(2)
          do idir=1,ndim
             ! The source is based on the time centered wCT
             lxmin1=ixOmin1+3*kr(idir,1);lxmin2=ixOmin2+3*kr(idir,2)
             lxmin3=ixOmin3+3*kr(idir,3);lxmax1=ixOmax1+3*kr(idir,1)
             lxmax2=ixOmax2+3*kr(idir,2);lxmax3=ixOmax3+3*kr(idir,3);
             kxmin1=ixOmin1+2*kr(idir,1);kxmin2=ixOmin2+2*kr(idir,2)
             kxmin3=ixOmin3+2*kr(idir,3);kxmax1=ixOmax1+2*kr(idir,1)
             kxmax2=ixOmax2+2*kr(idir,2);kxmax3=ixOmax3+2*kr(idir,3);
             jxmin1=ixOmin1+kr(idir,1);jxmin2=ixOmin2+kr(idir,2)
             jxmin3=ixOmin3+kr(idir,3);jxmax1=ixOmax1+kr(idir,1)
             jxmax2=ixOmax2+kr(idir,2);jxmax3=ixOmax3+kr(idir,3);
             hxmin1=ixOmin1-kr(idir,1);hxmin2=ixOmin2-kr(idir,2)
             hxmin3=ixOmin3-kr(idir,3);hxmax1=ixOmax1-kr(idir,1)
             hxmax2=ixOmax2-kr(idir,2);hxmax3=ixOmax3-kr(idir,3);
             gxmin1=ixOmin1-2*kr(idir,1);gxmin2=ixOmin2-2*kr(idir,2)
             gxmin3=ixOmin3-2*kr(idir,3);gxmax1=ixOmax1-2*kr(idir,1)
             gxmax2=ixOmax2-2*kr(idir,2);gxmax3=ixOmax3-2*kr(idir,3);
             fxmin1=ixOmin1-3*kr(idir,1);fxmin2=ixOmin2-3*kr(idir,2)
             fxmin3=ixOmin3-3*kr(idir,3);fxmax1=ixOmax1-3*kr(idir,1)
             fxmax2=ixOmax2-3*kr(idir,2);fxmax3=ixOmax3-3*kr(idir,3);
             ! 4th order centered difference for -\partial_xxx \rho 
             ! warning: needs 3 ghostcells, equidistant grid
             skdv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)=skdv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)+(-wCT(lxmin1:lxmax1,lxmin2:lxmax2,&
                lxmin3:lxmax3,iw_rho)+8.0d0*wCT(kxmin1:kxmax1,kxmin2:kxmax2,&
                kxmin3:kxmax3,iw_rho)-13.0d0*wCT(jxmin1:jxmax1,jxmin2:jxmax2,&
                jxmin3:jxmax3,iw_rho) +13.0d0*wCT(hxmin1:hxmax1,hxmin2:hxmax2,&
                hxmin3:hxmax3,iw_rho)-8.0d0*wCT(gxmin1:gxmax1,gxmin2:gxmax2,&
                gxmin3:gxmax3,iw_rho)+wCT(fxmin1:fxmax1,fxmin2:fxmax2,&
                fxmin3:fxmax3,iw_rho)) /(8.0d0 *dxlevel(idir)**3)
          enddo
        case default
          call mpistop('undefined kdv_order parameter: see mod_kdv.t')
      end select
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_rho) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_rho) - qdt*kdv_delta**2*skdv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    end if

  end subroutine kdv_add_source

  subroutine kdv_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim), w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: dxarr(ndim), max_sinefactor

    !> Time step constraint for leap-frog time stepping combined with central 2nd order FD 
    !> see e.g. Chun-Te Lee et al, Journal of Mathematics Research vol. 9, no.4, 2017
    !> ISSN 1916-9795
    dxarr(1)=dx1;dxarr(2)=dx2;dxarr(3)=dx3;
    max_sinefactor=3.0d0*dsqrt(3.0d0)/2.0d0
    dtnew=dtdiffpar*minval(dxarr(1:ndim))**3/(max_sinefactor*kdv_delta**2)

  end subroutine kdv_get_dt

end module mod_kdv
