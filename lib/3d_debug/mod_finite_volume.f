!> Module with finite volume methods for fluxes
module mod_finite_volume
#include "amrvac.h"
  implicit none
  private

  public :: finite_volume
  public :: hancock
  public :: reconstruct_LR

contains

  !> The non-conservative Hancock predictor for TVDLF
  !>
  !> on entry:
  !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
  !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
  !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
  subroutine hancock(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimsmin,idimsmax,&
     qtC,sCT,qt,snew,dxs,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idimsmin,idimsmax
    double precision, intent(in) :: qdt, dtfactor,qtC, qt, dxs(ndim),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    type(state) :: sCT, snew

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wprim, wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wLp, wRp
    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)      :: inv_volume
    double precision :: fLC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        nwflux), fRC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, nwflux)
    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3
    logical :: active

    associate(wCT=>sCT%w,wnew=>snew%w)
    ! Expand limits in each idims direction in which fluxes are added
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmin3=ixOmin3;ixmax1=ixOmax1
    ixmax2=ixOmax2;ixmax3=ixOmax3;
    do idims= idimsmin,idimsmax
       ixmin1=ixmin1-kr(idims,1);ixmin2=ixmin2-kr(idims,2)
       ixmin3=ixmin3-kr(idims,3);ixmax1=ixmax1+kr(idims,1)
       ixmax2=ixmax2+kr(idims,2);ixmax3=ixmax3+kr(idims,3);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImin3>ixmin3.or.ixImax1<ixmax1.or.&
       ixImax2<ixmax2.or.ixImax3<ixmax3) call &
       mpistop("Error in Hancock: Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wprim,x)

    dxinv=-qdt/dxs
    if(.not.slab_uniform) inv_volume = 1.d0/block%dvolume(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    do idims= idimsmin,idimsmax
      b0i=idims
      ! Calculate w_j+g_j/2 and w_j-g_j/2
      ! First copy all variables, then upwind wLC and wRC.
      ! wLC is to the left of ixO, wRC is to the right of wCT.
      hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
      hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
      hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);

      wRp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
         1:nwflux)=wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)=wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)

      ! apply limited reconstruction for left and right status at cell interfaces
      call reconstruct_LR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,hxOmin1,hxOmin2,&
         hxOmin3,hxOmax1,hxOmax2,hxOmax3,idims,wprim,wLC,wRC,wLp,wRp,x,&
         dxs(idims))

      ! Calculate the fLC and fRC fluxes
      call phys_get_flux(wRC,wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,idims,fRC)
      call phys_get_flux(wLC,wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,fLC)

      ! Advect w(iw)
      if (slab_uniform) then
        if(local_timestep) then
          do iw=1,nwflux
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw)-block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)*dtfactor/dxs(idims)* (fLC(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, iw)-fRC(hxOmin1:hxOmax1,&
               hxOmin2:hxOmax2,hxOmin3:hxOmax3, iw))
          end do
        else  
          do iw=1,nwflux
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw)+dxinv(idims)* (fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, iw)-fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3, iw))
          end do
        endif
      else
        if(local_timestep) then
          do iw=1,nwflux
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw) - block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)*dtfactor * inv_volume &
               *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,idims)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, iw) -block%surfaceC(hxOmin1:hxOmax1,&
               hxOmin2:hxOmax2,hxOmin3:hxOmax3,idims)*fRC(hxOmin1:hxOmax1,&
               hxOmin2:hxOmax2,hxOmin3:hxOmax3, iw))
          end do
        else
          do iw=1,nwflux
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw) - qdt * inv_volume *(block%surfaceC(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,idims)*fLC(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw) -block%surfaceC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3,idims)*fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3, iw))
          end do
        end if
      end if
    end do ! next idims
    b0i=0

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,dtfactor,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3,wCT,wnew,x)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
        dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim),ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       1,nw,qtC,wCT,wprim,qt,wnew,x,.false.,active)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImin3,&
          ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,'exit hancock finite_volume')
    endif
    end associate
  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimsmin,&
     idimsmax, qtC,sCT,qt,snew,fC,fE,dxs,x)
    use mod_physics
    use mod_variables
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods
    use mod_comm_lib, only: mpistop

    integer, intent(in)                                   :: method
    double precision, intent(in)                          :: qdt, dtfactor,&
        qtC, qt, dxs(ndim)
    integer, intent(in)                                   :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3, idimsmin,idimsmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), intent(in) :: x
    type(state)                                           :: sCT, snew
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux,1:ndim)    :: fC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)             :: fE

    ! primitive w at cell center
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wLp, wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux) :: fLC, fRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)      :: cmaxC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)      :: cminC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)      :: Hspeed
    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv
    integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)               :: patchf
    integer :: idims, iw, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ixCmin1,ixCmin2,ixCmin3,&
       ixCmax1,ixCmax2,ixCmax3, ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,&
       ixCRmax3, kxCmin1,kxCmin2,kxCmin3,kxCmax1,kxCmax2,kxCmax3, kxRmin1,&
       kxRmin2,kxRmin3,kxRmax1,kxRmax2,kxRmax3, ii
    logical :: active
    type(ct_velocity) :: vcts

    associate(wCT=>sCT%w, wnew=>snew%w)

    fC=0.d0
    fLC=0.d0
    fRC=0.d0

    ! The flux calculation contracts by one in the idims direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmin3=ixOmin3;ixmax1=ixOmax1
    ixmax2=ixOmax2;ixmax3=ixOmax3;
    do idims= idimsmin,idimsmax
       ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
       ixmin3=ixmin3-2*kr(idims,3);ixmax1=ixmax1+2*kr(idims,1)
       ixmax2=ixmax2+2*kr(idims,2);ixmax3=ixmax3+2*kr(idims,3);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImin3>ixmin3.or.ixImax1<ixmax1.or.&
       ixImax2<ixmax2.or.ixImax3<ixmax3) call &
       mpistop("Error in fv : Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,wprim,x)

    do idims= idimsmin,idimsmax
       ! use interface value of w0 at idims
       b0i=idims

       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
       hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);

       kxCmin1=ixImin1;kxCmin2=ixImin2;kxCmin3=ixImin3
       kxCmax1=ixImax1-kr(idims,1);kxCmax2=ixImax2-kr(idims,2)
       kxCmax3=ixImax3-kr(idims,3);
       kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
       kxRmin3=kxCmin3+kr(idims,3);kxRmax1=kxCmax1+kr(idims,1)
       kxRmax2=kxCmax2+kr(idims,2);kxRmax3=kxCmax3+kr(idims,3);

       if(stagger_grid) then
         ! ct needs all transverse cells
         ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
         ixCmax2=ixOmax2+nghostcells-nghostcells*kr(idims,2)
         ixCmax3=ixOmax3+nghostcells-nghostcells*kr(idims,3);
         ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1)
         ixCmin2=hxOmin2-nghostcells+nghostcells*kr(idims,2)
         ixCmin3=hxOmin3-nghostcells+nghostcells*kr(idims,3);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=hxOmin1
         ixCmin2=hxOmin2;ixCmin3=hxOmin3;
       end if

       ! wRp and wLp are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idims they are
       ! shifted half a cell towards the 'lower' direction.
       wRp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,&
          1:nw)=wprim(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,1:nw)
       wLp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,&
          1:nw)=wprim(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,1:nw)

       ! Determine stencil size
       ixCRmin1 = max(ixCmin1 - phys_wider_stencil,ixGlo1)
       ixCRmin2 = max(ixCmin2 - phys_wider_stencil,ixGlo2)
       ixCRmin3 = max(ixCmin3 - phys_wider_stencil,ixGlo3)
       ixCRmax1 = min(ixCmax1 + phys_wider_stencil,ixGhi1)
       ixCRmax2 = min(ixCmax2 + phys_wider_stencil,ixGhi2)
       ixCRmax3 = min(ixCmax3 + phys_wider_stencil,ixGhi3)

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,ixCRmin1,&
          ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,idims,wprim,wLC,wRC,wLp,&
          wRp,x,dxs(idims))

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,qt,wLC,wRC,wLp,&
          wRp,sCT,idims)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(wLC,wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fLC)
       call phys_get_flux(wRC,wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,fRC)

       if(H_correction) then
         call phys_get_H_speed(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
            ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,&
            Hspeed)
       end if
       ! estimating bounds for the minimum and maximum signal velocities
       if(method==fs_tvdlf.or.method==fs_tvdmu) then
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImin3,&
            ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
            ixCmax3,idims,Hspeed,cmaxC)
         ! index of var  velocity appears in the induction eq. 
         if(stagger_grid) call phys_get_ct_velocity(vcts,wLp,wRp,ixImin1,&
            ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,ixCmin3,&
            ixCmax1,ixCmax2,ixCmax3,idims,cmaxC(ixImin1:ixImax1,&
            ixImin2:ixImax2,ixImin3:ixImax3,index_v_mag))
       else
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImin3,&
            ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
            ixCmax3,idims,Hspeed,cmaxC,cminC)
         if(stagger_grid) call phys_get_ct_velocity(vcts,wLp,wRp,ixImin1,&
            ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,ixCmin3,&
            ixCmax1,ixCmax2,ixCmax3,idims,cmaxC(ixImin1:ixImax1,&
            ixImin2:ixImax2,ixImin3:ixImax3,index_v_mag),cminC(ixImin1:ixImax1,&
            ixImin2:ixImax2,ixImin3:ixImax3,index_v_mag))
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case(fs_hll)
         do ii=1,number_species
           call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
         end do
       case(fs_hllc,fs_hllcd)
         do ii=1,number_species
           call get_Riemann_flux_hllc(start_indices(ii),stop_indices(ii))
         end do
       case(fs_hlld)
         do ii=1,number_species
           if(ii==index_v_mag) then
             call get_Riemann_flux_hlld(start_indices(ii),stop_indices(ii))
           else
             call get_Riemann_flux_hll(start_indices(ii),stop_indices(ii))
           endif   
         end do
       case(fs_tvdlf)
         do ii=1,number_species
           call get_Riemann_flux_tvdlf(start_indices(ii),stop_indices(ii))
         end do
       case(fs_tvdmu)
         call get_Riemann_flux_tvdmu()
       case default
         call mpistop('unkown Riemann flux in finite volume')
       end select

    end do ! Next idims
    b0i=0
    if(stagger_grid) call phys_update_faces(ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,&
       wprim,fC,fE,sCT,snew,vcts)
    if(slab_uniform) then
      dxinv=-qdt/dxs
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
        hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
        ! TODO maybe put if outside loop idims: but too much code is copy pasted
        ! this is also done in hancock and fd, centdiff in mod_finite_difference
        if(local_timestep) then
          do iw=iwstart,nwflux
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)=-block%dt(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3)*dtfactor/dxs(idims)*fC(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImin3:ixImax3,iw,idims)
          end do
        else
          ! Multiply the fluxes by -dt/dx since Flux fixing expects this
          fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nwflux,&
             idims)=dxinv(idims)*fC(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,1:nwflux,idims)
        end if

        wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iwstart:nwflux)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iwstart:nwflux)+(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iwstart:nwflux,idims)-fC(hxOmin1:hxOmax1,&
           hxOmin2:hxOmax2,hxOmin3:hxOmax3,iwstart:nwflux,idims))

        ! For the MUSCL scheme apply the characteristic based limiter
        if(method==fs_tvdmu) call tvdlimit2(method,qdt,ixImin1,ixImin2,ixImin3,&
           ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
           ixCmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,wLC,&
           wRC,wnew,x,fC,dxs)

      end do ! Next idims
    else
      inv_volume = 1.d0/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      do idims= idimsmin,idimsmax
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
        hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);

        if(local_timestep) then
          do iw=iwstart,nwflux
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)=-block%dt(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3)*dtfactor*fC(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3,iw,idims)*block%surfaceC(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImin3:ixImax3,idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
               idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,&
               idims)) * inv_volume
          end do
        else
          do iw=iwstart,nwflux
            fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw,&
               idims)=-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw,idims)*block%surfaceC(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3,idims)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw,&
               idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,iw,&
               idims)) * inv_volume
          end do
        end if 
        ! For the MUSCL scheme apply the characteristic based limiter
        if (method==fs_tvdmu) call tvdlimit2(method,qdt,ixImin1,ixImin2,&
           ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
           ixCmax2,ixCmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           idims,wLC,wRC,wnew,x,fC,dxs)

      end do ! Next idims
    end if

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,dtfactor,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3,wCT,wnew,x)

    if(stagger_grid) call phys_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,snew)

    ! check and optionally correct unphysical values
    if(fix_small_values) then
       call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImin3,&
          ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,'multi-D finite_volume')
    end if
 
    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim),&
       dtfactor*dble(idimsmax-idimsmin+1)/dble(ndim),ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       1,nw,qtC,wCT,wprim,qt,wnew,x,.false.,active)

  end associate
  contains

    subroutine get_Riemann_flux_tvdmu()
      do iw=iwstart,nwflux
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=half*(fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3, iw))
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,&
            idims)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3, iw)
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf(iws,iwe)
      integer, intent(in) :: iws,iwe
      double precision :: fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)

      fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3) = -0.5d0*tvdlfeps*cmaxC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,ii)
      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=iws,iwe
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=0.5d0*(fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3, iw))
         ! Add TVDLF dissipation to the flux
         if (flux_type(idims, iw) /= flux_no_dissipation) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                iw)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                iw) + fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,iw))
         end if
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,&
            idims)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3, iw)
      end do ! Next iw

    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll(iws,iwe)
      integer, intent(in) :: iws,iwe
      integer :: ix1,ix2,ix3

      do iw=iws,iwe
        if(flux_type(idims, iw) == flux_tvdlf) then
          ! CT MHD does not need normal B flux
          if(stagger_grid) cycle
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,&
             idims) = -tvdlfeps*half*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ii),dabs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ii))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw))
        else
         do ix3=ixCmin3,ixCmax3
         do ix2=ixCmin2,ixCmax2
         do ix1=ixCmin1,ixCmax1
           if(cminC(ix1,ix2,ix3,ii) >= zero) then
             fC(ix1,ix2,ix3,iw,idims)=fLC(ix1,ix2,ix3,iw)
           else if(cmaxC(ix1,ix2,ix3,ii) <= zero) then
             fC(ix1,ix2,ix3,iw,idims)=fRC(ix1,ix2,ix3,iw)
           else
             ! Add hll dissipation to the flux
             fC(ix1,ix2,ix3,iw,idims)=(cmaxC(ix1,ix2,ix3,ii)*fLC(ix1,ix2,ix3,&
                 iw)-cminC(ix1,ix2,ix3,ii)*fRC(ix1,ix2,ix3,iw)+cminC(ix1,ix2,&
                ix3,ii)*cmaxC(ix1,ix2,ix3,ii)*(wRC(ix1,ix2,ix3,iw)-wLC(ix1,ix2,&
                ix3,iw)))/(cmaxC(ix1,ix2,ix3,ii)-cminC(ix1,ix2,ix3,ii))
           end if
         end do
         end do
         end do
       endif 
      end do

    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc(iws,iwe)
      integer, intent(in) :: iws, iwe  
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:nwflux)     :: whll, Fhll, fCD
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)              :: lambdaCD

      integer  :: rho_, p_, e_, mom(1:ndir)

      rho_ = iw_rho
      if (allocated(iw_mom)) mom(:) = iw_mom(:)
      e_ = iw_e 

      if(associated(phys_hllc_init_species)) then
       call phys_hllc_init_species(ii, rho_, mom(:), e_)
      endif  

      p_ = e_

      patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) =  1
      where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1) >= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) = -2
      elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1) <= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method==fs_hllcd) call phys_diffuse_hllcd(ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
         ixCmax3,idims,wLC,wRC,fLC,fRC,patchf)

      !---- calculate speed lambda at CD ----!
      if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)==1)) call phys_get_lCD(wLC,wRC,fLC,fRC,&
         cminC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii),&
         cmaxC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii),idims,&
         ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,&
         ixCmin3,ixCmax1,ixCmax2,ixCmax3, whll,Fhll,lambdaCD,patchf)

      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))== 1))then
         !======== flux at intermediate state ========!
         call phys_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
            cminC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii),&
            cmaxC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii),ixImin1,&
            ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,ixCmin2,ixCmin3,&
            ixCmax1,ixCmax2,ixCmax3,idims,fCD)
      endif ! Calculate the CD flux

      do iw=iws,iwe
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               iw)=-tvdlfeps*half*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,ii),abs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,ii))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,iw) - wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,iw))
         else
            where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)==-2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  iw)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
            elsewhere(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3))==1)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  iw)=fCD(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)==2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  iw)=fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)==3)
               ! fallback option, reducing to HLL flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  iw)=Fhll(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)==4)
               ! fallback option, reducing to TVDLF flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  iw) = half*((fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  ixCmin3:ixCmax3,iw)) -tvdlfeps * max(cmaxC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,ixCmin3:ixCmax3,ii),&
                   dabs(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  ii))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                  iw)))
            endwhere
         end if

         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw,&
            idims)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315 and Guo 2016 JCP, 327, 543
    subroutine get_Riemann_flux_hlld(iws,iwe)
      integer, intent(in) :: iws, iwe
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:nwflux) :: w1R,w1L,f1R,f1L,f2R,f2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:nwflux) :: w2R,w2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) :: sm,s1R,s1L,suR,suL,Bx
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) :: pts,ptR,ptL,signBx,r1L,r1R,tmp
      ! velocity from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,ndir) :: vRC, vLC
      ! magnetic field from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,ndir) :: BR, BL
      integer :: ip1,ip2,ip3,idir,ix1,ix2,ix3
      integer  :: rho_, p_, e_, mom(1:ndir), mag(1:ndir)

      associate (sR=>cmaxC,sL=>cminC)

      rho_ = iw_rho
      mom(:) = iw_mom(:)
      mag(:) = iw_mag(:) 
      e_ = iw_e 
      p_ = e_

      f1R=0.d0
      f1L=0.d0
      f2R=0.d0
      f2L=0.d0
      w1L=0.d0
      w1R=0.d0
      w2L=0.d0
      w2R=0.d0
      ip1=idims
      ip3=3
      vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         :)=wRp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))
      vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         :)=wLp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))
      if(B0field) then
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           :)=wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,:,&
           ip1)
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           :)=wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,:,&
           ip1)
      else
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           :)=wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:))
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           :)=wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:))
      end if
      if(stagger_grid) then
        Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=block%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1)
      else
        ! HLL estimation of normal magnetic field at cell interfaces
        ! Li, Shenghai, 2005 JCP, 203, 344, equation (33)
        Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ii)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ii)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ii)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ii))
      end if
      ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=wRp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,p_)+0.5d0*sum(BR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,:)**2,dim=ndim+1)
      ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=wLp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,p_)+0.5d0*sum(BL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,:)**2,dim=ndim+1)
      if(iw_equi_rho>0) then
        suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3) = wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,rho_)+ block%equi_vars(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw_equi_rho,ip1)
      else
        suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3) = wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,rho_)
      endif
      suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=(sR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,ii)-vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1))*suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      if(iw_equi_rho>0) then
        suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3) = wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,rho_)+ block%equi_vars(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw_equi_rho,ip1)
      else
        suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3) = wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,rho_)
      endif
      suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=(sL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,ii)-vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1))*suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      ! Miyoshi equation (38) and Guo euqation (20)
      sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,ip1)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1)-ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))/(suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))
      ! Miyoshi equation (39) and Guo euqation (28)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip1))=sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip1))=sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip1))=sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip1))=sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      ! Guo equation (22)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))=Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))=Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      if(B0field) then
        ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=wRp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,p_)+0.5d0*sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:))**2,dim=ndim+1)
        ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=wLp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,p_)+0.5d0*sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:))**2,dim=ndim+1)
      end if

      ! Miyoshi equation (43) and Guo equation (27)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         rho_)=suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         rho_)=suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))

      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)*(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)**2
      where(abs(r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))>smalldouble)
        r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=1.d0/r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
      else where
        r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=0.d0
      end where
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)*(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)**2
      where(abs(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))>smalldouble)
        r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=1.d0/r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
      else where
        r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=0.d0
      end where
      ! Miyoshi equation (44)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))=vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))=vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      ! partial solution for later usage
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=(suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ii)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)**2)*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=(suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ii)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)**2)*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        ! Miyoshi equation (46)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))=vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))=vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        ! Miyoshi equation (47)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))=BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))=BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip2))
      end if
      ! Miyoshi equation (45)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip2))
      if(B0field) then
        ! Guo equation (26)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,:,&
           ip1)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,:,&
           ip1)
      end if
      ! equation (48)
      if(phys_energy) then
        ! Guo equation (25) equivalent to Miyoshi equation (41)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           p_)=suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1))+ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        !w1L(ixC^S,p_)=suL(ixC^S)*(sm(ixC^S)-vLC(ixC^S,ip1))+ptL(ixC^S)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           p_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,p_)
        if(B0field) then
          ! Guo equation (32)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             p_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             :,ip1)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))),dim=ndim+1)
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             p_)=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             :,ip1)*(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))),dim=ndim+1)
        end if
        ! Miyoshi equation (48) and main part of Guo euqation (31)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)=((sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ii)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip1))*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)-ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1)+w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,p_)*sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)+Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(sum(vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,:)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:)),dim=ndim+1)-sum(w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))*w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:)),&
           dim=ndim+1)))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)=((sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ii)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip1))*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)-ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1)+w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,p_)*sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)+Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(sum(vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,:)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:)),dim=ndim+1)-sum(w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))*w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:)),&
           dim=ndim+1)))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
        if(B0field) then
          ! Guo equation (31)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)+(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             :,ip1),dim=ndim+1)*sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)-sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,:,ip1),dim=ndim+1)*vRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1))/(sR(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,ii)-sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3))
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)+(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             :,ip1),dim=ndim+1)*sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)-sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,:,ip1),dim=ndim+1)*vLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1))/(sL(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,ii)-sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3))
        end if
        if(iw_equi_p>0) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)= w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_) + 1d0/(phys_gamma - 1) * block%equi_vars(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw_equi_p,&
             ip1) * (sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ip1))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3))
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)= w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_) + 1d0/(phys_gamma - 1) * block%equi_vars(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw_equi_p,&
             ip1) * (sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ip1))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3))
#else
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)= w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_) + phys_gamma /(phys_gamma - 1) * &
             block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw_equi_p,ip1) * (sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ip1))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3))
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)= w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_) + phys_gamma /(phys_gamma - 1) * &
             block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw_equi_p,ip1) * (sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ip1))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,ii)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3))
#endif
        endif
      end if

      ! Miyoshi equation (49) and Guo equation (35)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         rho_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,rho_)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         rho_)=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,rho_)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))

      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=sqrt(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,rho_))
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=sqrt(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,rho_))
      tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=1.d0/(r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))
      signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sign(1.d0,&
         Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
      ! Miyoshi equation (51) and Guo equation (33)
      s1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)+abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3))/r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      s1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)-abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3))/r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      ! Miyoshi equation (59) and Guo equation (41)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))=(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2)))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))
      ! Miyoshi equation (61) and Guo equation (43)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2)))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))
      if(ndir==3) then
        ! Miyoshi equation (60) and Guo equation (42)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))=(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(ip3))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(ip3)))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))
        ! Miyoshi equation (62) and Guo equation (44)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))=(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(ip3))+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(ip3)))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))
      end if
      ! Miyoshi equation (63) and Guo equation (45)
      if(phys_energy) then
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:)),dim=ndim+1)-sum(w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))*w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)-r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:)),dim=ndim+1)-sum(w2L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))*w2L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
      end if

      ! convert velocity to momentum
      do idir=1,ndir
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))=w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_)
      end do
      if(iw_equi_rho>0) then
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_) = w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_) - block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,iw_equi_rho,ip1)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_) = w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_) - block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,iw_equi_rho,ip1)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_) = w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_) - block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,iw_equi_rho,ip1)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_) = w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_) - block%equi_vars(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,iw_equi_rho,ip1)
      endif
      ! get fluxes of intermedate states
      do iw=iws,iwe
        if(flux_type(idims, iw) == flux_special) then
          ! known flux (fLC=fRC) for normal B and psi_ in GLM method
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
        else if(flux_type(idims, iw) == flux_hll) then
          ! using hll flux for tracers
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             ii)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              iw)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             ii)*fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              iw) +sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             ii)*sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             ii)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             ii)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,ii))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
        else
          ! construct hlld flux
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             ii)*(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             ii)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)-wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw))
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+s1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)*(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw)-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw))
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+s1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)*(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw)-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw))
        end if
      end do

      ! Miyoshi equation (66) and Guo equation (46)
     do ix3=ixCmin3,ixCmax3
     do ix2=ixCmin2,ixCmax2
     do ix1=ixCmin1,ixCmax1
        if(sL(ix1,ix2,ix3,ii)>0.d0) then
          fC(ix1,ix2,ix3,iws:iwe,ip1)=fLC(ix1,ix2,ix3,iws:iwe)
        else if(s1L(ix1,ix2,ix3)>=0.d0) then
          fC(ix1,ix2,ix3,iws:iwe,ip1)=f1L(ix1,ix2,ix3,iws:iwe)
        else if(sm(ix1,ix2,ix3)>=0.d0) then
          fC(ix1,ix2,ix3,iws:iwe,ip1)=f2L(ix1,ix2,ix3,iws:iwe)
        else if(s1R(ix1,ix2,ix3)>=0.d0) then
          fC(ix1,ix2,ix3,iws:iwe,ip1)=f2R(ix1,ix2,ix3,iws:iwe)
        else if(sR(ix1,ix2,ix3,ii)>=0.d0) then
          fC(ix1,ix2,ix3,iws:iwe,ip1)=f1R(ix1,ix2,ix3,iws:iwe)
        else if(sR(ix1,ix2,ix3,ii)<0.d0) then
          fC(ix1,ix2,ix3,iws:iwe,ip1)=fRC(ix1,ix2,ix3,iws:iwe)
        end if
     end do
     end do
     end do

      end associate
    end subroutine get_Riemann_flux_hlld

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315 and Guo 2016 JCP, 327, 543
    !> https://arxiv.org/pdf/2108.04991.pdf
    subroutine get_Riemann_flux_hlld_mag2()
      !use mod_mhd_phys
      use mod_variables
      use mod_physics
      implicit none
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:nwflux) :: w1R,w1L,f1R,f1L,f2R,f2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:nwflux) :: w2R,w2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) :: sm,s1R,s1L,suR,suL,Bx
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) :: pts,ptR,ptL,signBx,r1L,r1R,tmp
      ! velocity from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,ndir) :: vRC, vLC
      ! magnetic field from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,ndir) :: BR, BL
      integer :: ip1,ip2,ip3,idir,ix1,ix2,ix3
      double precision :: phiPres, thetaSM, du, dv, dw
      integer :: ixVmin1,ixVmin2,ixVmin3,ixVmax1,ixVmax2,ixVmax3, ixVbmin1,&
         ixVbmin2,ixVbmin3,ixVbmax1,ixVbmax2,ixVbmax3, ixVcmin1,ixVcmin2,&
         ixVcmin3,ixVcmax1,ixVcmax2,ixVcmax3, ixVdmin1,ixVdmin2,ixVdmin3,&
         ixVdmax1,ixVdmax2,ixVdmax3, ixVemin1,ixVemin2,ixVemin3,ixVemax1,&
         ixVemax2,ixVemax3, ixVfmin1,ixVfmin2,ixVfmin3,ixVfmax1,ixVfmax2,&
         ixVfmax3
      integer  :: rho_, p_, e_, mom(1:ndir), mag(1:ndir)
      double precision, parameter :: aParam = 4d0

      rho_ = iw_rho
      mom(:) = iw_mom(:)
      mag(:) = iw_mag(:) 
      p_ = iw_e
      e_ = iw_e 

      associate (sR=>cmaxC,sL=>cminC)

      f1R=0.d0
      f1L=0.d0
      ip1=idims
      ip3=3

      vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         :)=wRp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))
      vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         :)=wLp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))

      ! reuse s1L s1R
      call get_hlld2_modif_c(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,s1L)
      call get_hlld2_modif_c(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,s1R)
      !phiPres = min(1, maxval(max(s1L(ixO^S),s1R(ixO^S))/cmaxC(ixO^S,1)))
      phiPres = min(1d0, maxval(max(s1L(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),s1R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))/maxval(cmaxC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)))
      phiPres = phiPres*(2D0 - phiPres)  

     !we use here not reconstructed velocity: wprim?       
     ixVmin1=ixOmin1;ixVmin2=ixOmin2;ixVmin3=ixOmin3;ixVmax1=ixOmax1
     ixVmax2=ixOmax2;ixVmax3=ixOmax3;
     !first dim
     ixVmin1=ixOmin1+1  
     ixVmax1=ixOmax1+1
     du = minval(wprim(ixVmin1:ixVmax1,ixVmin2:ixVmax2,ixVmin3:ixVmax3,&
        mom(1))-wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1)))
     if(du>0d0) du=0d0  
     dv = 0d0 
     dw = 0d0 

     
     !second dim
     !i,j-1,k 
     ixVmin1=ixOmin1;ixVmin2=ixOmin2;ixVmin3=ixOmin3;ixVmax1=ixOmax1
     ixVmax2=ixOmax2;ixVmax3=ixOmax3;
     ixVmin2=ixOmin2-1  
     ixVmax2=ixOmax2-1
 
     !i,j+1,k 
     ixVbmin1=ixOmin1;ixVbmin2=ixOmin2;ixVbmin3=ixOmin3;ixVbmax1=ixOmax1
     ixVbmax2=ixOmax2;ixVbmax3=ixOmax3;
     ixVbmin2=ixOmin2+1  
     ixVbmax2=ixOmax2+1

     !i+1,j,k 
     ixVcmin1=ixOmin1;ixVcmin2=ixOmin2;ixVcmin3=ixOmin3;ixVcmax1=ixOmax1
     ixVcmax2=ixOmax2;ixVcmax3=ixOmax3;
     ixVcmin1=ixOmin1+1  
     ixVcmax1=ixOmax1+1

     !i+1,j-1,k 
     ixVdmin1=ixOmin1;ixVdmin2=ixOmin2;ixVdmin3=ixOmin3;ixVdmax1=ixOmax1
     ixVdmax2=ixOmax2;ixVdmax3=ixOmax3;
     ixVdmin1=ixOmin1+1  
     ixVdmax1=ixOmax1+1
     ixVdmin2=ixOmin2-1  
     ixVdmax2=ixOmax2-1

     !i+1,j+1,k 
     ixVemin1=ixOmin1;ixVemin2=ixOmin2;ixVemin3=ixOmin3;ixVemax1=ixOmax1
     ixVemax2=ixOmax2;ixVemax3=ixOmax3;
     ixVemin1=ixOmin1+1  
     ixVemax1=ixOmax1+1
     ixVemin2=ixOmin2+1  
     ixVemax2=ixOmax2+1

     dv = minval(min(wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        mom(2))-wprim(ixVmin1:ixVmax1,ixVmin2:ixVmax2,ixVmin3:ixVmax3,mom(2)),&
        wprim(ixVbmin1:ixVbmax1,ixVbmin2:ixVbmax2,ixVbmin3:ixVbmax3,&
        mom(2))-wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2)),&
        wprim(ixVcmin1:ixVcmax1,ixVcmin2:ixVcmax2,ixVcmin3:ixVcmax3,&
        mom(2))-wprim(ixVdmin1:ixVdmax1,ixVdmin2:ixVdmax2,ixVdmin3:ixVdmax3,&
        mom(2)),wprim(ixVemin1:ixVemax1,ixVemin2:ixVemax2,ixVemin3:ixVemax3,&
        mom(2))-wprim(ixVcmin1:ixVcmax1,ixVcmin2:ixVcmax2,ixVcmin3:ixVcmax3,&
        mom(2))))
     if(dv>0d0) dv=0d0 

      
     !third dim
     !i,j,k-1 
     ixVmin1=ixOmin1;ixVmin2=ixOmin2;ixVmin3=ixOmin3;ixVmax1=ixOmax1
     ixVmax2=ixOmax2;ixVmax3=ixOmax3;
     ixVmin3=ixOmin3-1  
     ixVmax3=ixOmax3-1
 
     !i,j,k+1 
     ixVbmin1=ixOmin1;ixVbmin2=ixOmin2;ixVbmin3=ixOmin3;ixVbmax1=ixOmax1
     ixVbmax2=ixOmax2;ixVbmax3=ixOmax3;
     ixVbmin3=ixOmin3+1  
     ixVbmax3=ixOmax3+1

     !i+1,j,k 
     ixVcmin1=ixOmin1;ixVcmin2=ixOmin2;ixVcmin3=ixOmin3;ixVcmax1=ixOmax1
     ixVcmax2=ixOmax2;ixVcmax3=ixOmax3;
     ixVcmin1=ixOmin1+1  
     ixVcmax1=ixOmax1+1

     !i+1,j,k-1 
     ixVdmin1=ixOmin1;ixVdmin2=ixOmin2;ixVdmin3=ixOmin3;ixVdmax1=ixOmax1
     ixVdmax2=ixOmax2;ixVdmax3=ixOmax3;
     ixVdmin1=ixOmin1+1  
     ixVdmax1=ixOmax1+1
     ixVdmin3=ixOmin3-1  
     ixVdmax3=ixOmax3-1

     !i+1,j,k+1 
     ixVemin1=ixOmin1;ixVemin2=ixOmin2;ixVemin3=ixOmin3;ixVemax1=ixOmax1
     ixVemax2=ixOmax2;ixVemax3=ixOmax3;
     ixVemin1=ixOmin1+1  
     ixVemax1=ixOmax1+1
     ixVemin3=ixOmin3+1  
     ixVemax3=ixOmax3+1
     dw = minval(min(wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        mom(3))-wprim(ixVmin1:ixVmax1,ixVmin2:ixVmax2,ixVmin3:ixVmax3,mom(3)),&
        wprim(ixVbmin1:ixVbmax1,ixVbmin2:ixVbmax2,ixVbmin3:ixVbmax3,&
        mom(3))-wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3)),&
        wprim(ixVcmin1:ixVcmax1,ixVcmin2:ixVcmax2,ixVcmin3:ixVcmax3,&
        mom(3))-wprim(ixVdmin1:ixVdmax1,ixVdmin2:ixVdmax2,ixVdmin3:ixVdmax3,&
        mom(3)),wprim(ixVemin1:ixVemax1,ixVemin2:ixVemax2,ixVemin3:ixVemax3,&
        mom(3))-wprim(ixVcmin1:ixVcmax1,ixVcmin2:ixVcmax2,ixVcmin3:ixVcmax3,&
        mom(3))))
     if(dw>0d0) dw=0d0
     thetaSM = maxval(cmaxC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        1)) 
       
     thetaSM = (min(1d0, (thetaSM-du)/(thetaSM-min(dv,dw))))**aParam
     !print*, "HLLD2 ", du,dv,dw, thetaSM, phiPres 

      if(B0field) then
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           :)=wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,:,&
           ip1)
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           :)=wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,:,&
           ip1)
      else
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           :)=wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:))
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           :)=wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:))
      end if
      ! HLL estimation of normal magnetic field at cell interfaces
      Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=(sR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,index_v_mag)*BR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1)-sL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,index_v_mag)*BL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1)-fLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip1))-fRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip1)))/(sR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,index_v_mag)-sL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,index_v_mag))
      ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=wRp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,p_)+0.5d0*sum(BR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,:)**2,dim=ndim+1)
      ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=wLp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,p_)+0.5d0*sum(BL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,:)**2,dim=ndim+1)
      suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=(sR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,index_v_mag)-vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1))*wRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,rho_)
      suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=(sL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,index_v_mag)-vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1))*wLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3,rho_)
      ! Miyoshi equation (38) and Guo euqation (20)
      sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,ip1)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1)-thetaSM*(ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)-ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)) )/(suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))
      ! Miyoshi equation (39) and Guo euqation (28)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip1))=sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip1))=sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip1))=sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip1))=sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      ! Guo equation (22)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))=Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))=Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      if(B0field) then
        ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=wRp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,p_)+0.5d0*sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:))**2,dim=ndim+1)
        ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=wLp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,p_)+0.5d0*sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:))**2,dim=ndim+1)
      end if

      ! Miyoshi equation (43) and Guo equation (27)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         rho_)=suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         rho_)=suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))

      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)*(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)**2
      where(r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)/=0.d0)
        r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=1.d0/r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
      endwhere
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)*(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)**2
      where(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)/=0.d0)
        r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)=1.d0/r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
      endwhere
      ! Miyoshi equation (44)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))=vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))=vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
      ! partial solution for later usage
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=(suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         index_v_mag)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)**2)*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=(suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         index_v_mag)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)**2)*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        ! Miyoshi equation (46)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))=vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))=vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        ! Miyoshi equation (47)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))=BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))=BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip2))
      end if
      ! Miyoshi equation (45)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         ip2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(ip2))
      if(B0field) then
        ! Guo equation (26)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,:,&
           ip1)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,:,&
           ip1)
      end if
      ! equation (48)
      if(phys_energy) then
        ! Guo equation (25) equivalent to Miyoshi equation (41)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           p_)=(suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3) - suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3) +phiPres * suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1)))/(suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           p_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,p_)
        if(B0field) then
          ! Guo equation (32)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             p_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             :,ip1)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))),dim=ndim+1)
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             p_)=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             :,ip1)*(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))),dim=ndim+1)
        end if
        ! Miyoshi equation (48) and main part of Guo euqation (31)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)=((sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           index_v_mag)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip1))*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)-ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1)+w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,p_)*sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)+Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(sum(vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,:)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:)),dim=ndim+1)-sum(w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))*w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:)),&
           dim=ndim+1)))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)=((sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           index_v_mag)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           ip1))*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)-ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,ip1)+w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,p_)*sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)+Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(sum(vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,:)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:)),dim=ndim+1)-sum(w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))*w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:)),&
           dim=ndim+1)))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           index_v_mag)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
        if(B0field) then
          ! Guo equation (31)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)+(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             :,ip1),dim=ndim+1)*sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)-sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,:,ip1),dim=ndim+1)*vRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1))/(sR(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,index_v_mag)-sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3))
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             e_)+(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             :,ip1),dim=ndim+1)*sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)-sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,:,ip1),dim=ndim+1)*vLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,ip1))/(sL(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,index_v_mag)-sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3))
        end if
      end if

      ! Miyoshi equation (49) and Guo equation (35)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         rho_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,rho_)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         rho_)=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,rho_)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip1))

      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=sqrt(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,rho_))
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=sqrt(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,rho_))
      tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)=1.d0/(r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))
      signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sign(1.d0,&
         Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3))
      ! Miyoshi equation (51) and Guo equation (33)
      s1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)+abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3))/r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      s1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3)-abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ixCmin3:ixCmax3))/r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      ! Miyoshi equation (59) and Guo equation (41)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))=(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2)))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))
      ! Miyoshi equation (61) and Guo equation (43)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(ip2)))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mag(ip2))
      if(ndir==3) then
        ! Miyoshi equation (60) and Guo equation (42)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))=(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(ip3))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(ip3)))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(ip3))
        ! Miyoshi equation (62) and Guo equation (44)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))=(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(ip3))+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(ip3)))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mag(ip3))
      end if
      ! Miyoshi equation (63) and Guo equation (45)
      if(phys_energy) then
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:)),dim=ndim+1)-sum(w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))*w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           e_)-r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)*(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3,mag(:)),dim=ndim+1)-sum(w2L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mom(:))*w2L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ixCmin3:ixCmax3)
      end if

      ! convert velocity to momentum
      do idir=1,ndir
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))=w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))=w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))=w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))=w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           mom(idir))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           rho_)
      end do

      ! get fluxes of intermedate states
      do iw=1,nwflux
        ! CT MHD does not need normal B flux
        if(stagger_grid .and. flux_type(idims, iw) == flux_tvdlf) cycle
        if(flux_type(idims, iw) == flux_special) then
          ! known flux (fLC=fRC) for normal B and psi_ in GLM method
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
        else if(flux_type(idims, iw) == flux_hll) then
          ! using hll flux for tracers
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             index_v_mag)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              iw)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             index_v_mag)*fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              iw) +sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             index_v_mag)*sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             index_v_mag)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             index_v_mag)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             index_v_mag))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw)
        else
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             index_v_mag)*(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             index_v_mag)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)-wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iw))
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+s1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)*(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw)-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw))
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)=f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             iw)+s1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)*(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw)-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw))
        end if
      end do

      ! Miyoshi equation (66) and Guo equation (46)
     do ix3=ixCmin3,ixCmax3
     do ix2=ixCmin2,ixCmax2
     do ix1=ixCmin1,ixCmax1
        if(sL(ix1,ix2,ix3,index_v_mag)>0.d0) then
          fC(ix1,ix2,ix3,1:nwflux,ip1)=fLC(ix1,ix2,ix3,1:nwflux)
        else if(s1L(ix1,ix2,ix3)>=0.d0) then
          fC(ix1,ix2,ix3,1:nwflux,ip1)=f1L(ix1,ix2,ix3,1:nwflux)
        else if(sm(ix1,ix2,ix3)>=0.d0) then
          fC(ix1,ix2,ix3,1:nwflux,ip1)=f2L(ix1,ix2,ix3,1:nwflux)
        else if(s1R(ix1,ix2,ix3)>=0.d0) then
          fC(ix1,ix2,ix3,1:nwflux,ip1)=f2R(ix1,ix2,ix3,1:nwflux)
        else if(sR(ix1,ix2,ix3,index_v_mag)>=0.d0) then
          fC(ix1,ix2,ix3,1:nwflux,ip1)=f1R(ix1,ix2,ix3,1:nwflux)
        else if(sR(ix1,ix2,ix3,index_v_mag)<0.d0) then
          fC(ix1,ix2,ix3,1:nwflux,ip1)=fRC(ix1,ix2,ix3,1:nwflux)
        end if
     end do
     end do
     end do

      end associate
    end subroutine get_Riemann_flux_hlld_mag2

    !> Calculate fast magnetosonic wave speed
    subroutine get_hlld2_modif_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound)
      use mod_global_parameters

      integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
      double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3), AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3), b2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
          kmax
      double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3), gamma_A2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      integer  :: rho_, p_, e_, mom(1:ndir), mag(1:ndir)

        rho_ = iw_rho
        mom(:) = iw_mom(:)
        mag(:) = iw_mag(:) 
        p_ = iw_e
        e_ = iw_e 

      inv_rho=1.d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

      ! store |B|^2 in v

      if (B0field) then
        b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,:,b0i))**2, dim=ndim+1)
      else
        b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mag(:))**2, dim=ndim+1)
      end if


      if (B0field) then
        AvMinCs2= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(idims))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idims,b0i)
      else
        AvMinCs2= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(idims))
      end if


      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mom(:))**2, dim=ndim+1)

      cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)   = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2-4.0d0*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2 * inv_rho 

      where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
         AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
      end where

      AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))

     csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        ixOmin3:ixOmax3) = sqrt(half*(cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        ixOmin3:ixOmax3)))

    end subroutine get_hlld2_modif_c

  end subroutine finite_volume

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,&
     ixRmax1,ixRmax2,ixRmax3,idims,w,wLC,wRC,wLp,wRp,x,dxdim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3, ixRmin1,ixRmin2,&
       ixRmin3,ixRmax1,ixRmax2,ixRmax3, idims
    double precision, intent(in) :: dxdim
    ! cell center w in primitive form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: w
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wLp, wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) :: x

    integer            :: jxRmin1,jxRmin2,jxRmin3,jxRmax1,jxRmax2,jxRmax3,&
        ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, jxCmin1,jxCmin2,&
       jxCmin3,jxCmax1,jxCmax2,jxCmax3, iw
    double precision   :: ldw(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        rdw(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        dwC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision   :: a2max

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixLmin1,&
          ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,w,wLp,wRp)
    case (limiter_weno3)
       call WENO3limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp,1)
    case (limiter_wenoyc3)
       call WENO3limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp,2)
    case (limiter_weno5)
       call WENO5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp,1)
    case (limiter_weno5nm)
       call WENO5NMlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp,1)
    case (limiter_wenoz5)
       call WENO5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp,2)
    case (limiter_wenozp5)
       call WENO5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp,3)
    case (limiter_weno5cu6)
       call WENO5CU6limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,w,wLp,wRp)
    case (limiter_teno5ad)
       call TENO5ADlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp)
    case (limiter_weno7)
       call WENO7limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,w,wLp,wRp,1)
    case (limiter_mpweno7)
       call WENO7limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,w,wLp,wRp,2)
    case (limiter_venk)
       call venklimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idims,dxdim,w,wLp,&
          wRp) 
       if(fix_small_values) then
          call phys_handle_small_values(.true.,wLp,x,ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,&
             ixLmax3,'reconstruct left')
          call phys_handle_small_values(.true.,wRp,x,ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
             ixRmax3,'reconstruct right')
       end if
    case (limiter_ppm)
       call PPMlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixMlo1,&
          ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,idims,w,w,wLp,wRp)
       if(fix_small_values) then
          call phys_handle_small_values(.true.,wLp,x,ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,&
             ixLmax3,'reconstruct left')
          call phys_handle_small_values(.true.,wRp,x,ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
             ixRmax3,'reconstruct right')
       end if
    case default
       jxRmin1=ixRmin1+kr(idims,1);jxRmin2=ixRmin2+kr(idims,2)
       jxRmin3=ixRmin3+kr(idims,3);jxRmax1=ixRmax1+kr(idims,1)
       jxRmax2=ixRmax2+kr(idims,2);jxRmax3=ixRmax3+kr(idims,3);
       ixCmax1=jxRmax1;ixCmax2=jxRmax2;ixCmax3=jxRmax3
       ixCmin1=ixLmin1-kr(idims,1);ixCmin2=ixLmin2-kr(idims,2)
       ixCmin3=ixLmin3-kr(idims,3);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
       jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);
       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,ixCmin3:jxCmax3,&
                iw)=dlog10(w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,ixCmin3:jxCmax3,&
                iw))
             wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
                iw)=dlog10(wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
                iw))
             wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                iw)=dlog10(wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                iw))
          end if

          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             jxCmin3:jxCmax3,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw)
          if(need_global_a2max) then 
            a2max=a2max_global(idims)
          else
            select case(idims)
            case(1)
              a2max=schmid_rad1
            
            
            case(2)
              a2max=schmid_rad2
            case(3)
              a2max=schmid_rad3
            case default
              call mpistop("idims is wrong in mod_limiter")
            end select
          end if

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,&
             type_limiter(block%level),ldw,rdw,a2max=a2max)
          wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
             iw)=wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
             iw)+half*ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3)
          wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
             iw)=wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
             iw)-half*rdw(jxRmin1:jxRmax1,jxRmin2:jxRmax2,jxRmin3:jxRmax3)

          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,ixCmin3:jxCmax3,&
                iw)=10.0d0**w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,ixCmin3:jxCmax3,&
                iw)
             wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
                iw)=10.0d0**wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                ixLmin3:ixLmax3,iw)
             wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                iw)=10.0d0**wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                ixRmin3:ixRmax3,iw)
          end if
       end do
       if(fix_small_values) then
          call phys_handle_small_values(.true.,wLp,x,ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,&
             ixLmax3,'reconstruct left')
          call phys_handle_small_values(.true.,wRp,x,ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,&
             ixRmax3,'reconstruct right')
       end if
    end select

   wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
      1:nwflux) = wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
      1:nwflux)
   wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
      1:nwflux) = wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
      1:nwflux)
   call phys_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,wLC,x)
   call phys_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixRmin1,ixRmin2,ixRmin3,ixRmax1,ixRmax2,ixRmax3,wRC,x)
   if(nwaux>0)then
      wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
         nwflux+1:nwflux+nwaux) = wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
         ixLmin3:ixLmax3,nwflux+1:nwflux+nwaux)
      wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
         nwflux+1:nwflux+nwaux) = wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
         ixRmin3:ixRmax3,nwflux+1:nwflux+nwaux)
   endif

  end subroutine reconstruct_LR

end module mod_finite_volume
