!> Nicolas Moens
!> Module for including flux limited diffusion (FLD)-approximation in Radiation-hydrodynamics simulations using mod_rhd
!> Based on Turner and stone 2001. See
!> [1]Moens, N., Sundqvist, J. O., El Mellah, I., Poniatowski, L., Teunissen, J., and Keppens, R.,
!> “Radiation-hydrodynamics with MPI-AMRVAC . Flux-limited diffusion”,
!> <i>Astronomy and Astrophysics</i>, vol. 657, 2022. doi:10.1051/0004-6361/202141023.
!> For more information.

module mod_fld

    use mod_comm_lib, only: mpistop

    implicit none

    !> source split for energy interact and radforce:
    logical :: fld_Eint_split = .false.
    logical :: fld_Radforce_split = .false.

    !> Opacity per unit of unit_density
    double precision, public :: fld_kappa0 = 0.34d0

    !> mean particle mass
    double precision, public :: fld_mu = 0.6d0

    !> Tolerance for bisection method for Energy sourceterms
    !> This is a percentage of the minimum of gas- and radiation energy
    double precision, public :: fld_bisect_tol = 1.d-4

    !> Tolerance for adi method for radiative Energy diffusion
    double precision, public :: fld_diff_tol = 1.d-4

    !> Number for splitting the diffusion module
    double precision, public :: diff_crit

    !> Use constant Opacity?
    character(len=8)  :: fld_opacity_law = 'const'
    character(len=50) :: fld_opal_table = 'Y09800' !>'xxxxxx'

    !> Diffusion limit lambda = 0.33
    character(len=16) :: fld_fluxlimiter = 'Pomraning'

    !> diffusion coefficient for multigrid method
    integer :: i_diff_mg

    !> Which method to solve diffusion part
    character(len=8) :: fld_diff_scheme = 'mg'

    !> Which method to find the root for the energy interaction polynomial
    character(len=8) :: fld_interaction_method = 'Halley'

    !> Take running average for Diffusion coefficient
    logical :: diff_coef_filter = .false.
    integer :: size_D_filter = 1

    !> Take a running average over the fluxlimiter
    logical :: flux_lim_filter = .false.
    integer :: size_L_filter = 1

    !> Use or dont use lineforce opacities
    logical :: Lineforce_opacities = .false.

    !> Resume run when multigrid returns error
    logical :: diffcrash_resume = .true.

    !> Index for Flux weighted opacities
    integer, allocatable, public :: i_opf(:)

    !> A copy of rhd_Gamma
    double precision, private, protected :: fld_gamma

    !> running timestep for diffusion solver, initialised as zero
    double precision :: dt_diff = 0.d0

    !> public methods
    !> these are called in mod_rhd_phys
    public :: get_fld_rad_force
    public :: get_fld_energy_interact
    public :: fld_radforce_get_dt
    public :: fld_init
    public :: fld_get_radflux
    public :: fld_get_radpress
    public :: fld_get_fluxlimiter
    public :: fld_get_opacity

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! GENERAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reading in fld-list parameters from .par file
  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa0, fld_Eint_split, fld_Radforce_split,&
        fld_bisect_tol, fld_diff_tol,fld_opacity_law, fld_fluxlimiter,&
        fld_diff_scheme, fld_interaction_method, diff_coef_filter,&
        size_D_filter, flux_lim_filter, size_L_filter, lineforce_opacities,&
        diffcrash_resume, fld_opal_table

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, fld_list, end=111)
       111    close(unitpar)
    end do

  end subroutine fld_params_read

  !> Initialising FLD-module:
  !> Read opacities
  !> Initialise Multigrid
  !> adimensionalise kappa
  !> Add extra variables to w-array, flux, kappa, eddington Tensor
  !> Lambda and R
  !> ...
  subroutine fld_init(He_abundance, rhd_radiation_diffusion, rhd_gamma)
    use mod_global_parameters
    use mod_variables
    use mod_physics
    use mod_opal_opacity, only: init_opal_table
    use mod_multigrid_coupling

    double precision, intent(in) :: He_abundance, rhd_gamma
    logical, intent(in) :: rhd_radiation_diffusion
    double precision :: sigma_thomson
    integer :: idir,jdir

    character(len=1) :: ind_1
    character(len=1) :: ind_2
    character(len=2) :: cmp_f
    character(len=5) :: cmp_e

    !> read par files
    call fld_params_read(par_files)

    !> Set lineforce opacities as variable
    if (lineforce_opacities) then
      allocate(i_opf(ndim))
      do idir = 1,ndim
        write(ind_1,'(I1)') idir
        cmp_f = 'k' // ind_1
        i_opf(idir) = var_set_extravar(cmp_f,cmp_f)
      enddo
    endif

    if (rhd_radiation_diffusion) then
      if (fld_diff_scheme .eq. 'mg') then

        use_multigrid = .true.

        phys_implicit_update => Diffuse_E_rad_mg
        phys_evaluate_implicit => Evaluate_E_rad_mg

        mg%n_extra_vars = 1
        mg%operator_type = mg_vhelmholtz

      endif
    endif

    i_diff_mg = var_set_extravar("D", "D")

    !> Need mean molecular weight
    fld_mu = (1.+4*He_abundance)/(2.+3.*He_abundance)

    !> set rhd_gamma
    fld_gamma = rhd_gamma

    !> Read in opacity table if necesary
    if (fld_opacity_law .eq. 'opal') call init_opal_table(fld_opal_table)
    if ((fld_opacity_law .eq. 'thomson') .or. (fld_opacity_law .eq. &
       'fastwind'))  then
      sigma_thomson = 6.6524585d-25
      fld_kappa0 = sigma_thomson/const_mp * &
         (1.+2.*He_abundance)/(1.+4.*He_abundance)
    endif
  end subroutine fld_init

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  !> This subroutine handles the radiation force
  subroutine get_fld_rad_force(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,energy,&
     qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    double precision :: radiation_forceCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1:ndim)
    double precision :: kappaCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: rad_fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1:ndim)

    double precision :: div_v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim,1:ndim)
    double precision :: edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndim,1:ndim)
    double precision :: nabla_vP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        grad_v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        grad0_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: grad_E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    integer :: idir, jdir

    double precision :: fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_Radforce_split) then
      active = .true.

      call fld_get_opacity(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappaCT)
      call fld_get_radflux(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, rad_fluxCT)
      call fld_get_fluxlimiter(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, lambda,&
          fld_R)

      do idir = 1,ndim
        !> Radiation force = kappa*rho/c *Flux = lambda gradE
        radiation_forceCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir) = kappaCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*rad_fluxCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir)/(const_c/unit_velocity)

        ! call gradientO(wCT(ixI^S,iw_r_e),x,ixI^L,ixO^L,idir,grad_E,nghostcells)
        ! radiation_forceCT(ixO^S,idir) = lambda(ixO^S)*grad_E(ixO^S)

        !> Momentum equation source term
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(idir)) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw_rho)*radiation_forceCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
        ! w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) &
        !     + qdt *radiation_forceCT(ixO^S,idir)

        ! if (energy .and. .not. block%e_is_internal) then
          !> Energy equation source term (kinetic energy)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_e) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_mom(idir))*radiation_forceCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir)
          ! w(ixO^S,iw_e) = w(ixO^S,iw_e) &
          !     + qdt * wCT(ixO^S,iw_mom(idir))/wCT(ixO^S,iw_rho)*radiation_forceCT(ixO^S,idir)
        ! endif
      enddo

      !> Photon tiring
      !> calculate tensor div_v
      !> !$OMP PARALLEL DO
      do idir = 1,ndim
        do jdir = 1,ndim
          vel(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3) = wCt(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,iw_mom(jdir))/wCt(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,iw_rho)

          call gradient(vel,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,grad_v)
          div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,&
             jdir) = grad_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

          ! call gradientO(vel,x,ixI^L,ixO^L,idir,grad0_v,nghostcells)
          ! div_v(ixO^S,idir,jdir) = grad0_v(ixO^S)
        enddo
      enddo
      !> !$OMP END PARALLEL DO

      call fld_get_eddington(wCt, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, edd)

      !> VARIABLE NAMES DIV ARE ACTUALLY GRADIENTS
      

      

      
      nabla_vP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1,1)*edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1,1)  + div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1,2)*edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1,2)  + div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1,3)*edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1,3)  + div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2,1)*edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2,1)  + div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2,2)*edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2,2)  + div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2,3)*edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2,3)  + div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3,1)*edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3,1)  + div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3,2)*edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3,2)  + div_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3,3)*edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3,3)
     

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_r_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_r_e) - qdt * nabla_vP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*wCt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_r_e)

    end if

  end subroutine get_fld_rad_force


  subroutine fld_radforce_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,&
     x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim), w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radiation_force(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1:ndim)
    double precision :: rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1:ndim)
    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision                :: dxinv(1:ndim), max_grad
    integer                         :: idim

    dxinv(1)=one/dx1;dxinv(2)=one/dx2;dxinv(3)=one/dx3;

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)
    call fld_get_radflux(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, rad_flux)

    do idim = 1, ndim
      radiation_force(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_rho)*kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)/(const_c/unit_velocity)
      max_grad = maxval(abs(radiation_force(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)))
      max_grad = max(max_grad, epsilon(1.0d0))
      dtnew = min(dtnew, courantpar / sqrt(max_grad * dxinv(idim)))
    end do

  end subroutine fld_radforce_get_dt

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  !> This subroutine handles the energy exchange between gas and radiation
  subroutine get_fld_energy_interact(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
     energy,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods

    use mod_physics, only: phys_get_pthermal  !needed to get temp

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_Eint_split) then
      active = .true.
      !> Add energy sourceterms
      call Energy_interaction(w, w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    end if
  end subroutine get_fld_energy_interact

  !> Sets the opacity in the w-array
  !> by calling mod_opal_opacity
  subroutine fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, fld_kappa)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal
    use mod_physics, only: phys_get_tgas
    use mod_usr_methods
    use mod_opal_opacity

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out) :: fld_kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        pth(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: rho0,Temp0,n,sigma_b
    double precision :: akram, bkram
    double precision :: vth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        gradv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        t(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    integer :: i,j,ix1,ix2,ix3, idir

    select case (fld_opacity_law)
      case('const')
        fld_kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = fld_kappa0/unit_opacity
      case('thomson')
        fld_kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = fld_kappa0/unit_opacity
      case('kramers')
        rho0 = half !> Take lower value of rho in domain
        fld_kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = fld_kappa0/unit_opacity*((w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_rho)/rho0))
      case('bump')
        !> Opacity bump
        rho0 = 0.2d0 !0.5d-1
        n = 7.d0
        sigma_b = 2.d-2
        !fld_kappa(ixO^S) = fld_kappa0/unit_opacity*(one + n*dexp(-((rho0  - w(ixO^S,iw_rho))**two)/rho0))
        fld_kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = fld_kappa0/unit_opacity*(one + &
           n*dexp(-one/sigma_b*(dlog(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw_rho)/rho0))**two))
      case('non_iso')
        call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Temp)
        Temp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=Temp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_rho)

        rho0 = 0.5d0 !> Take lower value of rho in domain
        Temp0 = one
        n = -7.d0/two
        fld_kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = fld_kappa0/unit_opacity*(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_rho)/rho0)*(Temp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)/Temp0)**n
      case('fastwind')
        call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
        a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_rho)*unit_velocity**2.d0

        akram = 13.1351597305
        bkram = -4.5182188206

        fld_kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = fld_kappa0/unit_opacity * &
           (1.d0+10.d0**akram*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw_rho)*unit_density*(a2(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)/1.d12)**bkram)

        do ix1=ixOmin1,ixOmax1
         do ix2=ixOmin2,ixOmax2
         do ix3=ixOmin3,ixOmax3
        
          !> Hard limit on kappa
          fld_kappa(ix1,ix2,ix3) = min(fld_kappa(ix1,ix2,ix3),&
             2.3d0*fld_kappa0/unit_opacity)

          !> Limit kappa through T
          ! fld_kappa(ix^D) = fld_kappa0/unit_opacity &
          ! * (1.d0+10.d0**akram*w(ix^D,iw_rho)*unit_density &
          ! * (max(a2(ix^D),const_kB*5.9d4/(fld_mu*const_mp))/1.d12)**bkram)
        enddo
         enddo
         enddo
        

      case('opal')
        call phys_get_tgas(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Temp)
        do ix1=ixOmin1,ixOmax1
         do ix2=ixOmin2,ixOmax2
         do ix3=ixOmin3,ixOmax3
        
            rho0 = w(ix1,ix2,ix3,iw_rho)*unit_density
            Temp0 = Temp(ix1,ix2,ix3)*unit_temperature
            call set_opal_opacity(rho0,Temp0,n)
            fld_kappa(ix1,ix2,ix3) = n/unit_opacity
        enddo
         enddo
         enddo
        

      case('special')
        if (.not. associated(usr_special_opacity)) then
          call mpistop("special opacity not defined")
        endif
        call usr_special_opacity(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
            fld_kappa)

      case default
        call mpistop("Doesn't know opacity law")
      end select
  end subroutine fld_get_opacity

  !> Calculate fld flux limiter
  !> This subroutine calculates flux limiter lambda using the prescription
  !> stored in fld_fluxlimiter.
  !> It also calculates the ratio of radiation scaleheight and mean free path
  subroutine fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, fld_lambda,&
      fld_R)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out) :: fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision ::  normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: grad_r_e(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: rad_e(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: idir, ix1,ix2,ix3

    double precision :: tmp_L(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        filtered_L(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: filter, idim

    select case (fld_fluxlimiter)
    case('Diffusion')
      fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 1.d0/3.d0
      fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero

    case('FreeStream')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.d0 !smalldouble

      rad_e(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          iw_r_e)
      do idir = 1,ndim
        call gradientO(rad_e,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,grad_r_eO,&
           nghostcells)
        normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)**2

        ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)

      fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = dsqrt(normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_rho)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_r_e))

      !> Calculate the flux limiter, lambda
      fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = one/fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)

    case('Pomraning')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.d0 !smalldouble !*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_r_e)

      rad_e(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          iw_r_e)
      do idir = 1,ndim
        call gradientO(rad_e,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,grad_r_eO,&
           nghostcells)
        normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)**2

        ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)

      fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = dsqrt(normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_rho)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
      fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = (2.d0+fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/(6.d0+3*fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2.d0)

    case('Pomraning2')
      call mpistop("Pomraning2 is not quite working, use Pomraning or Minerbo")

      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.d0 !smalldouble

      rad_e(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          iw_r_e)
      do idir = 1,ndim
        call gradientO(rad_e,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,grad_r_eO,&
           nghostcells)
        normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)**2

        ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)

      fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = dsqrt(normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_rho)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = 1/R(coth(R)-1/R)
      fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = one/fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*(one/dtanh(fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)) - one/fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      ! fld_lambda(ixO^S) = one/fld_R(ixO^S)*(dcosh(fld_R(ixO^S))/dsinh(fld_R(ixO^S)) - one/fld_R(ixO^S))

      !>WHAT HAPPENS WHEN R=0 (full diffusion) => 1/R = NAN => dtanh(1/R) =????
      where(dtanh(fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)) .ne. dtanh(fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))
        fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = 1.d0/3.d0
      endwhere

    case('Minerbo')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.d0 !smalldouble

      rad_e(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          iw_r_e)
      do idir = 1,ndim
        call gradientO(rad_e,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,grad_r_eO,&
           nghostcells)
        normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)**2

        ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)

      fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = dsqrt(normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_rho)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Minerbo:
      do ix1=ixOmin1,ixOmax1
       do ix2=ixOmin2,ixOmax2
       do ix3=ixOmin3,ixOmax3
      
          if (fld_R(ix1,ix2,ix3) .lt. 3.d0/2.d0) then
            fld_lambda(ix1,ix2,ix3) = 2.d0/(3.d0 + dsqrt(9.d0 + &
               12.d0*fld_R(ix1,ix2,ix3)**2.d0))
          else
            fld_lambda(ix1,ix2,ix3) = 1.d0/(1.d0 + fld_R(ix1,ix2,&
               ix3) + dsqrt(1.d0 + 2.d0*fld_R(ix1,ix2,ix3)))
          endif
      enddo
      enddo
      enddo

    case('special')
      if (.not. associated(usr_special_fluxlimiter)) then
        call mpistop("special fluxlimiter not defined")
      endif
      call usr_special_fluxlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
          fld_lambda, fld_R)
    case default
      call mpistop('Fluxlimiter unknown')
    end select


    if (flux_lim_filter) then
      if (size_L_filter .lt. 1) call mpistop(&
         "D filter of size < 1 makes no sense")
      if (size_L_filter .gt. nghostcells) call &
         mpistop("D filter of size > nghostcells makes no sense")

      tmp_L(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      filtered_L(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero

      do filter = 1,size_L_filter
        do ix1 = ixOmin1+size_D_filter,ixOmax1-size_L_filter
        do ix2 = ixOmin2+size_D_filter,ixOmax2-size_L_filter
        do ix3 = ixOmin3+size_D_filter,ixOmax3-size_L_filter
          do idim = 1,ndim
            filtered_L(ix1,ix2,ix3) = filtered_L(ix1,ix2,&
               ix3) + tmp_L(ix1+filter*kr(idim,1),ix2+filter*kr(idim,2),&
               ix3+filter*kr(idim,3)) + tmp_L(ix1-filter*kr(idim,1),&
               ix2-filter*kr(idim,2),ix3-filter*kr(idim,3))
          enddo
        enddo
        enddo
        enddo
      enddo

      do ix1 = ixOmin1+size_D_filter,ixOmax1-size_D_filter
      do ix2 = ixOmin2+size_D_filter,ixOmax2-size_D_filter
      do ix3 = ixOmin3+size_D_filter,ixOmax3-size_D_filter
        tmp_L(ix1,ix2,ix3) = (tmp_L(ix1,ix2,ix3)+filtered_L(ix1,ix2,&
           ix3))/(1+2*size_L_filter*ndim)
      enddo
      enddo
      enddo

      fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = tmp_L(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    endif

  end subroutine fld_get_fluxlimiter

  !> Calculate Radiation Flux
  !> stores radiation flux in w-array
  subroutine fld_get_radflux(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, rad_flux)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out) :: rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3, 1:ndim)

    double precision :: grad_r_e(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: rad_e(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: ix1,ix2,ix3, idir

    rad_e(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, iw_r_e)

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, lambda,&
        fld_R)

    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndim
      call gradientO(rad_e,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,grad_r_eO,&
         nghostcells)
      rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir) = -(const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(kappa(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,iw_rho))*grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)

      ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
      ! rad_flux(ixO^S, idir) = -(const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,iw_rho))*grad_r_e(ixO^S)
    end do

  end subroutine fld_get_radflux

  !> Calculate Eddington-tensor
  !> Stores Eddington-tensor in w-array
  subroutine fld_get_eddington(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      eddington_tensor)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out) :: eddington_tensor(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndim,1:ndim)
    double precision :: tnsr2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndim,1:ndim)
    double precision :: normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: grad_r_e(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim), rad_e(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3, 1:ndim)
    double precision :: lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer :: i,j, idir,jdir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero

    rad_e(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, iw_r_e)
    do idir = 1,ndim
      call gradientO(rad_e,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,&
         grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir),&
         nghostcells)
      normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idir)**2

      ! call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e(ixI^S,idir))
      ! normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S,idir)**two
    end do

    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, lambda,&
        fld_R)

    !> Calculate radiation pressure
    !> P = (lambda + lambda^2 R^2)*E
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**two * fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**two
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       one/two*(one-f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) + one/two*(3.d0*f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) - one)

    

    
    do idir = 1,ndim
      eddington_tensor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,&
         idir) = half*(one-f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    enddo

    do idir = 1,ndim
      do jdir = 1,ndim
        if (idir .ne. jdir) eddington_tensor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir,jdir) = zero
        tnsr2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,&
           jdir) =  half*(3.d0*f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - 1)*grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir)*grad_r_eO(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,jdir)/normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      enddo
    enddo

    ! do idir = 1,ndim
    !   do jdir = 1,ndim
 !if (idir .ne. jdir) eddington_tensor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,jdir) = zero
 !tnsr2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,jdir) =  half*(3.d0*f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) - 1)!     *grad_r_e(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)*grad_r_e(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,jdir)/normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    !   enddo
    ! enddo

    do idir = 1,ndim
      do jdir = 1,ndim
        where ((tnsr2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,&
           jdir) .eq. tnsr2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir,jdir)) .and. (normgrad2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) .gt. smalldouble))
          eddington_tensor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir,jdir) = eddington_tensor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir,jdir) + tnsr2(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,jdir)
        endwhere
      enddo
    enddo
   

  end subroutine fld_get_eddington

  !> Calculate Radiation Pressure
  !> Returns Radiation Pressure as tensor
  subroutine fld_get_radpress(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, rad_pressure)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision             :: eddington_tensor(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndim,1:ndim)
    double precision, intent(out):: rad_pressure(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndim,1:ndim)

    integer i,j

    call fld_get_eddington(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
        eddington_tensor)

    do i=1,ndim
      do j=1,ndim
        rad_pressure(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,&
           j) = eddington_tensor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,i,j)* w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw_r_e)
      enddo
    enddo
  end subroutine fld_get_radpress

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! Multigrid diffusion
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calling all subroutines to perform the multigrid method
  !> Communicates rad_e and diff_coeff to multigrid library
  subroutine Diffuse_E_rad_mg(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_forest
    ! use mod_ghostcells_update
    use mod_multigrid_coupling
    ! use mod_physics, only: phys_set_mg_bounds, phys_req_diagonal

    type(state), target :: psa(max_blocks) !< Advance psa=psb+dtfactor*qdt*F_im(psa)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer                      :: n
    double precision             :: res, max_residual, lambda
    integer, parameter           :: max_its = 100

    integer                        :: iw_to,iw_from
    integer                        :: iigrid, igrid, id
    integer                        :: nc, lvl, i
    type(tree_node), pointer       :: pnode
    real(dp)                       :: fac, facD

    ! print*, '1'

    !> Set diffusion timestep, add previous timestep if mg did not converge:
    if (it == 0) dt_diff = 0
    dt_diff = dt_diff + qdt

    ! Avoid setting a very restrictive limit to the residual when the time step
    ! is small (as the operator is ~ 1/(D * qdt))
    if (qdt < dtmin) then
        if(mype==0)then
            print *,'skipping implicit solve: dt too small!'
            print *,'Currently at time=',global_time,' time step=',qdt,' dtmin=',dtmin
        endif
        return
    endif
    ! max_residual = 1d-7/qdt
    max_residual = fld_diff_tol !1d-7/qdt

    mg%operator_type = mg_vhelmholtz
    mg%smoother_type = mg_smoother_gs
    call mg_set_methods(mg)

    if (.not. mg%is_allocated) call mpistop(&
       "multigrid tree not allocated yet")

!   lambda = 1.d0/(dtfactor * qdt)
    lambda = 1.d0/(dtfactor * dt_diff)
    call vhelmholtz_set_lambda(lambda)

    call update_diffcoeff(psb)

    fac = 1.d0
    facD = 1.d0

    ! print*, '2'

    !This is mg_copy_to_tree from psb state
    call mg_copy_to_tree(i_diff_mg, mg_iveps, factor=facD, state_from=psb)
    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(su_, mg_irhs, factor=-lambda)
    ! iw_from=i_diff_mg
    ! iw_to=mg_iveps
    ! fac=facD
    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    pnode => igrid_to_node(igrid, mype)%node
    !    id    =  pnode%id
    !    lvl   =  mg%boxes(id)%lvl
    !    nc    =  mg%box_size_lvl(lvl)
    !    ! Include one layer of ghost cells on grid leaves
    !    {^IFONED
    !    mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
    !    }
    !    {^IFTWOD
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
    !    }
    !    {^IFTHREED
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
    !         ixMlo3-1:ixMhi3+1, iw_from)
    !    }
    ! end do

    ! print*, '3'

    !This is mg_copy_to_tree from psb state
    call mg_copy_to_tree(iw_r_e, mg_iphi, factor=fac, state_from=psb)
    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(su_, mg_irhs, factor=-lambda)
    ! iw_from=iw_r_e
    ! iw_to=mg_iphi
    ! fac=fac
    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    pnode => igrid_to_node(igrid, mype)%node
    !    id    =  pnode%id
    !    lvl   =  mg%boxes(id)%lvl
    !    nc    =  mg%box_size_lvl(lvl)
    !    ! Include one layer of ghost cells on grid leaves
    !    {^IFONED
    !    mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
    !    }
    !    {^IFTWOD
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
    !    }
    !    {^IFTHREED
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
    !         ixMlo3-1:ixMhi3+1, iw_from)
    !    }
    ! end do

    ! print*, '4'


    !>replace call set_rhs(mg, -1/dt, 0.0_dp)
    call mg_copy_to_tree(iw_r_e, mg_irhs, factor=-1/(dtfactor*dt_diff),&
        state_from=psb)
    !This is mg_copy_to_tree from psb state
    !!!  replaces::  call mg_copy_to_tree(su_, mg_irhs, factor=-lambda)
    ! iw_from=iw_r_e
    ! iw_to=mg_irhs
    ! fac=-1/(dtfactor*dt_diff)
    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    pnode => igrid_to_node(igrid, mype)%node
    !    id    =  pnode%id
    !    lvl   =  mg%boxes(id)%lvl
    !    nc    =  mg%box_size_lvl(lvl)
    !    ! Include one layer of ghost cells on grid leaves
    !    {^IFONED
    !    mg%boxes(id)%cc(0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, iw_from)
    !    }
    !    {^IFTWOD
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_from)
    !    }
    !    {^IFTHREED
    !    mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_to) = fac * &
    !         psb(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
    !         ixMlo3-1:ixMhi3+1, iw_from)
    !    }
    ! end do

    ! call phys_set_mg_bounds()

    ! print*, '5'


    if (time_advance)then
      call mg_restrict(mg, mg_iveps)
      call mg_fill_ghost_cells(mg, mg_iveps)
    endif

    ! print*, '6'


    call mg_fas_fmg(mg, .true., max_res=res)
    do n = 1, max_its
      !print*, n, res
      if (res < max_residual) exit
      call mg_fas_vcycle(mg, max_res=res)
    end do

    ! print*, '7'


    if (res .le. 0.d0) then
      if (diffcrash_resume) then
        if (mg%my_rank == 0) write(*,*) it, ' resiudal zero ', res
        return
      endif
      if (mg%my_rank == 0) then
        print*, res
        error stop "Diffusion residual to zero"
      endif
    endif

    ! print*, '8'


    if (n == max_its + 1) then
      if (diffcrash_resume) then
        if (mg%my_rank == 0) write(*,*) it, ' resiudal high ', res
        return
      endif
       if (mg%my_rank == 0) then
          print *, "Did you specify boundary conditions correctly?"
          print *, "Or is the variation in diffusion too large?"
          print*, n, res
          print *, mg%bc(1, mg_iphi)%bc_value, mg%bc(2, mg_iphi)%bc_value
       end if
       error stop "diffusion_solve: no convergence"
    end if

    ! print*, '9'


    !> Reset dt_diff when diffusion worked out
!0887 dt_diff = 0.d0
    dt_diff = 0.d0

    ! !This is mg_copy_from_tree_gc for psa state
    call mg_copy_from_tree_gc(mg_iphi, iw_r_e, state_to=psa)
    !This is mg_copy_from_tree_gc for psa state
    !!! replaces:: call mg_copy_from_tree_gc(mg_iphi, su_)
    ! iw_from=mg_iphi
    ! iw_to=iw_r_e
    ! do iigrid=1,igridstail; igrid=igrids(iigrid);
    !    pnode => igrid_to_node(igrid, mype)%node
    !    id    =  pnode%id
    !    lvl   =  mg%boxes(id)%lvl
    !    nc    =  mg%box_size_lvl(lvl)
    !    ! Include one layer of ghost cells on grid leaves
    !    {^IFONED
    !    psa(igrid)%w(ixMlo1-1:ixMhi1+1, iw_to) = &
    !         mg%boxes(id)%cc(0:nc+1, iw_from)
    !    }
    !    {^IFTWOD
    !    psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, iw_to) = &
    !         mg%boxes(id)%cc(0:nc+1, 0:nc+1, iw_from)
    !    }
    !    {^IFTHREED
    !    psa(igrid)%w(ixMlo1-1:ixMhi1+1, ixMlo2-1:ixMhi2+1, &
    !         ixMlo3-1:ixMhi3+1, iw_to) = &
    !         mg%boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iw_from)
    !    }
    ! end do

    ! print*, '10'


  end subroutine Diffuse_E_rad_mg


  !> inplace update of psa==>F_im(psa)
  subroutine Evaluate_E_rad_mg(qtC,psa)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_physics, only: phys_req_diagonal

    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: iigrid, igrid, level
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3

    call update_diffcoeff(psa)

    ixOmin1=ixGlo1+1;ixOmin2=ixGlo2+1;ixOmin3=ixGlo3+1;ixOmax1=ixGhi1-1
    ixOmax2=ixGhi2-1;ixOmax3=ixGhi3-1;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);
       ! call fld_get_diffcoef_central(psa(igrid)%w, psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L)
       call put_diffterm_onegrid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,psa(igrid)%w)
    end do
    !$OMP END PARALLEL DO

    ! enforce boundary conditions for psa
    call getbc(qtC,0.d0,psa,1,nwflux+nwaux,phys_req_diagonal)

  end subroutine Evaluate_E_rad_mg

  !> inplace update of psa==>F_im(psa)
  subroutine put_diffterm_onegrid(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)

    double precision :: gradE(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       divF(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: divF_h(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),divF_j(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: diff_term(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    integer                       :: idir, jxOmin1,jxOmin2,jxOmin3,jxOmax1,&
       jxOmax2,jxOmax3, hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3

    ! call mpistop("phys_evaluate_implicit not implemented for FLD")

    divF(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.d0
    do idir = 1,ndim
      hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
      hxOmin3=ixOmin3-kr(idir,3);hxOmax1=ixOmax1-kr(idir,1)
      hxOmax2=ixOmax2-kr(idir,2);hxOmax3=ixOmax3-kr(idir,3);
      jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
      jxOmin3=ixOmin3+kr(idir,3);jxOmax1=ixOmax1+kr(idir,1)
      jxOmax2=ixOmax2+kr(idir,2);jxOmax3=ixOmax3+kr(idir,3);

      divF_h(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         i_diff_mg)*w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
         i_diff_mg)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         i_diff_mg) + w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
         i_diff_mg))*(w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
         iw_r_e) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_r_e))
      divF_j(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         i_diff_mg)*w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
         i_diff_mg)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         i_diff_mg) + w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
         i_diff_mg))*(w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
         iw_r_e) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_r_e))
      divF(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = divF(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + 2.d0*(divF_h(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + divF_j(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/dxlevel(idir)**2
    enddo

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       iw_r_e) = divF(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine put_diffterm_onegrid


  !> Calculates cell-centered diffusion coefficient to be used in multigrid
  subroutine fld_get_diffcoef_central(w, wCT, x, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    double precision :: max_D(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        grad_r_e(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        rad_e(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: idir,i,j, ix1,ix2,ix3


    call fld_get_opacity(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)
    call fld_get_fluxlimiter(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, lambda,&
        fld_R)

    !> calculate diffusion coefficient
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_diff_mg) = (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       iw_rho))

    where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_diff_mg) .lt. 0.d0) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_diff_mg) = smalldouble

    if (diff_coef_filter) then
      !call mpistop('Hold your bloody horses, not implemented yet ')
      call fld_smooth_diffcoef(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    endif

    if (associated(usr_special_diffcoef)) call usr_special_diffcoef(w, wCT, x,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3)

  end subroutine fld_get_diffcoef_central

  subroutine update_diffcoeff(psa)
    use mod_global_parameters

    type(state), target :: psa(max_blocks)

    ! double precision :: wCT(ixG^LL,1:nw)
    integer :: iigrid, igrid, level
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3

    ixOmin1=ixGlo1+1;ixOmin2=ixGlo2+1;ixOmin3=ixGlo3+1;ixOmax1=ixGhi1-1
    ixOmax2=ixGhi2-1;ixOmax3=ixGhi3-1;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);

       ! wCT = psa(igrid)%w
        call fld_get_diffcoef_central(psa(igrid)%w, psa(igrid)%w, psa(igrid)%x,&
            ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, ixOmin1,ixOmin2,ixOmin3,&
           ixOmax1,ixOmax2,ixOmax3)
    end do
    !$OMP END PARALLEL DO

  end subroutine update_diffcoeff

  !> Use running average on Diffusion coefficient
  subroutine fld_smooth_diffcoef(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)

    double precision :: tmp_D(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        filtered_D(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: ix1,ix2,ix3, filter, idim

    if (size_D_filter .lt. 1) call mpistop(&
       "D filter of size < 1 makes no sense")
    if (size_D_filter .gt. nghostcells) call &
       mpistop("D filter of size > nghostcells makes no sense")

    tmp_D(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,i_diff_mg)
    filtered_D(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = zero

    do filter = 1,size_D_filter
      do ix1 = ixOmin1+size_D_filter,ixOmax1-size_D_filter
      do ix2 = ixOmin2+size_D_filter,ixOmax2-size_D_filter
      do ix3 = ixOmin3+size_D_filter,ixOmax3-size_D_filter
        do idim = 1,ndim
          filtered_D(ix1,ix2,ix3) = filtered_D(ix1,ix2,&
             ix3) + tmp_D(ix1+filter*kr(idim,1),ix2+filter*kr(idim,2),&
             ix3+filter*kr(idim,3)) + tmp_D(ix1-filter*kr(idim,1),&
             ix2-filter*kr(idim,2),ix3-filter*kr(idim,3))
        enddo
      enddo
      enddo
      enddo
    enddo

    do ix1 = ixOmin1+size_D_filter,ixOmax1-size_D_filter
    do ix2 = ixOmin2+size_D_filter,ixOmax2-size_D_filter
    do ix3 = ixOmin3+size_D_filter,ixOmax3-size_D_filter
      tmp_D(ix1,ix2,ix3) = (tmp_D(ix1,ix2,ix3)+filtered_D(ix1,ix2,&
         ix3))/(1+2*size_D_filter*ndim)
    enddo
    enddo
    enddo

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_diff_mg) = tmp_D(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
  end subroutine fld_smooth_diffcoef


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! Gas-Rad Energy interaction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> This subroutine calculates the radiation heating, radiation cooling
  !> and photon tiring using an implicit scheme.
  !> These sourceterms are applied using the root-finding of a 4th order polynomial
  !> This routine loops over every cell in the domain
  !> and computes the coefficients of the polynomials in every cell
  subroutine Energy_interaction(w, wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    use mod_global_parameters
    use mod_geometry
    use mod_physics
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: a1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: c0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        c1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: e_gas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        E_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: sigma_b

    integer :: i,j,idir,ix1,ix2,ix3

    !> e_gas is the INTERNAL ENERGY without KINETIC ENERGY
    ! if (.not. block%e_is_internal) then
      e_gas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,iw_e) - half*sum(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, iw_mom(:))**2, dim=ndim+1)/wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, iw_rho)
    ! else
    !   e_gas(ixO^S) = wCT(ixO^S,iw_e)
    ! endif

    do ix1=ixOmin1,ixOmax1
     do ix2=ixOmin2,ixOmax2
     do ix3=ixOmin3,ixOmax3
    
      e_gas(ix1,ix2,ix3) = max(e_gas(ix1,ix2,ix3),&
         small_pressure/(fld_gamma-1.d0))
    enddo
    enddo
    enddo

    E_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       iw_r_e)

    if (associated(usr_special_opacity_qdot)) then
      call usr_special_opacity_qdot(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,kappa)
    else
      call fld_get_opacity(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)
    endif

    sigma_b = const_rad_a*const_c/4.d0*(unit_temperature**&
       4.d0)/(unit_velocity*unit_pressure)

    if (fld_interaction_method .eq. 'Instant') then
      a1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = const_rad_a*(fld_mu*const_mp/const_kB*(fld_gamma-&
         1))**4/(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_rho)*unit_density)**4 /unit_pressure**3
      a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = e_gas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + E_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)

      c0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = a2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/a1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      c1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1.d0/a1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      !> Calculate coefficients for polynomial
      a1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 4*kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*sigma_b*(fld_gamma-one)**4/wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_rho)**3*dt
      a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = (const_c/unit_velocity)*kappa(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,iw_rho)*dt

      c0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = ((one + a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*e_gas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*E_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/a1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      c1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = (one + a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/a1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    endif

    !> Loop over every cell for rootfinding method
    do ix1=ixOmin1,ixOmax1
     do ix2=ixOmin2,ixOmax2
     do ix3=ixOmin3,ixOmax3
    
      select case(fld_interaction_method)
      case('Bisect')
        call Bisection_method(e_gas(ix1,ix2,ix3), E_rad(ix1,ix2,ix3), c0(ix1,&
           ix2,ix3), c1(ix1,ix2,ix3))
      case('Newton')
        call Newton_method(e_gas(ix1,ix2,ix3), E_rad(ix1,ix2,ix3), c0(ix1,ix2,&
           ix3), c1(ix1,ix2,ix3))
      case('Halley')
        call Halley_method(e_gas(ix1,ix2,ix3), E_rad(ix1,ix2,ix3), c0(ix1,ix2,&
           ix3), c1(ix1,ix2,ix3))
      case('Instant')
        call Halley_method(e_gas(ix1,ix2,ix3), E_rad(ix1,ix2,ix3), c0(ix1,ix2,&
           ix3), c1(ix1,ix2,ix3))
      case default
        call mpistop('root-method not known')
      end select
    enddo
    enddo
    enddo

    if (fld_interaction_method .eq. 'Instant') then
      E_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) - e_gas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      !> advance E_rad
      E_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = (a1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*e_gas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**4.d0 + E_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))/(one + a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
    endif

    !> new w = w + dt f(wCT)
    !> e_gas,E_rad = wCT + dt f(wCT)
    !> dt f(wCT) = e_gas,E_rad - wCT
    !> new w = w +  e_gas,E_rad - wCT

    !> WAIT A SECOND?! DOES THIS EVEN WORK WITH THESE KIND OF IMPLICIT METHODS?
    !> NOT QUITE SURE TO BE HONEST
    !> IS IT POSSIBLE TO SHUT DOWN SOURCE SPLITTING FOR JUST THIS TERM?
    !> FIX BY PERFORMING Energy_interaction on (w,w,...)

    ! call mpistop('This still has to be fixed somehow')

    !> Update gas-energy in w, internal + kinetic
    ! w(ixO^S,iw_e) = w(ixO^S,iw_e) + e_gas(ixO^S) - wCT(ixO^S,iw_e)
    ! {do ix^D=ixOmin^D,ixOmax^D\ }
    !   e_gas(ix^D) = max(e_gas(ix^D),small_pressure/(fld_gamma-1.d0))
    ! {enddo\}

    !{do ix^D=ixOmin^D,ixOmax^D\ }
    !  w(ix^D,iw_e) = max(e_gas(ix^D),small_pressure/(fld_gamma - 1))
    !{enddo\}
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       iw_e) = e_gas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    !> Beginning of module substracted wCT Ekin
    !> So now add wCT Ekin
    ! if (.not. block%e_is_internal) then
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       iw_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       iw_e) + half*sum(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        iw_mom(:))**2, dim=ndim+1)/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3, iw_rho)
    ! else
    !   w(ixO^S,iw_e) = w(ixO^S,iw_e)
    ! endif

    !> Update rad-energy in w
    ! w(ixO^S,iw_r_e) = w(ixO^S,iw_r_e) + E_rad(ixO^S) - wCT(ixO^S,iw_r_e)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       iw_r_e) = E_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine Energy_interaction


  !> Find the root of the 4th degree polynomial using the bisection method
  subroutine Bisection_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: bisect_a, bisect_b, bisect_c
    integer :: n, max_its

    n = 0
    max_its = 100

    bisect_a = zero
    bisect_b = min(abs(c0/c1),abs(c0)**(1.d0/4.d0))+smalldouble

    ! do while (abs(Polynomial_Bisection(bisect_b, c0, c1)-Polynomial_Bisection(bisect_a, c0, c1))&
    !    .ge. fld_bisect_tol*min(e_gas,E_rad))
    do while (abs(bisect_b-bisect_a) .ge. fld_bisect_tol*min(e_gas,E_rad))
      bisect_c = (bisect_a + bisect_b)/two

      n = n +1
      if (n .gt. max_its) then
        goto 2435
        call mpistop('No convergece in bisection scheme')
      endif

      if (Polynomial_Bisection(bisect_a, c0, c1)*Polynomial_Bisection(bisect_b,&
          c0, c1) .lt. zero) then

        if (Polynomial_Bisection(bisect_a, c0,&
            c1)*Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_b = bisect_c
        elseif (Polynomial_Bisection(bisect_b, c0,&
            c1)*Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_a = bisect_c
        elseif (Polynomial_Bisection(bisect_a, c0, c1) .eq. zero) then
          bisect_b = bisect_a
          bisect_c = bisect_a
          goto 2435
        elseif (Polynomial_Bisection(bisect_b, c0, c1) .eq. zero) then
          bisect_a = bisect_b
          bisect_c = bisect_b
          goto 2435
        elseif (Polynomial_Bisection(bisect_c, c0, c1) .eq. zero) then
          bisect_a = bisect_c
          bisect_b = bisect_c
          goto 2435
        else
          call mpistop("Problem with fld bisection method")
        endif
      elseif (Polynomial_Bisection(bisect_a, c0,&
          c1) - Polynomial_Bisection(bisect_b, c0,&
          c1) .lt. fld_bisect_tol*Polynomial_Bisection(bisect_a, c0, c1)) then
        goto 2435
      else
        bisect_a = e_gas
        bisect_b = e_gas
        print*, "IGNORING GAS-RAD ENERGY EXCHANGE ", c0, c1

        print*, Polynomial_Bisection(bisect_a, c0, c1), Polynomial_Bisection(bisect_b, c0, c1)

        if (Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_b = bisect_a
        elseif (Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_a = bisect_b
        endif

        goto 2435

      endif
    enddo

      2435 e_gas = (bisect_a + bisect_b)/two
  end subroutine Bisection_method

  !> Find the root of the 4th degree polynomial using the Newton-Ralphson method
  subroutine Newton_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: xval, yval, der, deltax

    integer :: ii

    yval = bigdouble
    xval = e_gas
    der = one
    deltax = one

    ii = 0
    !> Compare error with dx = dx/dy dy
    do while (abs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der = dPolynomial_Bisection(xval, c0, c1)
      deltax = -yval/der
      xval = xval + deltax
      ii = ii + 1
      if (ii .gt. 1d3) then
        print*, 'skip to bisection algorithm'
        call Bisection_method(e_gas, E_rad, c0, c1)
        return
      endif
    enddo

    e_gas = xval
  end subroutine Newton_method

  !> Find the root of the 4th degree polynomial using the Halley method
  subroutine Halley_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: xval, yval, der, dder, deltax

    integer :: ii

    yval = bigdouble
    xval = e_gas
    der = one
    dder = one
    deltax = one

    ii = 0
    !> Compare error with dx = dx/dy dy
    do while (abs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der = dPolynomial_Bisection(xval, c0, c1)
      dder = ddPolynomial_Bisection(xval, c0, c1)
      deltax = -two*yval*der/(two*der**2 - yval*dder)
      xval = xval + deltax
      ii = ii + 1
      if (ii .gt. 1d3) then
        ! call mpistop('Halley did not convergggge')
        call Newton_method(e_gas, E_rad, c0, c1)
        return
      endif
    enddo

    e_gas = xval
  end subroutine Halley_method

  !> Evaluate polynomial at argument e_gas
  function Polynomial_Bisection(e_gas, c0, c1) result(val)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: val

    val = e_gas**4.d0 + c1*e_gas - c0
  end function Polynomial_Bisection

  !> Evaluate first derivative of polynomial at argument e_gas
  function dPolynomial_Bisection(e_gas, c0, c1) result(der)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: der

    der = 4.d0*e_gas**3.d0 + c1
  end function dPolynomial_Bisection

  !> Evaluate second derivative of polynomial at argument e_gas
  function ddPolynomial_Bisection(e_gas, c0, c1) result(dder)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: dder

    dder = 4.d0*3.d0*e_gas**2.d0
  end function ddPolynomial_Bisection

  !> Calculate gradient of a scalar q within ixL in direction idir
  !> difference with gradient is gradq(ixO^S), NOT gradq(ixI^S)
  subroutine gradientO(q,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,gradq,n)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
    integer, intent(in)             :: n

    double precision, intent(in)    :: q(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out)   :: gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,&
       jxOmax3, hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3

    ! hxO^L=ixO^L-n*kr(idir,^D);
    ! jxO^L=ixO^L+n*kr(idir,^D);
    !
    ! if (n .gt. nghostcells) call mpistop("gradientO stencil too wide")
    !
    ! gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(2*n*dxlevel(idir))
    ! gradq(ixO^S)=(q(jxO^S)-q(hxO^S))/(x(jxO^S,idir)-x(hxO^S,idir))

    !> Using higher order derivatives with wider stencil according to:
    !> https://en.wikipedia.org/wiki/Finite_difference_coefficient

     if (n .gt. nghostcells) then
       call mpistop("gradientO stencil too wide")
     elseif (n .eq. 1) then
       hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
       hxOmin3=ixOmin3-kr(idir,3);hxOmax1=ixOmax1-kr(idir,1)
       hxOmax2=ixOmax2-kr(idir,2);hxOmax3=ixOmax3-kr(idir,3);
       jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
       jxOmin3=ixOmin3+kr(idir,3);jxOmax1=ixOmax1+kr(idir,1)
       jxOmax2=ixOmax2+kr(idir,2);jxOmax3=ixOmax3+kr(idir,3);
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
          jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
          hxOmin3:hxOmax3))/(x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
          idir)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,idir))
     elseif (n .eq. 2) then
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.d0
       !> coef 2/3
       hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
       hxOmin3=ixOmin3-kr(idir,3);hxOmax1=ixOmax1-kr(idir,1)
       hxOmax2=ixOmax2-kr(idir,2);hxOmax3=ixOmax3-kr(idir,3);
       jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
       jxOmin3=ixOmin3+kr(idir,3);jxOmax1=ixOmax1+kr(idir,1)
       jxOmax2=ixOmax2+kr(idir,2);jxOmax3=ixOmax3+kr(idir,3);
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) + 2.d0/3.d0*(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
          jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))
       !> coef -1/12
       hxOmin1=ixOmin1-2*kr(idir,1);hxOmin2=ixOmin2-2*kr(idir,2)
       hxOmin3=ixOmin3-2*kr(idir,3);hxOmax1=ixOmax1-2*kr(idir,1)
       hxOmax2=ixOmax2-2*kr(idir,2);hxOmax3=ixOmax3-2*kr(idir,3);
       jxOmin1=ixOmin1+2*kr(idir,1);jxOmin2=ixOmin2+2*kr(idir,2)
       jxOmin3=ixOmin3+2*kr(idir,3);jxOmax1=ixOmax1+2*kr(idir,1)
       jxOmax2=ixOmax2+2*kr(idir,2);jxOmax3=ixOmax3+2*kr(idir,3);
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) - 1.d0/12.d0*(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
          jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))
       !> divide by dx
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)/dxlevel(idir)
     elseif (n .eq. 3) then
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.d0
       !> coef 3/4
       hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
       hxOmin3=ixOmin3-kr(idir,3);hxOmax1=ixOmax1-kr(idir,1)
       hxOmax2=ixOmax2-kr(idir,2);hxOmax3=ixOmax3-kr(idir,3);
       jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
       jxOmin3=ixOmin3+kr(idir,3);jxOmax1=ixOmax1+kr(idir,1)
       jxOmax2=ixOmax2+kr(idir,2);jxOmax3=ixOmax3+kr(idir,3);
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) + 3.d0/4.d0*(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
          jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))
       !> coef -3/20
       hxOmin1=ixOmin1-2*kr(idir,1);hxOmin2=ixOmin2-2*kr(idir,2)
       hxOmin3=ixOmin3-2*kr(idir,3);hxOmax1=ixOmax1-2*kr(idir,1)
       hxOmax2=ixOmax2-2*kr(idir,2);hxOmax3=ixOmax3-2*kr(idir,3);
       jxOmin1=ixOmin1+2*kr(idir,1);jxOmin2=ixOmin2+2*kr(idir,2)
       jxOmin3=ixOmin3+2*kr(idir,3);jxOmax1=ixOmax1+2*kr(idir,1)
       jxOmax2=ixOmax2+2*kr(idir,2);jxOmax3=ixOmax3+2*kr(idir,3);
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) - 3.d0/20.d0*(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
          jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))
       !> coef 1/60
       hxOmin1=ixOmin1-3*kr(idir,1);hxOmin2=ixOmin2-3*kr(idir,2)
       hxOmin3=ixOmin3-3*kr(idir,3);hxOmax1=ixOmax1-3*kr(idir,1)
       hxOmax2=ixOmax2-3*kr(idir,2);hxOmax3=ixOmax3-3*kr(idir,3);
       jxOmin1=ixOmin1+3*kr(idir,1);jxOmin2=ixOmin2+3*kr(idir,2)
       jxOmin3=ixOmin3+3*kr(idir,3);jxOmax1=ixOmax1+3*kr(idir,1)
       jxOmax2=ixOmax2+3*kr(idir,2);jxOmax3=ixOmax3+3*kr(idir,3);
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) + 1.d0/60.d0*(q(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
          jxOmin3:jxOmax3)-q(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3))
       !> divide by dx
       gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = gradq(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)/dxlevel(idir)
     else
       call mpistop("gradientO stencil unknown")
     endif

  end subroutine gradientO

end module mod_fld
