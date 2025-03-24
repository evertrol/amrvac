!> Module to include CAK radiation line force in (magneto)hydrodynamic models
!> Computes both the force from free electrons and the force from an ensemble of
!> lines (various possibilities for the latter).
!> There is an option to only simulate the pure radial CAK force (with various
!> corrections applied) as well as the full vector CAK force. Depending on the
!> chosen option additional output are the CAK line force component(s) and,
!> when doing a 1-D radial force, the finite disc factor.
!>
!> USAGE:
!>
!>  1. Include a cak_list in the .par file and activate (m)hd_cak_force in the
!>     (m)hd_list
!>  2. Create a mod_usr.t file for the problem with appropriate initial and
!>     boundary conditions
!>  3. In the mod_usr.t header call the mod_cak_force module to have access to
!>     global variables from mod_cak_force, which may be handy for printing or
!>     the computation of other variables inside mod_usr.t
!>  4. In usr_init of mod_usr.t call the set_cak_force_norm routine and pass
!>     along the stellar radius and wind temperature---this is needed to
!>     correctly compute the (initial) force normalisation inside mod_cak_force
!>  5. Ensure that the order of calls in usr_init is similar as for test problem
!>     CAKwind_spherical_1D: first reading usr.par list; then set unit scales;
!>     then call (M)HD_activate; then call set_cak_force_norm. This order avoids
!>     an incorrect force normalisation and code crash
!>
!> Developed by Florian Driessen (2022)
module mod_cak_force
  use mod_physics, only: phys_get_pthermal, physics_type
  implicit none

  !> Line-ensemble parameters in the Gayley (1995) formalism
  real(8), public :: cak_alpha, gayley_qbar, gayley_q0

  !> Switch to choose between the 1-D CAK line force options
  integer :: cak_1d_opt

  ! Avoid magic numbers in code for 1-D CAK line force option
  integer, parameter, private :: radstream=0, fdisc=1, fdisc_cutoff=2

  !> To treat source term in split or unsplit (default) fashion
  logical :: cak_split=.false.

  !> To activate the original CAK 1-D line force computation
  logical :: cak_1d_force=.false.

  !> To activate the vector CAK line force computation
  logical :: cak_vector_force=.false.

  !> To activate the pure radial vector CAK line force computation
  logical :: fix_vector_force_1d=.false.

  !> Amount of rays in radiation polar and radiation azimuthal direction
  integer :: nthetaray, nphiray

  !> Ray positions + weights for impact parameter and azimuthal radiation angle
  real(8), allocatable, private :: ay(:), wy(:), aphi(:), wphi(:)

  !> The adiabatic index
  real(8), private :: cak_gamma

  !> Variables needed to compute force normalisation fnorm in initialisation
  real(8), private :: lum, dlum, drstar, dke, dclight

  !> To enforce a floor temperature when doing adiabatic (M)HD
  real(8), private :: tfloor

  !> Extra slots to store quantities in w-array
  integer :: gcak1_, gcak2_, gcak3_, fdf_

  !> Public method
  public :: set_cak_force_norm
  
contains

  !> Read this module's parameters from a file
  subroutine cak_params_read(files)
    use mod_global_parameters, only: unitpar

    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n

    namelist /cak_list/ cak_alpha, gayley_qbar, gayley_q0, cak_1d_opt,&
        cak_split, cak_1d_force, cak_vector_force, nphiray, nthetaray,&
        fix_vector_force_1d

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, cak_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine cak_params_read

  !> Initialize the module
  subroutine cak_init(phys_gamma)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    real(8), intent(in) :: phys_gamma

    cak_gamma = phys_gamma

    ! Set some defaults when user does not
    cak_alpha   = 0.65d0
    gayley_qbar = 2000.0d0
    gayley_q0   = 2000.0d0
    cak_1d_opt  = fdisc
    nthetaray   = 6
    nphiray     = 6

    call cak_params_read(par_files)

    if (cak_1d_force) then
      gcak1_ = var_set_extravar("gcak1", "gcak1")
      fdf_   = var_set_extravar("fdfac", "fdfac")
    endif

    if (cak_vector_force) then
      gcak1_ = var_set_extravar("gcak1", "gcak1")
      gcak2_ = var_set_extravar("gcak2", "gcak2")
      gcak3_ = var_set_extravar("gcak3", "gcak3")
      call rays_init(nthetaray,nphiray)
    endif

    ! Some sanity checks
    if ((cak_alpha <= 0.0d0) .or. (cak_alpha > 1.0d0)) then
      call mpistop('CAK error: choose alpha in [0,1[')
    endif

    if ((gayley_qbar < 0.0d0) .or. (gayley_q0 < 0.0d0)) then
      call mpistop('CAK error: chosen Qbar or Q0 is < 0')
    endif

    if (cak_1d_force .and. cak_vector_force) then
      call mpistop('CAK error: choose either 1-D or vector force')
    endif

  end subroutine cak_init

  !> Compute some (unitless) variables for CAK force normalisation
  subroutine set_cak_force_norm(rstar,twind)
    use mod_global_parameters
    use mod_constants

    real(8), intent(in) :: rstar, twind

    lum     = 4.0d0*dpi * rstar**2.0d0 * const_sigma * twind**4.0d0
    dke     = const_kappae * unit_density * unit_length
    dclight = const_c/unit_velocity
    dlum    = lum/(unit_density * unit_length**5.0d0 / unit_time**3.0d0)
    drstar  = rstar/unit_length
    tfloor  = twind/unit_temperature

  end subroutine set_cak_force_norm
  
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine cak_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,energy,&
     qsourcesplit,active)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in)    :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    real(8), intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    real(8), intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in)    :: energy, qsourcesplit
    logical, intent(inout) :: active

    ! Local variables
    integer :: idir
    real(8) :: gl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3),&
        ge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        ptherm(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        pmin(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    ! By default add source in unsplit fashion together with the fluxes
    if (qsourcesplit .eqv. cak_split) then

      ! Thomson force
      call get_gelectron(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,x,ge)

      ! CAK line force
      gl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3) = 0.0d0

      if (cak_1d_force) then
        call get_cak_force_radial(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,gl)
      elseif (cak_vector_force) then
        call get_cak_force_vector(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,gl)
      else
        call mpistop("No valid force option")
      endif

      ! Update conservative vars: w = w + qdt*gsource
      do idir = 1,ndir
        if (idir == 1) gl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir) = gl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir) + ge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mom(idir)) + qdt * gl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw_rho)
                                
        if (energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_e) + qdt * gl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_mom(idir))
        endif
      enddo

      ! Impose fixed floor temperature to mimic stellar heating
      if (energy) then
        call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ptherm)
        pmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,iw_rho) * tfloor

        where (ptherm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) < pmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_e) + (pmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) - ptherm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))/(cak_gamma - 1.0d0)
        endwhere
      endif
    endif

  end subroutine cak_add_source

  !> 1-D CAK line force in the Gayley line-ensemble distribution parametrisation
  subroutine get_cak_force_radial(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,gcak)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in)    :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    real(8), intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    real(8), intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    real(8), intent(inout) :: gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1:3)
  
    ! Local variables
    real(8) :: vr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        dvrdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    real(8) :: beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        fdfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        taus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        ge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
  
    vr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,iw_mom(1)) / wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,iw_rho)
    call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vr,x,1,dvrdr)

    if (physics_type == 'hd') then
      ! Monotonic flow to avoid multiple resonances and radiative coupling
      dvrdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = abs(dvrdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
    else
      ! Allow material to fallback to the star in a magnetosphere model
      dvrdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = max(dvrdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3), smalldouble)
    endif
  
    ! Thomson force
    call get_gelectron(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,x,ge)

    ! Sobolev optical depth for line ensemble (tau = Qbar * t_r) and the force
    select case (cak_1d_opt)
    case(radstream, fdisc)
      taus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)   = gayley_qbar * dke * dclight * &
         wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         iw_rho)/dvrdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1) = gayley_qbar/(1.0d0 - cak_alpha) * ge(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/taus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**cak_alpha

    case(fdisc_cutoff)
      taus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)   = gayley_q0 * dke * dclight * wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_rho)/dvrdr(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1) = gayley_qbar * ge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)                                  * ( (1.0d0 + &
         taus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))**(1.0d0 - cak_alpha) - 1.0d0 ) / ( (1.0d0 - &
         cak_alpha) * taus(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) )
    case default
      call mpistop("Error in force computation.")
    end select

    ! Finite disk factor parameterisation (Owocki & Puls 1996)
    beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = ( 1.0d0 - vr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) * dvrdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) ) * (drstar/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1))**2.0d0

    select case (cak_1d_opt)
    case(radstream)
      fdfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 1.0d0
    case(fdisc, fdisc_cutoff)
      where (beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) >= 1.0d0)
        fdfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = 1.0d0/(1.0d0 + cak_alpha)
      elsewhere (beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) < -1.0d10)
        fdfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = abs(beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))**cak_alpha / (1.0d0 + cak_alpha)
      elsewhere (abs(beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)) > 1.0d-3)
        fdfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = (1.0d0 - (1.0d0 - beta_fd(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3))**(1.0d0 + cak_alpha)) / &
           (beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(1.0d0 + cak_alpha))
      elsewhere
        fdfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = 1.0d0 - 0.5d0*cak_alpha*beta_fd(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3) * (1.0d0 + 1.0d0/3.0d0 * (1.0d0 - &
           cak_alpha)*beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
      endwhere
    end select

    ! Correct radial line force for finite disc (if applicable)
    gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) = gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1) * fdfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) = 0.0d0
    gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) = 0.0d0
      
    ! Fill the nwextra slots for output
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       gcak1_) = gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       fdf_)   = fdfac(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
  end subroutine get_cak_force_radial

  !> Vector CAK line force in the Gayley line-ensemble distribution parametrisation
  subroutine get_cak_force_vector(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,gcak)
    use mod_global_parameters
    use mod_usr_methods

    ! Subroutine arguments
    integer, intent(in)    :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    real(8), intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    real(8), intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    real(8), intent(inout) :: gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1:3)

    ! Local variables
    integer :: ix1,ix2,ix3, itray, ipray
    real(8) :: a1, a2, a3, wyray, y, wpray, phiray, wtot, mustar, dvndn
    real(8) :: costp, costp2, sintp, cospp, sinpp, cott0
    real(8) :: vr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        vt(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        vp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real(8) :: vrr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        vtr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        vpr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    real(8) :: dvrdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        dvtdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        dvpdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    real(8) :: dvrdt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        dvtdt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        dvpdt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    real(8) :: dvrdp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        dvtdp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        dvpdp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    
    ! Initialisation to have full velocity strain tensor expression at all times
    vt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    vtr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    vp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    vpr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    cott0 = 0.0d0
    dvrdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    dvtdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    dvpdr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    dvrdt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    dvtdt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    dvpdt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    dvrdp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    dvtdp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    dvpdp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0

    ! Populate velocity field(s) depending on dimensions and directions
    vr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)  = wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,iw_mom(1)) / wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,iw_rho)
    vrr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = vr(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3) / x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1)

    
    vt(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)  = wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,iw_mom(2)) / wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,iw_rho)
    vtr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = vt(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3) / x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1)
    
    if (ndir > 2) then
      vp(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)  = wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,iw_mom(3)) / wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,iw_rho)
      vpr(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) = vp(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) / x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1)
    endif
   
    
    ! Derivatives of velocity field in each coordinate direction (r=1,t=2,p=3)
    call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vr,x,1,dvrdr)
    
    
    call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vr,x,2,dvrdt)
    call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vt,x,1,dvtdr)
    call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vt,x,2,dvtdt)

    if (ndir > 2) then
      call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vp,x,1,dvpdr)
      call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vp,x,2,dvpdt)
    endif
   
    
    call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vr,x,3,dvrdp)
    call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vt,x,3,dvtdp)
    call get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vp,x,3,dvpdp)
   

    ! Get total acceleration from all rays at a certain grid point
    do ix3=ixOmin3,ixOmax3
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
      ! Loop over the rays; first theta then phi radiation angle
      ! Get weights from current ray and their position
      do itray = 1,nthetaray
        wyray  = wy(itray)
        y      = ay(itray)

        do ipray = 1,nphiray
          wpray = wphi(ipray)
          phiray  = aphi(ipray)

          ! Redistribute the phi rays by a small offset
          ! if (mod(itp,3) == 1) then
          !   phip = phip + dphi/3.0d0
          ! elseif (mod(itp,3) == 2) then
          !   phip = phip - dphi/3.0d0
          ! endif
          
          ! === Geometrical factors ===
          ! Make y quadrature linear in mu, not mu**2; better for gtheta,gphi
          ! y -> mu quadrature is preserved; y=0 <=> mu=1; y=1 <=> mu=mustar
          mustar = sqrt(max(1.0d0 - (drstar/x(ix1,ix2,ix3,1))**2.0d0, 0.0d0))
          costp  = 1.0d0 - y*(1.0d0 - mustar)
          costp2 = costp*costp
          sintp  = sqrt(max(1.0d0 - costp2, 0.0d0))
          sinpp  = sin(phiray)
          cospp  = cos(phiray)
           cott0  = cos(x(ix1,ix2,ix3,2))/sin(x(ix1,ix2,ix3,2))

          ! More weight close to star, less farther away
          wtot  = wyray * wpray * (1.0d0 - mustar)

          ! Convenients a la Cranmer & Owocki (1995)
          a1 = costp
          a2 = sintp * cospp
          a3 = sintp * sinpp

          ! Get total velocity gradient for one ray with given (theta', phi')
          dvndn = a1*a1 * dvrdr(ix1,ix2,ix3) + a2*a2 * (dvtdt(ix1,ix2,&
             ix3) + vrr(ix1,ix2,ix3))  + a3*a3 * (dvpdp(ix1,ix2,&
             ix3) + cott0 * vtr(ix1,ix2,ix3) + vrr(ix1,ix2,&
             ix3))   + a1*a2 * (dvtdr(ix1,ix2,ix3) + dvrdt(ix1,ix2,&
             ix3) - vtr(ix1,ix2,ix3))         + a1*a3 * (dvpdr(ix1,ix2,&
             ix3) + dvrdp(ix1,ix2,ix3) - vpr(ix1,ix2,&
             ix3))         + a2*a3 * (dvpdt(ix1,ix2,ix3) + dvtdp(ix1,ix2,&
             ix3) - cott0 * vpr(ix1,ix2,ix3))

          ! No multiple resonances in CAK
          dvndn = abs(dvndn)

          ! Convert gradient back from wind coordinates (r',theta',phi') to
          ! stellar coordinates (r,theta,phi)
          gcak(ix1,ix2,ix3,1) = gcak(ix1,ix2,ix3,1) + (dvndn/wCT(ix1,ix2,ix3,&
             iw_rho))**cak_alpha * a1 * wtot
          gcak(ix1,ix2,ix3,2) = gcak(ix1,ix2,ix3,2) + (dvndn/wCT(ix1,ix2,ix3,&
             iw_rho))**cak_alpha * a2 * wtot
          gcak(ix1,ix2,ix3,3) = gcak(ix1,ix2,ix3,3) + (dvndn/wCT(ix1,ix2,ix3,&
             iw_rho))**cak_alpha * a3 * wtot
        enddo
      enddo
    enddo
    enddo
    enddo

    ! Normalisation for line force
    ! NOTE: extra 1/pi factor comes from integration in radiation Phi angle
    gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       :) = (dke*gayley_qbar)**(1.0d0 - cak_alpha)/(1.0d0 - cak_alpha)    * &
       dlum/(4.0d0*dpi*drstar**2.0d0 * dclight**(1.0d0+cak_alpha)) * &
       gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)/dpi

    if (fix_vector_force_1d) then
      gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2) = 0.0d0
      gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3) = 0.0d0
    endif
          
    ! Fill the nwextra slots for output
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       gcak1_) = gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       gcak2_) = gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       gcak3_) = gcak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)

  end subroutine get_cak_force_vector
  
  !> Compute continuum radiation force from Thomson scattering
  subroutine get_gelectron(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,ge)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    real(8), intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    real(8), intent(out):: ge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = dke * &
       dlum/(4.0d0*dpi * dclight * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)**2.0d0)

  end subroutine get_gelectron

  !> Check time step for total radiation contribution
  subroutine cak_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    use mod_global_parameters

    integer, intent(in)    :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    real(8), intent(in)    :: dx1,dx2,dx3, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    real(8), intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    real(8), intent(inout) :: dtnew
    
    ! Local variables
    real(8) :: ge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3), max_gr,&
        dt_cak

    call get_gelectron(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,ge)

    dtnew = bigdouble

    ! Get dt from line force that is saved in the w-array in nwextra slot
    max_gr = max( maxval(abs(ge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       gcak1_))), epsilon(1.0d0) )
    dt_cak = minval( sqrt(block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)/max_gr) )
    dtnew  = min(dtnew, courantpar*dt_cak)

    
    if (cak_vector_force) then
      max_gr = max( maxval(abs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,gcak2_))), epsilon(1.0d0) )
      dt_cak = minval( sqrt(block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1) * block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)/max_gr) )
      dtnew  = min(dtnew, courantpar*dt_cak)

      
      max_gr = max( maxval(abs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,gcak3_))), epsilon(1.0d0) )
      dt_cak = minval( sqrt(block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1) * sin(block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3))/max_gr) )
      dtnew  = min(dtnew, courantpar*dt_cak)
     
    endif
   

  end subroutine cak_get_dt

  !> Compute velocity gradient in direction 'idir' on a non-uniform grid
  subroutine get_velocity_gradient(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vfield,x,idir,&
     grad_vn)
    use mod_global_parameters

    integer, intent(in)  :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
    real(8), intent(in)  :: vfield(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    real(8), intent(out) :: grad_vn(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    ! Local variables
    real(8) :: forw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        backw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        cent(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: jrxmin1,jrxmin2,jrxmin3,jrxmax1,jrxmax2,jrxmax3, hrxmin1,&
       hrxmin2,hrxmin3,hrxmax1,hrxmax2,hrxmax3,jtxmin1,jtxmin2,jtxmin3,jtxmax1,&
       jtxmax2,jtxmax3, htxmin1,htxmin2,htxmin3,htxmax1,htxmax2,htxmax3,&
       jpxmin1,jpxmin2,jpxmin3,jpxmax1,jpxmax2,jpxmax3, hpxmin1,hpxmin2,&
       hpxmin3,hpxmax1,hpxmax2,hpxmax3

    ! Index +1 (j) and index -1 (h) in radial direction; kr(dir,dim)=1, dir=dim
    jrxmin1=ixOmin1+kr(1,1);jrxmin2=ixOmin2+kr(1,2);jrxmin3=ixOmin3+kr(1,3)
    jrxmax1=ixOmax1+kr(1,1);jrxmax2=ixOmax2+kr(1,2);jrxmax3=ixOmax3+kr(1,3);
    hrxmin1=ixOmin1-kr(1,1);hrxmin2=ixOmin2-kr(1,2);hrxmin3=ixOmin3-kr(1,3)
    hrxmax1=ixOmax1-kr(1,1);hrxmax2=ixOmax2-kr(1,2);hrxmax3=ixOmax3-kr(1,3);

    
    ! Index +1 (j) and index -1 (h) in polar direction
    jtxmin1=ixOmin1+kr(2,1);jtxmin2=ixOmin2+kr(2,2);jtxmin3=ixOmin3+kr(2,3)
    jtxmax1=ixOmax1+kr(2,1);jtxmax2=ixOmax2+kr(2,2);jtxmax3=ixOmax3+kr(2,3);
    htxmin1=ixOmin1-kr(2,1);htxmin2=ixOmin2-kr(2,2);htxmin3=ixOmin3-kr(2,3)
    htxmax1=ixOmax1-kr(2,1);htxmax2=ixOmax2-kr(2,2);htxmax3=ixOmax3-kr(2,3);
   

    
    ! Index +1 (j) and index -1 (h) in azimuthal direction
    jpxmin1=ixOmin1+kr(3,1);jpxmin2=ixOmin2+kr(3,2);jpxmin3=ixOmin3+kr(3,3)
    jpxmax1=ixOmax1+kr(3,1);jpxmax2=ixOmax2+kr(3,2);jpxmax3=ixOmax3+kr(3,3);
    hpxmin1=ixOmin1-kr(3,1);hpxmin2=ixOmin2-kr(3,2);hpxmin3=ixOmin3-kr(3,3)
    hpxmax1=ixOmax1-kr(3,1);hpxmax2=ixOmax2-kr(3,2);hpxmax3=ixOmax3-kr(3,3);
   

    ! grad(v.n) on non-uniform grid according to Sundqvist & Veronis (1970)
    select case (idir)
    case(1) ! Radial forward, backward, and central derivatives
      forw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)  = (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1) - x(hrxmin1:hrxmax1,hrxmin2:hrxmax2,&
         hrxmin3:hrxmax3,1)) * vfield(jrxmin1:jrxmax1,jrxmin2:jrxmax2,&
         jrxmin3:jrxmax3) / ((x(jrxmin1:jrxmax1,jrxmin2:jrxmax2,&
         jrxmin3:jrxmax3,1) - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)) * (x(jrxmin1:jrxmax1,jrxmin2:jrxmax2,&
         jrxmin3:jrxmax3,1) - x(hrxmin1:hrxmax1,hrxmin2:hrxmax2,&
         hrxmin3:hrxmax3,1)))

      backw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = -(x(jrxmin1:jrxmax1,jrxmin2:jrxmax2,&
         jrxmin3:jrxmax3,1) - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)) * vfield(hrxmin1:hrxmax1,hrxmin2:hrxmax2,&
         hrxmin3:hrxmax3) / ((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1) - x(hrxmin1:hrxmax1,hrxmin2:hrxmax2,&
         hrxmin3:hrxmax3,1)) * (x(jrxmin1:jrxmax1,jrxmin2:jrxmax2,&
         jrxmin3:jrxmax3,1) - x(hrxmin1:hrxmax1,hrxmin2:hrxmax2,&
         hrxmin3:hrxmax3,1)))

      cent(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)  = (x(jrxmin1:jrxmax1,jrxmin2:jrxmax2,&
         jrxmin3:jrxmax3,1) + x(hrxmin1:hrxmax1,hrxmin2:hrxmax2,&
         hrxmin3:hrxmax3,1) - 2.0d0*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)) * vfield(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) / ((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1) - x(hrxmin1:hrxmax1,hrxmin2:hrxmax2,&
         hrxmin3:hrxmax3,1)) * (x(jrxmin1:jrxmax1,jrxmin2:jrxmax2,&
         jrxmin3:jrxmax3,1) - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)))
    
    case(2) ! Polar forward, backward, and central derivatives
      forw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)  = (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2) - x(htxmin1:htxmax1,htxmin2:htxmax2,&
         htxmin3:htxmax3,2)) * vfield(jtxmin1:jtxmax1,jtxmin2:jtxmax2,&
         jtxmin3:jtxmax3) / (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1) * (x(jtxmin1:jtxmax1,jtxmin2:jtxmax2,jtxmin3:jtxmax3,&
         2) - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)) * (x(jtxmin1:jtxmax1,jtxmin2:jtxmax2,jtxmin3:jtxmax3,&
         2) - x(htxmin1:htxmax1,htxmin2:htxmax2,htxmin3:htxmax3,2)))

      backw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = -(x(jtxmin1:jtxmax1,jtxmin2:jtxmax2,&
         jtxmin3:jtxmax3,2) - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)) * vfield(htxmin1:htxmax1,htxmin2:htxmax2,&
         htxmin3:htxmax3) / ( x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1) * (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2) - x(htxmin1:htxmax1,htxmin2:htxmax2,&
         htxmin3:htxmax3,2)) * (x(jtxmin1:jtxmax1,jtxmin2:jtxmax2,&
         jtxmin3:jtxmax3,2) - x(htxmin1:htxmax1,htxmin2:htxmax2,&
         htxmin3:htxmax3,2)))

      cent(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)  = (x(jtxmin1:jtxmax1,jtxmin2:jtxmax2,&
         jtxmin3:jtxmax3,2) + x(htxmin1:htxmax1,htxmin2:htxmax2,&
         htxmin3:htxmax3,2) - 2.0d0*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)) * vfield(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) / ( x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1) * (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2) - x(htxmin1:htxmax1,htxmin2:htxmax2,&
         htxmin3:htxmax3,2)) * (x(jtxmin1:jtxmax1,jtxmin2:jtxmax2,&
         jtxmin3:jtxmax3,2) - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)))
   
    
    case(3) ! Azimuthal forward, backward, and central derivatives
      forw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)  = (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3) - x(hpxmin1:hpxmax1,hpxmin2:hpxmax2,&
         hpxmin3:hpxmax3,3)) *  vfield(jpxmin1:jpxmax1,jpxmin2:jpxmax2,&
         jpxmin3:jpxmax3) / ( x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)*sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)) * (x(jpxmin1:jpxmax1,jpxmin2:jpxmax2,&
         jpxmin3:jpxmax3,3) - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)) * (x(jpxmin1:jpxmax1,jpxmin2:jpxmax2,&
         jpxmin3:jpxmax3,3) - x(hpxmin1:hpxmax1,hpxmin2:hpxmax2,&
         hpxmin3:hpxmax3,3)))

      backw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = -(x(jpxmin1:jpxmax1,jpxmin2:jpxmax2,&
         jpxmin3:jpxmax3,3) - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)) *  vfield(hpxmin1:hpxmax1,hpxmin2:hpxmax2,&
         hpxmin3:hpxmax3) / ( x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)*sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)) * (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3) - x(hpxmin1:hpxmax1,hpxmin2:hpxmax2,&
         hpxmin3:hpxmax3,3)) * (x(jpxmin1:jpxmax1,jpxmin2:jpxmax2,&
         jpxmin3:jpxmax3,3) - x(hpxmin1:hpxmax1,hpxmin2:hpxmax2,&
         hpxmin3:hpxmax3,3)))

      cent(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)  = (x(jpxmin1:jpxmax1,jpxmin2:jpxmax2,&
         jpxmin3:jpxmax3,3) + x(hpxmin1:hpxmax1,hpxmin2:hpxmax2,&
         hpxmin3:hpxmax3,3) - 2.0d0*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)) *  vfield(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) / ( x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)*sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2)) * (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3) - x(hpxmin1:hpxmax1,hpxmin2:hpxmax2,&
         hpxmin3:hpxmax3,3)) * (x(jpxmin1:jpxmax1,jpxmin2:jpxmax2,&
         jpxmin3:jpxmax3,3) - x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)))
   
    end select

    ! Total gradient for given velocity field
    grad_vn(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = backw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + cent(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + forw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

  end subroutine get_velocity_gradient

  !> Initialise (theta',phi') radiation angles coming from stellar disc
  subroutine rays_init(ntheta_point,nphi_point)
    use mod_global_parameters

    ! Subroutine arguments
    integer, intent(in) :: ntheta_point, nphi_point

    ! Local variables
    real(8) :: ymin, ymax, phipmin, phipmax, adum
    integer :: ii

    ! Minimum and maximum range of theta and phi rays
    ! NOTE: theta points are cast into y-space
    ymin    = 0.0d0
    ymax    = 1.0d0
    phipmin = -dpi !0.0d0
    phipmax = dpi !2.0d0*dpi
    ! dphi    = (phipmax - phipmin) / nphi_point

    if (mype == 0) then
      allocate(ay(ntheta_point))
      allocate(wy(ntheta_point))
      allocate(aphi(nphi_point))
      allocate(wphi(nphi_point))

      ! theta and phi ray positions and weights: Gauss-Legendre
      call gauss_legendre_quadrature(ymin,ymax,ntheta_point,ay,wy)
      call gauss_legendre_quadrature(phipmin,phipmax,nphi_point,aphi,wphi)

      ! theta rays and weights: uniform
      ! dth = 1.0d0 / nthetap
      ! adum = ymin + 0.5d0*dth
      ! do ip = 1,nthetap
      !   ay(ip) = adum
      !   wy(ip) = 1.0d0/nthetap
      !   adum = adum + dth
      !   !print*,'phipoints'
      !   !print*,ip,aphi(ip),wphi(ip),dphi
      ! enddo

      ! phi ray position and weights: uniform
      ! adum = phipmin + 0.5d0*dphi
      ! do ii = 1,nphi_point
      !   aphi(ii) = adum
      !   wphi(ii) = 1.0d0/nphi_point
      !   adum     = adum + dphi
      ! enddo

      print*, '==========================='
      print*, '    Radiation ray setup    '
      print*, '==========================='
      print*, 'Theta ray points + weights '
      do ii = 1,ntheta_point
        print*,ii,ay(ii),wy(ii)
      enddo
      print*
      print*, 'Phi ray points + weights   '
      do ii = 1,nphi_point
        print*,ii,aphi(ii),wphi(ii)
      enddo
      print*
    endif

    call MPI_BARRIER(icomm,ierrmpi)

    !===========================
    ! Broadcast what mype=0 read
    !===========================
    if (npe > 1) then
      call MPI_BCAST(ntheta_point,1,MPI_INTEGER,0,icomm,ierrmpi)
      call MPI_BCAST(nphi_point,1,MPI_INTEGER,0,icomm,ierrmpi)

      if (mype /= 0) then
        allocate(ay(ntheta_point))
        allocate(wy(ntheta_point))
        allocate(aphi(nphi_point))
        allocate(wphi(nphi_point))
      endif

      call MPI_BCAST(ay,ntheta_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(wy,ntheta_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(aphi,nphi_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(wphi,nphi_point,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

  end subroutine rays_init
  
  !> Fast Gauss-Legendre N-point quadrature algorithm by G. Rybicki
  subroutine gauss_legendre_quadrature(xlow,xhi,n,x,w)
    ! Given the lower and upper limits of integration xlow and xhi, and given n,
    ! this routine returns arrays x and w of length n, containing the abscissas
    ! and weights of the Gauss-Legendre N-point quadrature
    use mod_global_parameters

    ! Subroutine arguments
    real(8), intent(in)  :: xlow, xhi
    integer, intent(in)  :: n
    real(8), intent(out) :: x(n), w(n)

    ! Local variables
    integer :: i, j, m
    real(8) :: p1, p2, p3, pp, xl, xm, z, z1
    real(8), parameter :: error=3.0d-14

    m = (n + 1)/2
    xm = 0.5d0*(xhi + xlow)
    xl = 0.5d0*(xhi - xlow)

    do i = 1,m
      z = cos( dpi * (i - 0.25d0)/(n + 0.5d0) )
      z1 = 2.0d0 * z

      do while (abs(z1 - z) > error)
        p1 = 1.0d0
        p2 = 0.0d0

        do j = 1,n
          p3 = p2
          p2 = p1
          p1 = ( (2.0d0*j - 1.0d0)*z*p2 - (j - 1.0d0)*p3 )/j
        enddo

        pp = n*(z*p1 - p2) / (z*z - 1.0d0)
        z1 = z
        z = z1 - p1/pp
      enddo

      x(i)     = xm - xl*z
      x(n+1-i) = xm + xl*z
      w(i)     = 2.0d0*xl / ((1.0d0 - z*z) * pp*pp)
      w(n+1-i) = w(i)
    enddo

  end subroutine gauss_legendre_quadrature

end module mod_cak_force
