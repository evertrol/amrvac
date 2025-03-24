!> Magneto-hydrodynamics module
module mod_mhd_phys

#include "amrvac.h"

  use mod_global_parameters, only: std_len, const_c
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  use mod_physics
  use mod_comm_lib, only: mpistop
  use mod_functions_bfield, only: get_divb, mag

  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: mhd_energy = .true.

  !> Whether thermal conduction is used
  logical, public, protected              :: mhd_thermal_conduction = .false.
  !> type of fluid for thermal conduction
  type(tc_fluid), public, allocatable     :: tc_fl
  !> type of fluid for thermal emission synthesis
  type(te_fluid), public, allocatable     :: te_fl_mhd

  !> Whether radiative cooling is added
  logical, public, protected              :: mhd_radiative_cooling = .false.
  !> type of fluid for radiative cooling
  type(rc_fluid), public, allocatable     :: rc_fl

  !> Whether viscosity is added
  logical, public, protected              :: mhd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: mhd_gravity = .false.

  !> Whether rotating frame is activated
  logical, public, protected              :: mhd_rotating_frame = .false.

  !> Whether Hall-MHD is used
  logical, public, protected              :: mhd_Hall = .false.

  !> Whether Ambipolar term is used
  logical, public, protected              :: mhd_ambipolar = .false.

  !> Whether Ambipolar term is implemented using supertimestepping
  logical, public, protected              :: mhd_ambipolar_sts = .false.

  !> Whether Ambipolar term is implemented explicitly
  logical, public, protected              :: mhd_ambipolar_exp = .false.

  !> Whether particles module is added
  logical, public, protected              :: mhd_particles = .false.

  !> Whether magnetofriction is added
  logical, public, protected              :: mhd_magnetofriction = .false.

  !> Whether GLM-MHD is used to control div B
  logical, public, protected              :: mhd_glm = .false.

  !> Whether extended GLM-MHD is used with additional sources
  logical, public, protected              :: mhd_glm_extended = .true.

  !> Whether TRAC method is used
  logical, public, protected              :: mhd_trac = .false.

  !> Which TRAC method is used
  integer, public, protected              :: mhd_trac_type=1

  !> Height of the mask used in the TRAC method
  double precision, public, protected     :: mhd_trac_mask = 0.d0

  !> Distance between two adjacent traced magnetic field lines (in finest cell size)
  integer, public, protected              :: mhd_trac_finegrid=4

  !> Whether internal energy is solved instead of total energy
  logical, public, protected              :: mhd_internal_e = .false.

  !TODO this does not work with the splitting: check mhd_check_w_hde and mhd_handle_small_values_hde
  !> Whether hydrodynamic energy is solved instead of total energy
  logical, public, protected              :: mhd_hydrodynamic_e = .false.

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: mhd_glm_alpha = 0.5d0

  !TODO this does not work with the splitting: check mhd_check_w_semirelati and mhd_handle_small_values_semirelati
  !> Whether semirelativistic MHD equations (Gombosi 2002 JCP) are solved
  logical, public, protected              :: mhd_semirelativistic = .false.

  !> Whether boris simplified semirelativistic MHD equations (Gombosi 2002 JCP) are solved
  logical, public, protected              :: mhd_boris_simplification = &
     .false.

  !> Reduced speed of light for semirelativistic MHD: 2% of light speed
  double precision, public, protected     :: mhd_reduced_c = 0.02d0*const_c

  !> Whether plasma is partially ionized
  logical, public, protected              :: mhd_partial_ionization = .false.

  !> Whether CAK radiation line force is activated
  logical, public, protected              :: mhd_cak_force = .false.

  !> MHD fourth order
  logical, public, protected              :: mhd_4th_order = .false.

  !> whether split off equilibrium density
  logical, public :: has_equi_rho0 = .false.
  !> whether split off equilibrium thermal pressure
  logical, public :: has_equi_pe0 = .false.
  logical, public :: mhd_equi_thermal = .false.

  !> equi vars indices in the state%equi_vars array
  integer, public :: equi_rho0_ = -1
  integer, public :: equi_pe0_ = -1

  !> whether dump full variables (when splitting is used) in a separate dat file
  logical, public, protected              :: mhd_dump_full_vars = .false.

  !> Number of tracer species
  integer, public, protected              :: mhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_


  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> Indices of temperature
  integer, public, protected :: Te_

  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_
  integer, public, protected              :: Tweight_

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> The adiabatic index
  double precision, public                :: mhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: mhd_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public                :: mhd_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public                :: mhd_eta_hyper = 0.0d0

  !> TODO: what is this?
  double precision, public                :: mhd_etah = 0.0d0

  !> The MHD ambipolar coefficient
  double precision, public                :: mhd_eta_ambi = 0.0d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> The number of waves
  integer :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'linde'

  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'uct_contact'

  !> Whether divB is computed with a fourth order approximation
  logical, public, protected :: mhd_divb_4thorder = .false.

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'

  !> Use a compact way to add resistivity
  logical :: compactres   = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> clean initial divB
  logical, public :: clean_initial_divb     = .false.

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0
  !> Ionization fraction of H
  !> H_ion_fr = H+/(H+ + H)
  double precision, public, protected  :: H_ion_fr=1d0
  !> Ionization fraction of He
  !> He_ion_fr = (He2+ + He+)/(He2+ + He+ + He)
  double precision, public, protected  :: He_ion_fr=1d0
  !> Ratio of number He2+ / number He+ + He2+
  !> He_ion_fr2 = He2+/(He2+ + He+)
  double precision, public, protected  :: He_ion_fr2=1d0
  ! used for eq of state when it is not defined by units,
  ! the units do not contain terms related to ionization fraction
  ! and it is p = RR * rho * T
  double precision, public, protected  :: RR=1d0
  ! remove the below flag  and assume default value = .false.
  ! when eq state properly implemented everywhere
  ! and not anymore through units
  logical, public, protected :: eq_state_units = .true.

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*3)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*3)=0

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.

  !> Whether an total energy equation is used
  logical :: total_energy = .true.

  !> Whether an internal or hydrodynamic energy equation is used
  logical :: partial_energy = .false.

  !> Whether gravity work is included in energy equation
  logical :: gravity_energy

  !> gravity work is calculated use density times velocity or conservative momentum
  logical :: gravity_rhov = .false.

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

  !> inverse of squared speed of light c0 and reduced speed of light c
  double precision :: inv_squared_c0, inv_squared_c

  ! DivB cleaning methods
  integer, parameter :: divb_none          = 0
  integer, parameter :: divb_multigrid     = -1
  integer, parameter :: divb_glm           = 1
  integer, parameter :: divb_powel         = 2
  integer, parameter :: divb_janhunen      = 3
  integer, parameter :: divb_linde         = 4
  integer, parameter :: divb_lindejanhunen = 5
  integer, parameter :: divb_lindepowel    = 6
  integer, parameter :: divb_lindeglm      = 7
  integer, parameter :: divb_ct            = 8

  !define the subroutine interface for the ambipolar mask
  abstract interface

    subroutine mask_subroutine(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,res)
      use mod_global_parameters
      integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:nw)
      double precision, intent(inout) :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
    end subroutine mask_subroutine

    function fun_kin_en(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, inv_rho) result(ke)
      use mod_global_parameters, only: nw, ndim,block
      integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3, nw)
      double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end function fun_kin_en

  end interface

  procedure(mask_subroutine), pointer  :: usr_mask_ambipolar => null()
  procedure(sub_convert), pointer      :: mhd_to_primitive  => null()
  procedure(sub_convert), pointer      :: mhd_to_conserved  => null()
  procedure(sub_small_values), pointer :: mhd_handle_small_values => null()
  procedure(sub_get_pthermal), pointer :: mhd_get_pthermal  => null()
  procedure(sub_get_pthermal), pointer :: mhd_get_Rfactor   => null()
  procedure(sub_get_pthermal), pointer :: mhd_get_temperature=> null()
  procedure(sub_get_v), pointer        :: mhd_get_v         => null()
  procedure(fun_kin_en), pointer       :: mhd_kin_en        => null()
  ! Public methods
  public :: usr_mask_ambipolar
  public :: mhd_phys_init
  public :: mhd_kin_en
  public :: mhd_get_pthermal
  public :: mhd_get_temperature
  public :: mhd_get_v
  public :: mhd_get_rho
  public :: mhd_get_v_idim
  public :: mhd_to_conserved
  public :: mhd_to_primitive
  public :: mhd_get_csound2
  public :: mhd_e_to_ei
  public :: mhd_ei_to_e
  public :: mhd_face_to_center
  public :: get_divb
  public :: get_current
  !> needed  public if we want to use the ambipolar coefficient in the user file
  public :: multiplyAmbiCoef
  public :: get_normalized_divb
  public :: b_from_vector_potential
  public :: mhd_mag_en_all
  
  public :: mhd_clean_divb_multigrid
 

contains

  !> Read this module"s parameters from a file
  subroutine mhd_read_params(files)
    use mod_global_parameters
    use mod_particles, only: particles_eta, particles_etah
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_energy, mhd_n_tracer, mhd_gamma, mhd_adiab,mhd_eta,&
        mhd_eta_hyper, mhd_etah, mhd_eta_ambi, mhd_glm_alpha, mhd_glm_extended,&
        mhd_magnetofriction,mhd_thermal_conduction, mhd_radiative_cooling,&
        mhd_Hall, mhd_ambipolar, mhd_ambipolar_sts, mhd_gravity,&
       mhd_rotating_frame,mhd_viscosity, mhd_4th_order, typedivbfix,&
        source_split_divb, divbdiff,typedivbdiff, type_ct, compactres,&
        divbwave, He_abundance, H_ion_fr, He_ion_fr, He_ion_fr2,&
        eq_state_units, SI_unit, B0field ,mhd_dump_full_vars,B0field_forcefree,&
        Bdip, Bquad, Boct, Busr, mhd_particles, mhd_partial_ionization,&
       particles_eta, particles_etah,has_equi_rho0, has_equi_pe0,&
       mhd_equi_thermal,boundary_divbfix, boundary_divbfix_skip,&
        mhd_divb_4thorder, mhd_semirelativistic,mhd_boris_simplification,&
        mhd_reduced_c, clean_initial_divb, mhd_internal_e, mhd_hydrodynamic_e,&
        mhd_trac, mhd_trac_type, mhd_trac_mask, mhd_trac_finegrid,&
        mhd_cak_force

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mhd_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine mhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine mhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = mhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine mhd_write_info

  subroutine mhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init, particles_eta, particles_etah
    use mod_magnetofriction, only: magnetofriction_init
    use mod_rotating_frame, only: rotating_frame_init
    use mod_supertimestepping, only: sts_init, add_sts_method,&
       set_conversion_methods_to_head, set_error_handling_to_head
    use mod_cak_force, only: cak_init
    use mod_ionization_degree
    use mod_usr_methods, only: usr_Rfactor
    
    use mod_multigrid_coupling
   

    integer :: itr, idir

    call mhd_read_params(par_files)

    if(mhd_internal_e) then
      if(mhd_hydrodynamic_e) then
        mhd_hydrodynamic_e=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_hydrodynamic_e=F when mhd_internal_e=T'
      end if
    end if

    if(mhd_semirelativistic) then
      if(mhd_boris_simplification) then
        mhd_boris_simplification=.false.
        if(mype==0) write(*,*)&
'WARNING: set mhd_boris_simplification=F when mhd_semirelativistic=T'
      end if
      if(B0field) B0fieldAllocCoarse=.true.
    end if

    if(.not. mhd_energy) then
      if(mhd_internal_e) then
        mhd_internal_e=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_internal_e=F when mhd_energy=F'
      end if
      if(mhd_hydrodynamic_e) then
        mhd_hydrodynamic_e=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_hydrodynamic_e=F when mhd_energy=F'
      end if
      if(mhd_thermal_conduction) then
        mhd_thermal_conduction=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_thermal_conduction=F when mhd_energy=F'
      end if
      if(mhd_radiative_cooling) then
        mhd_radiative_cooling=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_radiative_cooling=F when mhd_energy=F'
      end if
      if(mhd_trac) then
        mhd_trac=.false.
        if(mype==0) write(*,*) 'WARNING: set mhd_trac=F when mhd_energy=F'
      end if
      if(mhd_partial_ionization) then
        mhd_partial_ionization=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_partial_ionization=F when mhd_energy=F'
      end if
    end if
    if(.not.eq_state_units) then
      if(mhd_partial_ionization) then
        mhd_partial_ionization=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_partial_ionization=F when eq_state_units=F'
      end if
    end if

    physics_type = "mhd"
    phys_energy=mhd_energy
    phys_internal_e=mhd_internal_e
    phys_trac=mhd_trac
    phys_trac_type=mhd_trac_type
    phys_partial_ionization=mhd_partial_ionization

    phys_gamma = mhd_gamma

    if(mhd_energy) then
      if(mhd_internal_e.or.mhd_hydrodynamic_e) then
        partial_energy=.true.
        total_energy=.false.
      else
        partial_energy=.false.
        total_energy=.true.
      end if
    else
      total_energy=.false.
    end if
    phys_total_energy=total_energy
    if(mhd_energy) then
      if(mhd_internal_e) then
        gravity_energy=.false.
      else
        gravity_energy=.true.
      end if
      if(has_equi_rho0.or.mhd_boris_simplification) then
        gravity_rhov=.true.
      end if
      if(mhd_semirelativistic.and..not.mhd_hydrodynamic_e) then
        gravity_rhov=.true.
      end if
    else
      gravity_energy=.false.
    end if

    
    if(mhd_trac .and. mhd_trac_type .le. 4) then
      mhd_trac_mask=bigdouble
      if(mype==0) write(*,*)&
          'WARNING: set mhd_trac_mask==bigdouble for global TRAC method'
    end if
    phys_trac_mask=mhd_trac_mask

    ! set default gamma for polytropic/isothermal process
    use_particles=mhd_particles
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    
    case ('multigrid')
       type_divb = divb_multigrid
       use_multigrid = .true.
       mg%operator_type = mg_laplacian
       phys_global_source_after => mhd_clean_divb_multigrid
   
    case ('glm')
      mhd_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_glm
    case ('powel', 'powell')
      type_divb = divb_powel
    case ('janhunen')
      type_divb = divb_janhunen
    case ('linde')
      type_divb = divb_linde
    case ('lindejanhunen')
      type_divb = divb_lindejanhunen
    case ('lindepowel')
      type_divb = divb_lindepowel
    case ('lindeglm')
      mhd_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_lindeglm
    case ('ct')
      type_divb = divb_ct
      stagger_grid = .true.
    case default
      call mpistop('Unknown divB fix')
    end select

    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=1
    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (mhd_energy) then
      nwwave = 8
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 7
      e_     = -1
      p_     = -1
    end if

    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)

    if (mhd_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    allocate(tracer(mhd_n_tracer))
    ! Set starting index of tracers
    do itr = 1, mhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    !  set temperature as an auxiliary variable to get ionization degree
    if(mhd_partial_ionization) then
      Te_ = var_set_auxvar('Te','Te')
    else
      Te_ = -1
    end if

    ! set cutoff temperature when using the TRAC method, as well as an auxiliary weight
    Tweight_ = -1
    if(mhd_trac) then
      Tcoff_ = var_set_wextra()
      iw_Tcoff=Tcoff_
      if(mhd_trac_type .ge. 3) then
        Tweight_ = var_set_wextra()
      endif
    else
      Tcoff_ = -1
    end if

    ! set indices of equi vars and update number_equi_vars
    number_equi_vars = 0
    if(has_equi_rho0) then
      number_equi_vars = number_equi_vars + 1
      equi_rho0_ = number_equi_vars
      iw_equi_rho = equi_rho0_
    endif
    if(has_equi_pe0) then
      number_equi_vars = number_equi_vars + 1
      equi_pe0_ = number_equi_vars
      iw_equi_p = equi_pe0_
      phys_equi_pe=.true.
    endif
    ! determine number of stagger variables
    nws=ndim

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nwflux))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nwflux])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    if(nwflux>mag(ndir)) then
      ! for flux of tracers, using hll flux
      flux_type(:,mag(ndir)+1:nwflux)=flux_hll
    end if

    if(ndim>1) then
      if(mhd_glm) then
        flux_type(:,psi_)=flux_special
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_special
        end do
      else
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_tvdlf
        end do
      end if
    end if

    phys_get_dt              => mhd_get_dt
    if(mhd_semirelativistic) then
      phys_get_cmax            => mhd_get_cmax_semirelati
    else
      phys_get_cmax            => mhd_get_cmax_origin
    end if
    phys_get_a2max           => mhd_get_a2max
    phys_get_tcutoff         => mhd_get_tcutoff
    phys_get_H_speed         => mhd_get_H_speed
    if(has_equi_rho0) then
      phys_get_cbounds         => mhd_get_cbounds_split_rho
    else if(mhd_semirelativistic) then
      phys_get_cbounds         => mhd_get_cbounds_semirelati
    else
      phys_get_cbounds         => mhd_get_cbounds
    end if
    if(mhd_hydrodynamic_e) then
      phys_to_primitive        => mhd_to_primitive_hde
      mhd_to_primitive         => mhd_to_primitive_hde
      phys_to_conserved        => mhd_to_conserved_hde
      mhd_to_conserved         => mhd_to_conserved_hde
    else if(mhd_semirelativistic) then
      phys_to_primitive        => mhd_to_primitive_semirelati
      mhd_to_primitive         => mhd_to_primitive_semirelati
      phys_to_conserved        => mhd_to_conserved_semirelati
      mhd_to_conserved         => mhd_to_conserved_semirelati
    else
      if(has_equi_rho0) then
        phys_to_primitive        => mhd_to_primitive_split_rho
        mhd_to_primitive         => mhd_to_primitive_split_rho
        phys_to_conserved        => mhd_to_conserved_split_rho
        mhd_to_conserved         => mhd_to_conserved_split_rho
      else if(mhd_internal_e) then
        phys_to_primitive        => mhd_to_primitive_inte
        mhd_to_primitive         => mhd_to_primitive_inte
        phys_to_conserved        => mhd_to_conserved_inte
        mhd_to_conserved         => mhd_to_conserved_inte
      else
        phys_to_primitive        => mhd_to_primitive_origin
        mhd_to_primitive         => mhd_to_primitive_origin
        phys_to_conserved        => mhd_to_conserved_origin
        mhd_to_conserved         => mhd_to_conserved_origin
      end if
    end if
    if(mhd_hydrodynamic_e) then
      phys_get_flux            => mhd_get_flux_hde
    else if(mhd_semirelativistic) then
      phys_get_flux            => mhd_get_flux_semirelati
    else
      if(B0field.or.has_equi_rho0.or.has_equi_pe0) then
        phys_get_flux            => mhd_get_flux_split
      else
        phys_get_flux            => mhd_get_flux
      end if
    end if
    if(mhd_boris_simplification) then
      phys_get_v                 => mhd_get_v_boris
      mhd_get_v                  => mhd_get_v_boris
      mhd_kin_en                 => mhd_kin_en_boris
    else
      phys_get_v                 => mhd_get_v_origin
      mhd_get_v                  => mhd_get_v_origin
      mhd_kin_en                 => mhd_kin_en_origin
    end if
    if(B0field.or.has_equi_rho0) then
      phys_add_source_geom     => mhd_add_source_geom_split
    else
      phys_add_source_geom     => mhd_add_source_geom
    end if
    phys_add_source          => mhd_add_source
    phys_check_params        => mhd_check_params
    phys_write_info          => mhd_write_info
 
    if(mhd_hydrodynamic_e) then
      phys_handle_small_values => mhd_handle_small_values_hde
      mhd_handle_small_values  => mhd_handle_small_values_hde
      phys_check_w             => mhd_check_w_hde
    else if(mhd_semirelativistic) then
      phys_handle_small_values => mhd_handle_small_values_semirelati
      mhd_handle_small_values  => mhd_handle_small_values_semirelati
      phys_check_w             => mhd_check_w_semirelati
    else
      if(mhd_internal_e) then
        phys_handle_small_values => mhd_handle_small_values_inte
        mhd_handle_small_values  => mhd_handle_small_values_inte
        phys_check_w             => mhd_check_w_inte
      else
        phys_handle_small_values => mhd_handle_small_values_origin
        mhd_handle_small_values  => mhd_handle_small_values_origin
        phys_check_w             => mhd_check_w_origin
      end if
    end if
 
    if(.not.mhd_energy) then
      phys_get_pthermal        => mhd_get_pthermal_iso
      mhd_get_pthermal         => mhd_get_pthermal_iso
    else
      if(mhd_internal_e) then
        phys_get_pthermal        => mhd_get_pthermal_eint
        mhd_get_pthermal         => mhd_get_pthermal_eint
      else if(mhd_hydrodynamic_e) then
        phys_get_pthermal        => mhd_get_pthermal_hde
        mhd_get_pthermal         => mhd_get_pthermal_hde
      else if(mhd_semirelativistic) then
        phys_get_pthermal        => mhd_get_pthermal_semirelati
        mhd_get_pthermal         => mhd_get_pthermal_semirelati
      else
        phys_get_pthermal        => mhd_get_pthermal_origin
        mhd_get_pthermal         => mhd_get_pthermal_origin
      end if
    end if
    if(number_equi_vars>0) then
      phys_set_equi_vars => set_equi_vars_grid
    endif

    if(type_divb==divb_glm) then
      phys_modify_wLR => mhd_modify_wLR
    end if

    ! choose Rfactor in ideal gas law
    if(mhd_partial_ionization) then
      mhd_get_Rfactor=>Rfactor_from_temperature_ionization
      phys_update_temperature => mhd_update_temperature
    else if(associated(usr_Rfactor)) then
      mhd_get_Rfactor=>usr_Rfactor
    else
      mhd_get_Rfactor=>Rfactor_from_constant_ionization
    end if

    if(mhd_partial_ionization) then
      mhd_get_temperature => mhd_get_temperature_from_Te
    else
      if(mhd_internal_e) then
        if(has_equi_pe0 .and. has_equi_rho0) then
          mhd_get_temperature => mhd_get_temperature_from_eint_with_equi
        else
          mhd_get_temperature => mhd_get_temperature_from_eint
        end if
      else
        if(has_equi_pe0 .and. has_equi_rho0) then
          mhd_get_temperature => mhd_get_temperature_from_etot_with_equi
        else
          mhd_get_temperature => mhd_get_temperature_from_etot
        end if
      end if
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_get_ct_velocity => mhd_get_ct_velocity
      phys_update_faces => mhd_update_faces
      phys_face_to_center => mhd_face_to_center
      phys_modify_wLR => mhd_modify_wLR
    else if(ndim>1) then
      phys_boundary_adjust => mhd_boundary_adjust
    end if

    
    ! clean initial divb
    if(clean_initial_divb) phys_clean_divb => mhd_clean_divb_multigrid
   

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call mhd_physical_units()

    if(.not. mhd_energy .and. mhd_thermal_conduction) then
      call mpistop("thermal conduction needs mhd_energy=T")
    end if
    if(.not. mhd_energy .and. mhd_radiative_cooling) then
      call mpistop("radiative cooling needs mhd_energy=T")
    end if

    ! resistive MHD needs diagonal ghost cells
    if(mhd_eta/=0.d0) phys_req_diagonal = .true.

    ! initialize thermal conduction module
    if (mhd_thermal_conduction) then
      phys_req_diagonal = .true.

      call sts_init()
      call tc_init_params(mhd_gamma)

      allocate(tc_fl)
      call tc_get_mhd_params(tc_fl,tc_params_read_mhd)
      call add_sts_method(mhd_get_tc_dt_mhd,mhd_sts_set_source_tc_mhd,e_,1,e_,&
         1,.false.)
      if(phys_internal_e) then
        if(has_equi_pe0 .and. has_equi_rho0) then
          tc_fl%get_temperature_from_conserved => &
             mhd_get_temperature_from_eint_with_equi
        else
          tc_fl%get_temperature_from_conserved => &
             mhd_get_temperature_from_eint
        end if
      else
        if(has_equi_pe0 .and. has_equi_rho0) then
          tc_fl%get_temperature_from_conserved => &
             mhd_get_temperature_from_etot_with_equi
        else
          tc_fl%get_temperature_from_conserved => &
             mhd_get_temperature_from_etot
        end if
      end if
      if(has_equi_pe0 .and. has_equi_rho0) then
        tc_fl%get_temperature_from_eint => &
           mhd_get_temperature_from_eint_with_equi
        if(mhd_equi_thermal) then
          tc_fl%has_equi = .true.
          tc_fl%get_temperature_equi => mhd_get_temperature_equi
          tc_fl%get_rho_equi => mhd_get_rho_equi
        else
          tc_fl%has_equi = .false.
        end if
      else
        tc_fl%get_temperature_from_eint => mhd_get_temperature_from_eint
      end if
      if(.not.mhd_internal_e) then
        if(mhd_hydrodynamic_e) then
          call set_conversion_methods_to_head(mhd_e_to_ei_hde,&
              mhd_ei_to_e_hde)
        else if(mhd_semirelativistic) then
          call set_conversion_methods_to_head(mhd_e_to_ei_semirelati,&
              mhd_ei_to_e_semirelati)
        else
          call set_conversion_methods_to_head(mhd_e_to_ei, mhd_ei_to_e)
        end if
      end if
      call set_error_handling_to_head(mhd_tc_handle_small_e)
      tc_fl%get_rho => mhd_get_rho
      tc_fl%e_ = e_
      tc_fl%Tcoff_ = Tcoff_
    end if

    ! Initialize radiative cooling module
    if (mhd_radiative_cooling) then
      call radiative_cooling_init_params(mhd_gamma,He_abundance)
      allocate(rc_fl)
      call radiative_cooling_init(rc_fl,rc_params_read)
      rc_fl%get_rho => mhd_get_rho
      rc_fl%get_pthermal => mhd_get_pthermal
      rc_fl%get_var_Rfactor => mhd_get_Rfactor
      rc_fl%e_ = e_
      rc_fl%Tcoff_ = Tcoff_
      if(has_equi_pe0 .and. has_equi_rho0 .and. mhd_equi_thermal) then
        rc_fl%has_equi = .true.
        rc_fl%get_rho_equi => mhd_get_rho_equi
        rc_fl%get_pthermal_equi => mhd_get_pe_equi
      else
        rc_fl%has_equi = .false.
      end if
    end if
    allocate(te_fl_mhd)
    te_fl_mhd%get_rho=> mhd_get_rho
    te_fl_mhd%get_pthermal=> mhd_get_pthermal
    te_fl_mhd%get_var_Rfactor => mhd_get_Rfactor

    phys_te_images => mhd_te_images

    ! Initialize viscosity module
    if (mhd_viscosity) call viscosity_init(phys_wider_stencil,&
       phys_req_diagonal)

    ! Initialize gravity module
    if(mhd_gravity) then
      call gravity_init()
    end if

    ! Initialize rotating frame module
    if(mhd_rotating_frame) call rotating_frame_init()

    ! Initialize particles module
    if(mhd_particles) then
      call particles_init()
      if (particles_eta  < zero) particles_eta = mhd_eta
      if (particles_etah < zero) particles_eta = mhd_etah
      phys_req_diagonal = .true.
      if(mype==0) then
         write(*,*) '*****Using particles:        with mhd_eta, mhd_etah :',&
             mhd_eta, mhd_etah
         write(*,*) '*****Using particles: particles_eta, particles_etah :',&
             particles_eta, particles_etah
      end if
    end if

    ! initialize magnetofriction module
    if(mhd_magnetofriction) then
      phys_req_diagonal = .true.
      call magnetofriction_init()
    end if

    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in mhd_get_flux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if(mhd_hall) then
      phys_req_diagonal = .true.
      if(mhd_4th_order) then
        phys_wider_stencil = 2
      else
        phys_wider_stencil = 1
      end if
    end if

    if(mhd_ambipolar) then
      phys_req_diagonal = .true.
      if(mhd_ambipolar_sts) then
        call sts_init()
        if(mhd_internal_e) then
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,mag(1),&
             ndir,mag(1),ndir,.true.)
        else
          call add_sts_method(get_ambipolar_dt,sts_set_source_ambipolar,&
             mom(ndir)+1,mag(ndir)-mom(ndir),mag(1),ndir,.true.)
        end if
      else
        mhd_ambipolar_exp=.true.
        ! For flux ambipolar term, we need one more reconstructed layer since currents are computed
        ! in mhd_get_flux: assuming one additional ghost layer (two for FOURTHORDER) was
        ! added in nghostcells.
        if(mhd_4th_order) then
          phys_wider_stencil = 2
        else
          phys_wider_stencil = 1
        end if
      end if
    end if

    ! initialize ionization degree table
    if(mhd_partial_ionization) call ionization_degree_init()

    ! Initialize CAK radiation force module
    if (mhd_cak_force) call cak_init(mhd_gamma)

  end subroutine mhd_phys_init


  subroutine mhd_te_images
    use mod_global_parameters
    use mod_thermal_emission

    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_mhd)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_mhd)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_mhd)
      case('WIvtiCCmpi','WIvtuCCmpi')
        call get_whitelight_image(unitconvert,te_fl_mhd)
      case default
        call mpistop("Error in synthesize emission: Unknown convert_type")
      end select
  end subroutine mhd_te_images


!!start th cond
  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  mhd_sts_set_source_tc_mhd(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
     fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_mhd
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, igrid, nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_mhd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
       fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine mhd_sts_set_source_tc_mhd

  function mhd_get_tc_dt_mhd(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,&
     x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_mhd
 
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_mhd(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,x,tc_fl)
  end function mhd_get_tc_dt_mhd

  subroutine mhd_tc_handle_small_e(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, step)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    integer, intent(in)    :: step
    character(len=140) :: error_msg

    write(error_msg,"(a,i3)") "Thermal conduction step ", step
    call mhd_handle_small_ei(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,e_,error_msg)
  end subroutine mhd_tc_handle_small_e

  ! fill in tc_fluid fields from namelist
  subroutine tc_params_read_mhd(fl)
    use mod_global_parameters, only: unitpar,par_files
    type(tc_fluid), intent(inout) :: fl

    integer                      :: n

    ! list parameters
    logical :: tc_perpendicular=.false.
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0
    double precision :: tc_k_perp=0d0
    character(len=std_len)  :: tc_slope_limiter="MC"

    namelist /tc_list/ tc_perpendicular, tc_saturate, tc_slope_limiter,&
        tc_k_para, tc_k_perp

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, tc_list, end=111)
111     close(unitpar)
    end do

    fl%tc_perpendicular = tc_perpendicular
    fl%tc_saturate = tc_saturate
    fl%tc_k_para = tc_k_para
    fl%tc_k_perp = tc_k_perp
    select case(tc_slope_limiter)
     case ('no','none')
       fl%tc_slope_limiter = 0
     case ('MC')
       ! montonized central limiter Woodward and Collela limiter (eq.3.51h), a factor of 2 is pulled out
       fl%tc_slope_limiter = 1
     case('minmod')
       ! minmod limiter
       fl%tc_slope_limiter = 2
     case ('superbee')
       ! Roes superbee limiter (eq.3.51i)
       fl%tc_slope_limiter = 3
     case ('koren')
       ! Barry Koren Right variant
       fl%tc_slope_limiter = 4
     case default
       call mpistop("Unknown tc_slope_limiter, choose MC, minmod")
    end select
  end subroutine tc_params_read_mhd
!!end th cond

!!rad cool
  subroutine rc_params_read(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_constants, only: bigdouble
    type(rc_fluid), intent(inout) :: fl
    integer                      :: n
    ! list parameters
    integer :: ncool = 4000
    double precision :: cfrac=0.1d0
  
    !> Name of cooling curve
    character(len=std_len)  :: coolcurve='JCcorona'
  
    !> Name of cooling method
    character(len=std_len)  :: coolmethod='exact'
  
    !> Fixed temperature not lower than tlow
    logical    :: Tfix=.false.
  
    !> Lower limit of temperature
    double precision   :: tlow=bigdouble
  
    !> Add cooling source in a split way (.true.) or un-split way (.false.)
    logical    :: rc_split=.false.


    namelist /rc_list/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix,&
        rc_split

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, rc_list, end=111)
111     close(unitpar)
    end do

    fl%ncool=ncool
    fl%coolcurve=coolcurve
    fl%coolmethod=coolmethod
    fl%tlow=tlow
    fl%Tfix=Tfix
    fl%rc_split=rc_split
    fl%cfrac=cfrac

  end subroutine rc_params_read
!! end rad cool

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid_faces(igrid,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in) :: igrid, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    double precision :: delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: xC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim),xshift1,xshift2,xshift3
    integer :: idims, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ix, idims2

    if(slab_uniform)then
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1)=rnode(rpdx1_,&
         igrid)
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2)=rnode(rpdx2_,&
         igrid)
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3)=rnode(rpdx3_,&
         igrid)
    else
      ! for all non-cartesian and stretched cartesian coordinates
      delx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)=ps(igrid)%dx(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)
    endif

    do idims=1,ndim
      hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
      hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
      hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
      if(stagger_grid) then
        ! ct needs all transverse cells
        ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
        ixCmax2=ixOmax2+nghostcells-nghostcells*kr(idims,2)
        ixCmax3=ixOmax3+nghostcells-nghostcells*kr(idims,3)
        ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1)
        ixCmin2=hxOmin2-nghostcells+nghostcells*kr(idims,2)
        ixCmin3=hxOmin3-nghostcells+nghostcells*kr(idims,3);
      else
        ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
        ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=hxOmin1
        ixCmin2=hxOmin2;ixCmin3=hxOmin3;
      end if
      ! always xshift=0 or 1/2
      xshift1=half*(one-kr(1,idims));xshift2=half*(one-kr(2,idims))
      xshift3=half*(one-kr(3,idims));
      do idims2=1,ndim
        select case(idims2)
        case(1)
          do ix = ixCmin1,ixCmax1
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ix,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1)=x(ix,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,1)+(half-xshift1)*delx(ix,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,1)
          end do
        case(2)
          do ix = ixCmin2,ixCmax2
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ixCmin1:ixCmax1,ix,ixCmin3:ixCmax3,2)=x(ixCmin1:ixCmax1,ix,&
               ixCmin3:ixCmax3,2)+(half-xshift2)*delx(ixCmin1:ixCmax1,ix,&
               ixCmin3:ixCmax3,2)
          end do
        case(3)
          do ix = ixCmin3,ixCmax3
            ! xshift=half: this is the cell center coordinate
            ! xshift=0: this is the cell edge i+1/2 coordinate
            xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ix,3)=x(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ix,3)+(half-xshift3)*delx(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ix,3)
          end do
        end select
      end do
      call usr_set_equi_vars(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,xC,&
         ps(igrid)%equi_vars(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:number_equi_vars,idims))
    end do

  end subroutine set_equi_vars_grid_faces

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid(igrid)
    use mod_global_parameters
    use mod_usr_methods
  
    integer, intent(in) :: igrid

    !values at the center
    call usr_set_equi_vars(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixGlo1,&
       ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ps(igrid)%x,&
       ps(igrid)%equi_vars(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
       1:number_equi_vars,0))

    !values at the interfaces
    call set_equi_vars_grid_faces(igrid,ps(igrid)%x,ixGlo1,ixGlo2,ixGlo3,&
       ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3)
 
  end subroutine set_equi_vars_grid

  ! w, wnew conserved, add splitted variables back to wnew
  function convert_vars_splitting(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
      nwc) result(wnew)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nwc
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision   :: wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        1:nwc)
    double precision   :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call  mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3))
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_) = rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(:)) =  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:))

    if (B0field) then
      ! add background magnetic field B0 to B
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0)
    else
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))
    end if

    if(mhd_energy) then
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
      if(has_equi_pe0) then
        wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) = wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe0_,0)* inv_gamma_1
      end if
      if(B0field .and. .not. mhd_internal_e) then
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)+0.5d0*sum(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,:,0)**2,dim=ndim+1) + sum(w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))*block%B0(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0),dim=ndim+1)
      end if
    end if

  end function convert_vars_splitting

  subroutine mhd_check_params
    use mod_global_parameters
    use mod_usr_methods
    use mod_convert, only: add_convert_method

    ! after user parameter setting
    gamma_1=mhd_gamma-1.d0
    if (.not. mhd_energy) then
       if (mhd_gamma <= 0.0d0) call mpistop ("Error: mhd_gamma <= 0")
       if (mhd_adiab < 0.0d0) call mpistop ("Error: mhd_adiab < 0")
       small_pressure = mhd_adiab*small_density**mhd_gamma
    else
       if (mhd_gamma <= 0.0d0 .or. mhd_gamma == 1.0d0) call mpistop &
          ("Error: mhd_gamma <= 0 or mhd_gamma == 1")
       inv_gamma_1=1.d0/gamma_1
       small_e = small_pressure * inv_gamma_1
    end if

    if (number_equi_vars > 0 .and. .not. associated(usr_set_equi_vars)) then
      call mpistop("usr_set_equi_vars has to be implemented in the user file")
    endif
    if(convert .or. autoconvert) then
      if(convert_type .eq. 'dat_generic_mpi') then
        if(mhd_dump_full_vars) then
          if(mype .eq. 0) print*, " add conversion method: split -> full "
          call add_convert_method(convert_vars_splitting, nw, cons_wnames,&
              "new")
        endif
      endif
    endif
  end subroutine mhd_check_params

  subroutine mhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0,c_lightspeed
    double precision :: a,b
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
      c_lightspeed=c_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi ! G2,G2,G2 cm2,cm2,cm2 dyne-1,dyne-1,dyne-1
      c_lightspeed=const_c
    end if
    if(eq_state_units) then
      a = 1d0 + 4d0 * He_abundance
      if(mhd_partial_ionization) then
        b = 2+.3d0
      else
        b = 1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + 1d0)+1d0)
      end if
      RR = 1d0
    else
      a = 1d0
      b = 1d0
      RR = (1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + &
         1d0)+1d0))/(1d0 + 4d0 * He_abundance)
    end if
    if(unit_density/=1.d0) then
      unit_numberdensity=unit_density/(a*mp)
    else
      ! unit of numberdensity is independent by default
      unit_density=a*mp*unit_numberdensity
    end if
    if(unit_velocity/=1.d0) then
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    else if(unit_pressure/=1.d0) then
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    else if(unit_magneticfield/=1.d0) then
      unit_pressure=unit_magneticfield**2/miu0
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
    else if(unit_temperature/=1.d0) then
      unit_pressure=b*unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
    end if
    if(unit_time/=1.d0) then
      unit_length=unit_time*unit_velocity
    else
      ! unit of length is independent by default
      unit_time=unit_length/unit_velocity
    end if
    ! Additional units needed for the particles
    c_norm=c_lightspeed/unit_velocity
    unit_charge=unit_magneticfield*unit_length**2/unit_velocity/miu0
    if (.not. SI_unit) unit_charge = unit_charge*const_c
    unit_mass=unit_density*unit_length**3

    if(mhd_semirelativistic.or.mhd_boris_simplification) then
      if(mhd_reduced_c<1.d0) then
        ! dimensionless speed
        inv_squared_c0=1.d0
        inv_squared_c=1.d0/mhd_reduced_c**2
      else
        inv_squared_c0=(unit_velocity/c_lightspeed)**2
        inv_squared_c=(unit_velocity/mhd_reduced_c)**2
      end if
    end if

  end subroutine mhd_physical_units

  subroutine mhd_check_w_semirelati(primitive,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision :: tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir),Ba(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)
    double precision :: v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndir),gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    logical, intent(inout) :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    integer :: idir, jdir, kdir

    flag=.false.
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_) = .true.

    if(mhd_energy) then
      if(primitive) then
        where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           p_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) = .true.
      else
        if(mhd_internal_e) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_)
        else
          if(B0field) then
            Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(1:ndir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,1:ndir,b0i)
          else
            Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(1:ndir))
          end if
          inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,rho_)
          b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=sum(Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,:)**2,dim=ndim+1)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)>smalldouble)
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=1.d0/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
          else where
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=0.d0
          end where
          do idir=1,ndir
            b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               idir)=Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               idir)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=sum(b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,:)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mom(:)),dim=ndim+1)
          ! Va^2/c^2
          b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*inv_squared_c
          ! equation (15)
          gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=1.d0/(1.d0+b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3))
          ! Convert momentum to velocity
          do idir = 1, ndir
             v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                idir) = gamma2*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mom(idir))+b2*b(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)*tmp)*inv_rho
          end do
          ! E=Bxv
          b=0.d0
          do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
            if(lvc(idir,jdir,kdir)==1)then
              b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)=b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)+Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 jdir)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
            else if(lvc(idir,jdir,kdir)==-1)then
              b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)=b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idir)-Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 jdir)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
            end if
          end do; end do; end do
          ! Calculate internal e = e-eK-eB-eE
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_)-half*(sum(v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,:)**2,dim=ndim+1)*w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)+sum(w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))**2,&
             dim=ndim+1)+sum(b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             :)**2,dim=ndim+1)*inv_squared_c)
        end if
        where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) = .true.
      end if
    end if

  end subroutine mhd_check_w_semirelati

  subroutine mhd_check_w_origin(primitive,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    logical, intent(inout) :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    flag=.false.
    if(has_equi_rho0) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,0)
      where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) < small_density) flag(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) = .true.
    else
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_) = .true.
    end if

    if(mhd_energy) then
      if(primitive) then
        if(has_equi_pe0) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_)+block%equi_vars(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_pressure) flag(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_) = .true.
        else
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_) = .true.
        end if
      else
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-mhd_kin_en(w,ixImin1,ixImin2,&
           ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
           ixOmax2,ixOmax3)-mhd_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
        if(has_equi_pe0) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe0_,0)*inv_gamma_1
        end if
        where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) = .true.
      end if
    end if

  end subroutine mhd_check_w_origin

  subroutine mhd_check_w_inte(primitive,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    logical, intent(inout) :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    flag=.false.
    if(has_equi_rho0) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,0)
      where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) < small_density) flag(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) = .true.
    else
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_) = .true.
    end if

    if(mhd_energy) then
      if(primitive) then
        if(has_equi_pe0) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_)+block%equi_vars(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_pressure) flag(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_) = .true.
        else
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_) = .true.
        end if
      else
        if(has_equi_pe0) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_)+block%equi_vars(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)*inv_gamma_1
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_) = .true.
        else
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_) = .true.
        end if
      end if
    end if

  end subroutine mhd_check_w_inte

  subroutine mhd_check_w_hde(primitive,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,flag)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    logical, intent(inout) :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    flag=.false.
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_) = .true.

    if(mhd_energy) then
      if(primitive) then
        where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) = .true.
      else
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-mhd_kin_en(w,ixImin1,ixImin2,&
           ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
           ixOmax2,ixOmax3)
        where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) = .true.
      end if
    end if

  end subroutine mhd_check_w_hde

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_origin(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision :: inv_gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    !if (fix_small_values) then
    !  call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_conserved')
    !end if

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1)*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)+mhd_mag_en(w, ixImin1,ixImin2,&
         ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
         ixOmax2,ixOmax3)
    end if

    if(mhd_boris_simplification) then
      ! rho*(1+B^2/rho/c^2)
      inv_gamma2=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = inv_gamma2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(idir))
      end do
    else
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idir))
      end do
    end if
  end subroutine mhd_to_conserved_origin

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_hde(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision :: inv_gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1)*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    end if

    ! Convert velocity to momentum
    if(mhd_boris_simplification) then
      ! rho*(1+B^2/rho/c^2)
      inv_gamma2=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = inv_gamma2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(idir))
      end do
    else
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idir))
      end do
    end if
  end subroutine mhd_to_conserved_hde

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_inte(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision :: inv_gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)*inv_gamma_1
    end if

    if(mhd_boris_simplification) then
      ! rho*(1+B^2/rho/c^2)
      inv_gamma2=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = inv_gamma2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(idir))
      end do
    else
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idir))
      end do
    end if
  end subroutine mhd_to_conserved_inte

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_split_rho(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision                :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), inv_gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    !if (fix_small_values) then
    !  call mhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'mhd_to_conserved')
    !end if

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,b0i)
    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      if(mhd_internal_e) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           p_)*inv_gamma_1
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1)*rho(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)+mhd_mag_en(w, ixImin1,ixImin2,&
           ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
           ixOmax2,ixOmax3)
      end if
    end if

    if(mhd_boris_simplification) then
      ! rho*(1+B^2/rho/c^2)
      inv_gamma2=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)*inv_squared_c
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idir)) = inv_gamma2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(idir))
      end do
    else
      ! Convert velocity to momentum
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, mom(idir))
      end do
    end if
  end subroutine mhd_to_conserved_split_rho

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved_semirelati(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision :: E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndir), B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir),&
        S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)
    integer                         :: idir, jdir, kdir

    if(B0field) then
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,b0i)
    else
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))
    end if
    E=0.d0
    do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
      if(lvc(idir,jdir,kdir)==1)then
        E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)+B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(kdir))
      else if(lvc(idir,jdir,kdir)==-1)then
        E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)-B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(kdir))
      end if
    end do; end do; end do

    if(mhd_internal_e) then
      ! internal energy
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)*inv_gamma_1
    else if(mhd_energy) then
      ! equation (9)
      ! Calculate total energy from internal, kinetic and magnetic energy
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)*inv_gamma_1+half*(sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1)*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)+sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))**2,&
         dim=ndim+1)+sum(E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)**2,dim=ndim+1)*inv_squared_c)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          rho_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))
    end do
    ! equation (5) Poynting vector
    S=0.d0
    do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
      if(lvc(idir,jdir,kdir)==1)then
        S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)+E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           jdir)*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
      else if(lvc(idir,jdir,kdir)==-1)then
        S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)-E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           jdir)*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
      end if
    end do; end do; end do
    ! equation (9)
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))+S(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir)*inv_squared_c
    end do
  end subroutine mhd_to_conserved_semirelati

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_origin(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call mhd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3, 'mhd_to_primitive_origin')
    end if

    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)

    if(mhd_boris_simplification) then
      gamma2=inv_rho/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*inv_rho
      end do
    end if

    ! Calculate pressure = (gamma-1) * (e-ek-eb)
    if(mhd_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)-0.5d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(:))**2,dim=ndim+1)-mhd_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3))
    end if

  end subroutine mhd_to_primitive_origin

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_hde(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call mhd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3, 'mhd_to_primitive_hde')
    end if

    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)

    if(mhd_boris_simplification) then
      gamma2=inv_rho/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*inv_rho
      end do
    end if

    ! Calculate pressure = (gamma-1) * (e-ek)
    if(mhd_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)-0.5d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(:))**2,dim=ndim+1))
    end if

  end subroutine mhd_to_primitive_hde

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_inte(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call mhd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3, 'mhd_to_primitive_inte')
    end if

    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)

    ! Calculate pressure = (gamma-1) * e_internal
    if(mhd_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)*gamma_1
    end if

    if(mhd_boris_simplification) then
      gamma2=inv_rho/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*inv_rho
      end do
    end if

  end subroutine mhd_to_primitive_inte

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_split_rho(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                         :: idir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call mhd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3, 'mhd_to_primitive_split_rho')
    end if

    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,equi_rho0_,b0i))

    if(mhd_boris_simplification) then
      gamma2=inv_rho/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*gamma2
      end do
    else
      ! Convert momentum to velocity
      do idir = 1, ndir
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*inv_rho
      end do
    end if

    ! Calculate pressure = (gamma-1) * (e-ek-eb)
    if(mhd_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)-0.5d0*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         equi_rho0_,b0i))*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1)-mhd_mag_en(w,ixImin1,ixImin2,&
         ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
         ixOmax2,ixOmax3))
    end if

  end subroutine mhd_to_primitive_split_rho

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive_semirelati(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    double precision                :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndir), Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir),&
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer                         :: idir, jdir, kdir

    if (fix_small_values) then
      ! fix small values preventing NaN numbers in the following converting
      call mhd_handle_small_values_semirelati(.false., w, x, ixImin1,ixImin2,&
         ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
         ixOmax2,ixOmax3, 'mhd_to_primitive_semirelati')
    end if

    if(B0field) then
      Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,b0i)
    else
      Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))
    end if

    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)

    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sum(Ba(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,dim=ndim+1)
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sqrt(b2(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)>smalldouble)
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else where
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=0.d0
    end where
    do idir=1,ndir
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)=Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=sum(b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mom(:)),dim=ndim+1)

    ! Va^2/c^2
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=b2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*inv_squared_c
    ! equation (15)
    gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.d0/(1.d0+b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir)) = gamma2*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, mom(idir))+b2*b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,idir)*tmp)*inv_rho
    end do

    if(mhd_internal_e) then
      ! internal energy to pressure
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)=gamma_1*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
    else if(mhd_energy) then
      ! E=Bxv
      b=0.d0
      do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
        if(lvc(idir,jdir,kdir)==1)then
          b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)+Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(kdir))
        else if(lvc(idir,jdir,kdir)==-1)then
          b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)-Ba(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(kdir))
        end if
      end do; end do; end do
      ! Calculate pressure = (gamma-1) * (e-eK-eB-eE)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)-half*(sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(:))**2,dim=ndim+1)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)+sum(b(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,dim=ndim+1)*inv_squared_c))
    end if

  end subroutine mhd_to_primitive_semirelati

  !> Transform internal energy to total energy
  subroutine mhd_ei_to_e(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    ! Calculate total energy from internal, kinetic and magnetic energy
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)+mhd_kin_en(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3)+mhd_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3)

  end subroutine mhd_ei_to_e

  !> Transform internal energy to hydrodynamic energy
  subroutine mhd_ei_to_e_hde(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    ! Calculate hydrodynamic energy from internal and kinetic
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)+mhd_kin_en(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3)

  end subroutine mhd_ei_to_e_hde

  !> Transform internal energy to total energy and velocity to momentum
  subroutine mhd_ei_to_e_semirelati(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,p_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)*gamma_1
    call mhd_to_conserved_semirelati(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x)

  end subroutine mhd_ei_to_e_semirelati

  !> Transform total energy to internal energy
  subroutine mhd_e_to_ei(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    ! Calculate ei = e - ek - eb
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)-mhd_kin_en(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3)-mhd_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3)

    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,e_,&
         'mhd_e_to_ei')
    end if

  end subroutine mhd_e_to_ei

  !> Transform hydrodynamic energy to internal energy
  subroutine mhd_e_to_ei_hde(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    ! Calculate ei = e - ek
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)-mhd_kin_en(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3)

    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,e_,&
         'mhd_e_to_ei_hde')
    end if

  end subroutine mhd_e_to_ei_hde

  !> Transform total energy to internal energy and momentum to velocity
  subroutine mhd_e_to_ei_semirelati(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    call mhd_to_primitive_semirelati(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x)
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,p_)*inv_gamma_1

  end subroutine mhd_e_to_ei_semirelati

  subroutine mhd_handle_small_values_semirelati(primitive, w, x, ixImin1,&
     ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
     ixOmax2,ixOmax3, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: pressure(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), inv_rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), b2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        gamma2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir), v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir),&
        Ba(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)
    integer :: idir, jdir, kdir, ix1,ix2,ix3
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    flag=.false.
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_) = .true.

    if(mhd_energy) then
      if(primitive) then
        where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           p_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) = .true.
      else
        if(mhd_internal_e) then
          pressure(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)=gamma_1*w(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,e_)
          where(pressure(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_pressure) flag(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_) = .true.
        else
          if(B0field) then
            Ba(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               1:ndir)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               mag(1:ndir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3,1:ndir,b0i)
          else
            Ba(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               1:ndir)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               mag(1:ndir))
          end if
          inv_rho(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3) = 1d0/w(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,rho_)
          b2(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)=sum(Ba(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,:)**2,dim=ndim+1)
          tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)=sqrt(b2(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3))
          where(tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)>smalldouble)
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3)=1.d0/tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3)
          else where
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=0.d0
          end where
          do idir=1,ndir
            b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               idir)=Ba(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               idir)*tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
          end do
          tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)=sum(b(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,:)*w(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,mom(:)),dim=ndim+1)
          ! Va^2/c^2
          b2(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)=b2(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)*inv_rho(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)*inv_squared_c
          ! equation (15)
          gamma2(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)=1.d0/(1.d0+b2(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3))
          ! Convert momentum to velocity
          do idir = 1, ndir
             v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                idir) = gamma2*(w(ixImin1:ixImax1,ixImin2:ixImax2,&
                ixImin3:ixImax3, mom(idir))+b2*b(ixImin1:ixImax1,&
                ixImin2:ixImax2,ixImin3:ixImax3,idir)*tmp(ixImin1:ixImax1,&
                ixImin2:ixImax2,ixImin3:ixImax3))*inv_rho(ixImin1:ixImax1,&
                ixImin2:ixImax2,ixImin3:ixImax3)
          end do
          ! E=Bxv
          b=0.d0
          do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
            if(lvc(idir,jdir,kdir)==1)then
              b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 idir)=b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 idir)+Ba(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 jdir)*v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,kdir)
            else if(lvc(idir,jdir,kdir)==-1)then
              b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 idir)=b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 idir)-Ba(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 jdir)*v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,kdir)
            end if
          end do; end do; end do
          ! Calculate pressure p = (gamma-1)(e-eK-eB-eE)
          pressure(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3)=gamma_1*(w(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,e_)-half*(sum(v(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,:)**2,dim=ndim+1)*w(ixImin1:ixImax1,&
             ixImin2:ixImax2,ixImin3:ixImax3,rho_)+sum(w(ixImin1:ixImax1,&
             ixImin2:ixImax2,ixImin3:ixImax3,mag(:))**2,&
             dim=ndim+1)+sum(b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
             :)**2,dim=ndim+1)*inv_squared_c))
          where(pressure(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_pressure) flag(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_) = .true.
        end if
      end if
    end if

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_) = small_density
        if(mhd_energy) then
          if(primitive) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               e_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               p_) = small_pressure
          else
            if(mhd_internal_e) then
              ! internal energy
              do ix3=ixOmin3,ixOmax3
              do ix2=ixOmin2,ixOmax2
              do ix1=ixOmin1,ixOmax1
                if(flag(ix1,ix2,ix3,e_)) then
                  w(ix1,ix2,ix3,e_)=small_pressure*inv_gamma_1
                end if
              end do
              end do
              end do
            else
              do ix3=ixOmin3,ixOmax3
              do ix2=ixOmin2,ixOmax2
              do ix1=ixOmin1,ixOmax1
                if(flag(ix1,ix2,ix3,e_)) then
                  w(ix1,ix2,ix3,e_)=small_pressure*inv_gamma_1+half*(sum(v(ix1,&
                     ix2,ix3,:)**2)*w(ix1,ix2,ix3,rho_)+sum(w(ix1,ix2,ix3,&
                     mag(:))**2)+sum(b(ix1,ix2,ix3,:)**2)*inv_squared_c)
                end if
              end do
              end do
              end do
            end if
          end if
        end if
      case ("average")
        ! do averaging of density
        call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
            flag, rho_)
        if(mhd_energy) then
          if(primitive) then
            call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
                flag, p_)
          else
            if(mhd_internal_e) then
              ! internal energy
              w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 e_)=pressure(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
              call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
                 ixOmax3, w, x, flag, p_)
              w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 e_)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 p_)*inv_gamma_1
            else
              w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 e_)=pressure(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
              call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
                 ixOmax3, w, x, flag, p_)
              w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 e_)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                 p_)*inv_gamma_1+half*(sum(v(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3,:)**2,dim=ndim+1)*w(ixImin1:ixImax1,&
                 ixImin2:ixImax2,ixImin3:ixImax3,rho_)+sum(w(ixImin1:ixImax1,&
                 ixImin2:ixImax2,ixImin3:ixImax3,mag(:))**2,&
                 dim=ndim+1)+sum(b(ixImin1:ixImax1,ixImin2:ixImax2,&
                 ixImin3:ixImax3,:)**2,dim=ndim+1)*inv_squared_c)
            end if
          end if
        end if
      case default
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if
  end subroutine mhd_handle_small_values_semirelati

  subroutine mhd_handle_small_values_origin(primitive, w, x, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision :: tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call phys_check_w(primitive, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_rho0) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_) = small_density-block%equi_vars(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,0)
        else
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_) = small_density
        end if
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom(idir)) = 0.0d0
          end if
        end do
        if(mhd_energy) then
          if(primitive) then
            if(has_equi_pe0) then
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = small_pressure - &
                 block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,equi_pe0_,0)
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 p_) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            else
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 p_) = small_pressure
            end if
          else
            ! conserved
            if(has_equi_pe0) then
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = small_e - block%equi_vars(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)*inv_gamma_1
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_))
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)+mhd_kin_en(w,ixImin1,ixImin2,ixImin3,&
                   ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                   ixOmax2,ixOmax3)+mhd_mag_en(w,ixImin1,ixImin2,ixImin3,&
                   ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                   ixOmax2,ixOmax3)
              end where
            else
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_))
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_) = small_e+mhd_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
                   ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
                   ixOmax3)+mhd_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
                   ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
                   ixOmax3)
              end where
            end if
          end if
        end if
      case ("average")
        ! do averaging of density
        call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
            flag, rho_)
        if(mhd_energy) then
          if(primitive)then
            call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
                flag, p_)
          else
            ! do averaging of internal energy
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               e_)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               e_)-mhd_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3)-mhd_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3)
            call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
                flag, e_)
            ! convert back
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               e_)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               e_)+mhd_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3)+mhd_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3)
          end if
        end if
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek-eb)
          if(mhd_energy) then
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               e_)-mhd_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
               ixOmax3)-mhd_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3))
          end if
          ! Convert momentum to velocity
          if(has_equi_rho0) then
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,0)
          else
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_)
          end if
          do idir = 1, ndir
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mom(idir))/tmp2(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do
        end if
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if

  end subroutine mhd_handle_small_values_origin

  subroutine mhd_handle_small_values_inte(primitive, w, x, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision :: tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call phys_check_w(primitive, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_rho0) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_) = small_density-block%equi_vars(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,0)
        else
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_) = small_density
        end if
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom(idir)) = 0.0d0
          end if
        end do
        if(mhd_energy) then
          if(primitive) then
            if(has_equi_pe0) then
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = small_pressure - &
                 block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,equi_pe0_,0)
               where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  e_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  p_) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            else
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 p_) = small_pressure
            end if
          else
            ! conserved
            if(has_equi_pe0) then
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = small_e - block%equi_vars(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)*inv_gamma_1
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_))
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_)=tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
              end where
            else
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_))
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=small_e
              end where
            end if
          end if
        end if
      case ("average")
        ! do averaging of density
        call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
            flag, rho_)
        if(mhd_energy) then
          if(primitive)then
            call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
                flag, p_)
          else
            ! do averaging of internal energy
            call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
                flag, e_)
          end if
        end if
      case default
        if(.not.primitive) then
          !convert w to primitive
          if(mhd_energy) then
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               p_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               e_)*gamma_1
          end if
          ! Convert momentum to velocity
          if(has_equi_rho0) then
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,0)
          else
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_)
          end if
          do idir = 1, ndir
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mom(idir))/tmp2(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do
        end if
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if

  end subroutine mhd_handle_small_values_inte

  subroutine mhd_handle_small_values_hde(primitive, w, x, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision :: tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call phys_check_w(primitive, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_) = small_density
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               rho_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom(idir)) = 0.0d0
          end if
        end do

        if(mhd_energy) then
          if(primitive) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               e_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               p_) = small_pressure
          else
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_))
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_) = small_e+mhd_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
                 ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
                 ixOmax3)
            end where
          end if
        end if
      case ("average")
        ! do averaging of density
        call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
            flag, rho_)
        if(mhd_energy) then
          if(primitive) then
            call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
                flag, p_)
          else
            ! do averaging of internal energy
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               e_)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               e_)-mhd_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3)
            call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
                flag, e_)
            ! convert back
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               e_)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               e_)+mhd_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3)
          end if
        end if
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek)
          if(mhd_energy) then
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               p_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               e_)-mhd_kin_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3))
          end if
          ! Convert momentum to velocity
          do idir = 1, ndir
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mom(idir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,rho_)
          end do
        end if
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if

  end subroutine mhd_handle_small_values_hde

  !> Calculate v vector
  subroutine mhd_get_v_origin(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndir)

    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: idir

    call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1.d0/rho(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    ! Convert momentum to velocity
    do idir = 1, ndir
       v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

  end subroutine mhd_get_v_origin

  !> Calculate v vector
  subroutine mhd_get_v_boris(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndir)

    double precision              :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer :: idir

    call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1.d0/rho(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    gamma2=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mag(:))**2,dim=ndim+1)*rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*inv_squared_c)
    ! Convert momentum to velocity
    do idir = 1, ndir
       v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*gamma2
    end do

  end subroutine mhd_get_v_boris

  !> Calculate v component
  subroutine mhd_get_v_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision              :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)

    if(mhd_boris_simplification) then
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idim)) / rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) /(1.d0+sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))**2,&
         dim=ndim+1)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*inv_squared_c)
    else
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idim)) / rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if

  end subroutine mhd_get_v_idim

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mhd_get_cmax_origin(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_csound(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,cmax)
    call mhd_get_v_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,vel)

    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=abs(vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))+cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine mhd_get_cmax_origin

  !> Calculate cmax_idim for semirelativistic MHD
  subroutine mhd_get_cmax_semirelati(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(inout):: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: wprim(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       nw)
    double precision :: csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), idim_Alfven_speed2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), Alfven_speed2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndir)

    if(B0field) then
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,b0i)
    else
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))
    end if
    inv_rho = 1.d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

    Alfven_speed2=sum(B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,&
       dim=ndim+1)*inv_rho
    gamma2 = 1.0d0/(1.d0+Alfven_speed2*inv_squared_c)

    wprim=w
    call mhd_to_primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wprim,x)
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.d0-gamma2*wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mom(idim))**2*inv_squared_c
    ! equation (69)
    Alfven_speed2=Alfven_speed2*cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    ! squared sound speed
    if(mhd_energy) then
      csound=mhd_gamma*wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)*inv_rho
    else
      csound=mhd_gamma*mhd_adiab*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)**gamma_1
    end if

    idim_Alfven_speed2=B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idim)**2*inv_rho

    ! Va_hat^2+a_hat^2 equation (57)
    Alfven_speed2=Alfven_speed2+csound*(1.d0+idim_Alfven_speed2*inv_squared_c)

    AvMinCs2=(gamma2*Alfven_speed2)**2&
       -4.0d0*gamma2*csound*idim_Alfven_speed2*cmax(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    where(AvMinCs2<zero)
       AvMinCs2=zero
    end where

    ! equation (68) fast magnetosonic wave speed
    csound = sqrt(half*(gamma2*Alfven_speed2+sqrt(AvMinCs2)))
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=gamma2*abs(wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mom(idim)))+csound

  end subroutine mhd_get_cmax_semirelati

  subroutine mhd_get_a2max(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a2max)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(inout) :: a2max(ndim)
    double precision :: a2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ndim,nw)
    integer :: gxOmin1,gxOmin2,gxOmin3,gxOmax1,gxOmax2,gxOmax3,hxOmin1,hxOmin2,&
       hxOmin3,hxOmax1,hxOmax2,hxOmax3,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,&
       jxOmax3,kxOmin1,kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3,i,j

    a2=zero
    do i = 1,ndim
      !> 4th order
      hxOmin1=ixOmin1-kr(i,1);hxOmin2=ixOmin2-kr(i,2);hxOmin3=ixOmin3-kr(i,3)
      hxOmax1=ixOmax1-kr(i,1);hxOmax2=ixOmax2-kr(i,2);hxOmax3=ixOmax3-kr(i,3);
      gxOmin1=hxOmin1-kr(i,1);gxOmin2=hxOmin2-kr(i,2);gxOmin3=hxOmin3-kr(i,3)
      gxOmax1=hxOmax1-kr(i,1);gxOmax2=hxOmax2-kr(i,2);gxOmax3=hxOmax3-kr(i,3);
      jxOmin1=ixOmin1+kr(i,1);jxOmin2=ixOmin2+kr(i,2);jxOmin3=ixOmin3+kr(i,3)
      jxOmax1=ixOmax1+kr(i,1);jxOmax2=ixOmax2+kr(i,2);jxOmax3=ixOmax3+kr(i,3);
      kxOmin1=jxOmin1+kr(i,1);kxOmin2=jxOmin2+kr(i,2);kxOmin3=jxOmin3+kr(i,3)
      kxOmax1=jxOmax1+kr(i,1);kxOmax2=jxOmax2+kr(i,2);kxOmax3=jxOmax3+kr(i,3);
      a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,&
         1:nw)=abs(-w(kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3,&
         1:nw)+16.d0*w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
         1:nw)-30.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nw)+16.d0*w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
         1:nw)-w(gxOmin1:gxOmax1,gxOmin2:gxOmax2,gxOmin3:gxOmax3,1:nw))
      a2max(i)=maxval(a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,&
         1:nw))/12.d0/dxlevel(i)**2
    end do
  end subroutine mhd_get_a2max

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine mhd_get_tcutoff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,Tco_local,Tmax_local)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(out) :: Tco_local,Tmax_local

    double precision, parameter :: trac_delta=0.25d0
    double precision :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       Te(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),lts(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir) :: bunitvec
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) :: gradT
    double precision :: Bdir(ndim)
    double precision :: ltrc,ltrp,altr(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer :: idims,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,ixA1,ixA2,ixA3,ixB1,ixB2,ixB3
    integer :: jxPmin1,jxPmin2,jxPmin3,jxPmax1,jxPmax2,jxPmax3,hxPmin1,hxPmin2,&
       hxPmin3,hxPmax1,hxPmax2,hxPmax3,ixPmin1,ixPmin2,ixPmin3,ixPmax1,ixPmax2,&
       ixPmax3,ixQmin1,ixQmin2,ixQmin3,ixQmax1,ixQmax2,ixQmax3
    logical :: lrlt(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_temperature(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,Te)
    Tco_local=zero
    Tmax_local=maxval(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

    
    
    select case(mhd_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      block%wextra(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         Tcoff_)=2.5d5/unit_temperature
    case(1,4,6)
      ! temperature gradient at cell centers
      do idims=1,ndim
        call gradient(Te,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,tmp1)
        gradT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do
      ! B vector
      if(B0field) then
        bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           :)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           iw_mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           :,0)
      else
        bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           :)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mag(:))
      end if
      if(mhd_trac_type .gt. 1) then
        ! B direction at cell center
        Bdir=zero
        do ixA1=0,1
    do ixA2=0,1
    do ixA3=0,1
          ixB1=(ixOmin1+ixOmax1-1)/2+ixA1;ixB2=(ixOmin2+ixOmax2-1)/2+ixA2
          ixB3=(ixOmin3+ixOmax3-1)/2+ixA3;
          Bdir(1:ndim)=Bdir(1:ndim)+bunitvec(ixB1,ixB2,ixB3,1:ndim)
        end do
    end do
    end do
        if(sum(Bdir(:)**2) .gt. zero) then
          Bdir(1:ndim)=Bdir(1:ndim)/dsqrt(sum(Bdir(:)**2))
        end if
        block%special_values(3:ndim+2)=Bdir(1:ndim)
      end if
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=dsqrt(sum(bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)**2,dim=ndim+1))
      where(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/=0.d0)
        tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=1.d0/tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      elsewhere
        tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=bigdouble
      end where
      ! b unit vector: magnetic field direction vector
      do idims=1,ndim
        bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)=bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do
      ! temperature length scale inversed
      lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=abs(sum(gradT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndim)*bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndim),dim=ndim+1))/Te(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      ! fraction of cells size to temperature length scale
      if(slab_uniform) then
        lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=minval(dxlevel)*lts(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      else
        lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=minval(block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,:),dim=ndim+1)*lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
      lrlt=.false.
      where(lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) > trac_delta)
        lrlt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=.true.
      end where
      if(any(lrlt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))) then
        block%special_values(1)=maxval(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3), mask=lrlt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
      else
        block%special_values(1)=zero
      end if
      block%special_values(2)=Tmax_local
    case(2)
      !> iijima et al. 2021, LTRAC method
      ltrc=1.5d0
      ltrp=4.d0
      ixPmin1=ixOmin1-2;ixPmin2=ixOmin2-2;ixPmin3=ixOmin3-2;ixPmax1=ixOmax1+2
      ixPmax2=ixOmax2+2;ixPmax3=ixOmax3+2;
      ! temperature gradient at cell centers
      do idims=1,ndim
        ixQmin1=ixPmin1;ixQmin2=ixPmin2;ixQmin3=ixPmin3;ixQmax1=ixPmax1
        ixQmax2=ixPmax2;ixQmax3=ixPmax3;
        hxPmin1=ixPmin1;hxPmin2=ixPmin2;hxPmin3=ixPmin3;hxPmax1=ixPmax1
        hxPmax2=ixPmax2;hxPmax3=ixPmax3;
        jxPmin1=ixPmin1;jxPmin2=ixPmin2;jxPmin3=ixPmin3;jxPmax1=ixPmax1
        jxPmax2=ixPmax2;jxPmax3=ixPmax3;
        select case(idims)
          case(1)
           ixQmin1=ixQmin1+1
           ixQmax1=ixQmax1-1
           hxPmax1=ixPmin1
           jxPmin1=ixPmax1
          
    case(2)
           ixQmin2=ixQmin2+1
           ixQmax2=ixQmax2-1
           hxPmax2=ixPmin2
           jxPmin2=ixPmax2
          
    case(3)
           ixQmin3=ixQmin3+1
           ixQmax3=ixQmax3-1
           hxPmax3=ixPmin3
           jxPmin3=ixPmax3
          
        end select
        call gradient(Te,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixQmin1,ixQmin2,ixQmin3,ixQmax1,ixQmax2,ixQmax3,idims,&
           gradT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims))
        call gradientx(Te,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           hxPmin1,hxPmin2,hxPmin3,hxPmax1,hxPmax2,hxPmax3,idims,&
           gradT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims),&
           .false.)
        call gradientq(Te,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           jxPmin1,jxPmin2,jxPmin3,jxPmax1,jxPmax2,jxPmax3,idims,&
           gradT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims))
      end do
      ! B vector
      if(B0field) then
        bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           :)=w(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           iw_mag(:))+block%B0(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           :,0)
      else
        bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           :)=w(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,iw_mag(:))
      end if
      tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3)=1.d0/(dsqrt(sum(bunitvec(ixPmin1:ixPmax1,&
         ixPmin2:ixPmax2,ixPmin3:ixPmax3,:)**2,dim=ndim+1))+smalldouble)
      ! b unit vector: magnetic field direction vector
      do idims=1,ndim
        bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           idims)=bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           idims)*tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)
      end do
      ! temperature length scale inversed
      lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3)=abs(sum(gradT(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3,1:ndim)*bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3,1:ndim),dim=ndim+1))/Te(ixPmin1:ixPmax1,&
         ixPmin2:ixPmax2,ixPmin3:ixPmax3)
      ! fraction of cells size to temperature length scale
      if(slab_uniform) then
        lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)=minval(dxlevel)*lts(ixPmin1:ixPmax1,&
           ixPmin2:ixPmax2,ixPmin3:ixPmax3)
      else
        lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)=minval(block%ds(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3,:),dim=ndim+1)*lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)
      end if
      lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)=max(one,&
          (exp(lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3))/ltrc)**ltrp)
  
      altr=zero
      ! need one ghost layer for thermal conductivity
      ixPmin1=ixOmin1-1;ixPmin2=ixOmin2-1;ixPmin3=ixOmin3-1;ixPmax1=ixOmax1+1
      ixPmax2=ixOmax2+1;ixPmax3=ixOmax3+1;
      do idims=1,ndim
        hxOmin1=ixPmin1-kr(idims,1);hxOmin2=ixPmin2-kr(idims,2)
        hxOmin3=ixPmin3-kr(idims,3);hxOmax1=ixPmax1-kr(idims,1)
        hxOmax2=ixPmax2-kr(idims,2);hxOmax3=ixPmax3-kr(idims,3);
        jxOmin1=ixPmin1+kr(idims,1);jxOmin2=ixPmin2+kr(idims,2)
        jxOmin3=ixPmin3+kr(idims,3);jxOmax1=ixPmax1+kr(idims,1)
        jxOmax2=ixPmax2+kr(idims,2);jxOmax3=ixPmax3+kr(idims,3);
        altr(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)=altr(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)+0.25d0*(lts(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3)+two*lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)+lts(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3))*bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3,idims)**2
      end do
      block%wextra(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
         Tcoff_)=Te(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3)*altr(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3)**0.4d0
    case(3,5)
      !> do nothing here
    case default
      call mpistop("unknown mhd_trac_type")
    end select
   
  end subroutine mhd_get_tcutoff

  !> get H speed for H-correction to fix the carbuncle problem at grid-aligned shock front
  subroutine mhd_get_H_speed(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,Hspeed)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wprim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    double precision :: csound(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ndim),tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3, ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, ixAmin1,ixAmin2,ixAmin3,&
       ixAmax1,ixAmax2,ixAmax3, id, ix1,ix2,ix3

    Hspeed=0.d0
    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    do id=1,ndim
      call mhd_get_csound_prim(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,id,tmp)
      csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)=tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)
    end do
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
    ixCmin1=ixOmin1+kr(idim,1)-1;ixCmin2=ixOmin2+kr(idim,2)-1
    ixCmin3=ixOmin3+kr(idim,3)-1;
    jxCmax1=ixCmax1+kr(idim,1);jxCmax2=ixCmax2+kr(idim,2)
    jxCmax3=ixCmax3+kr(idim,3);
    jxCmin1=ixCmin1+kr(idim,1);jxCmin2=ixCmin2+kr(idim,2)
    jxCmin3=ixCmin3+kr(idim,3);
    Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       1)=0.5d0*abs(wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
       mom(idim))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
       idim)-wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       mom(idim))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       idim))

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=ixCmax1+kr(id,1);ixAmax2=ixCmax2+kr(id,2)
      ixAmax3=ixCmax3+kr(id,3);
      ixAmin1=ixCmin1+kr(id,1);ixAmin2=ixCmin2+kr(id,2)
      ixAmin3=ixCmin3+kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)-wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(id))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,id)))
      ixAmax1=ixCmax1-kr(id,1);ixAmax2=ixCmax2-kr(id,2)
      ixAmax3=ixCmax3-kr(id,3);
      ixAmin1=ixCmin1-kr(id,1);ixAmin2=ixCmin2-kr(id,2)
      ixAmin3=ixCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom(id))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         id)-wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,id)))
    end do

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=jxCmax1+kr(id,1);ixAmax2=jxCmax2+kr(id,2)
      ixAmax3=jxCmax3+kr(id,3);
      ixAmin1=jxCmin1+kr(id,1);ixAmin2=jxCmin2+kr(id,2)
      ixAmin3=jxCmin3+kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)-wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         mom(id))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,id)))
      ixAmax1=jxCmax1-kr(id,1);ixAmax2=jxCmax2-kr(id,2)
      ixAmax3=jxCmax3-kr(id,3);
      ixAmin1=jxCmin1-kr(id,1);ixAmin2=jxCmin2-kr(id,2)
      ixAmin3=jxCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         mom(id))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         id)-wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,id)))
    end do

  end subroutine mhd_get_H_speed

  !> Estimating bounds for the minimum and maximum signal velocities without split
  subroutine mhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     Hspeed,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    double precision :: wmean(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       nw)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix1,ix2,ix3

    select case (boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_))
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_))
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call mhd_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+0.5d0*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim))-wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim)))**2
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=abs(umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (2)
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux))
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim))/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)
      call mhd_get_csound(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call mhd_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=max(csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))-csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    end select

  end subroutine mhd_get_cbounds

  !> Estimating bounds for the minimum and maximum signal velocities without split
  subroutine mhd_get_cbounds_semirelati(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3,idim,Hspeed,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) :: csoundL, csoundR, gamma2L, gamma2R

    ! Miyoshi 2005 JCP 208, 315 equation (67)
    call mhd_get_csound_semirelati(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
       csoundL,gamma2L)
    call mhd_get_csound_semirelati(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
       csoundR,gamma2R)
    csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=max(csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    if(present(cmin)) then
      cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)=min(gamma2L*wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim)),gamma2R*wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim)))-csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)=max(gamma2L*wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim)),gamma2R*wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)=max(gamma2L*wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim)),gamma2R*wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    end if

  end subroutine mhd_get_cbounds_semirelati

  !> Estimating bounds for the minimum and maximum signal velocities with rho split
  subroutine mhd_get_cbounds_split_rho(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3,idim,Hspeed,cmax,cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    double precision :: wmean(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       nw)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix1,ix2,ix3
    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    select case (boundspeed)
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,b0i))
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,b0i))
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call mhd_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+0.5d0*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**2*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim))-wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim)))**2
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=abs(umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (2)
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux))
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom(idim))/(wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,b0i))
      call mhd_get_csound(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call mhd_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=max(csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))-csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,1)=sign(one,cmin(ix1,ix2,ix3,1))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
            cmax(ix1,ix2,ix3,1)=sign(one,cmax(ix1,ix2,ix3,1))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,1))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    end select

  end subroutine mhd_get_cbounds_split_rho

  !> prepare velocities for ct methods
  subroutine mhd_get_ct_velocity(vcts,wLp,wRp,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,cmax,&
     cmin)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(in), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    type(ct_velocity), intent(inout):: vcts

    integer                         :: idimE,idimN

    ! calculate velocities related to different UCT schemes
    select case(type_ct)
    case('average')
    case('uct_contact')
      if(.not.allocated(vcts%vnorm)) allocate(vcts%vnorm(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1:ndim))
      ! get average normal velocity at cell faces
      vcts%vnorm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)=0.5d0*(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim)))
    case('uct_hll')
      if(.not.allocated(vcts%vbarC)) then
        allocate(vcts%vbarC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1:ndir,2),vcts%vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndir,2),vcts%vbarRC(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,1:ndir,2))
        allocate(vcts%cbarmin(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1:ndim),vcts%cbarmax(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndim))
      end if
      ! Store magnitude of characteristics
      if(present(cmin)) then
        vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)=max(-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           zero)
        vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)=max( cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           zero)
      else
        vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)=max( cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           zero)
        vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)=vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)
      end if

      idimN=mod(idim,ndir)+1 ! 'Next' direction
      idimE=mod(idim+1,ndir)+1 ! Electric field direction
      ! Store velocities
      vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         1)=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idimN))
      vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         1)=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idimN))
      vcts%vbarC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         1)=(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,1) +vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,1))/(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)+vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim))

      vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         2)=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idimE))
      vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         2)=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idimE))
      vcts%vbarC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         2)=(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,2) +vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,1))/(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)+vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim))
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine mhd_get_ct_velocity

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), b2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        kmax
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if(has_equi_rho0) then
      inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,b0i))
    else
      inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)
    endif

    call mhd_get_csound2(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound)

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = mhd_mag_en_all(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)

    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2-4.0d0*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * mhd_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       idim)**2 * inv_rho

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    if (.not. MHD_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = sqrt(half*(cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))
       if (mhd_boris_simplification) then
          ! equation (88)
          csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = mhd_gamma_alfven(w, ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
             ixOmax3) * csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
       end if
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),dxlevel(3),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3))),&
           mhd_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))*inv_rho*kmax)
    end if

  end subroutine mhd_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound_prim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), b2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        kmax
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)


    if(mhd_energy) then
      if(has_equi_pe0) then
        csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_) + block%equi_vars(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,b0i)
      else
        csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_)
      endif
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_gamma*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*inv_rho
    else
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_gamma*mhd_adiab*tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)**gamma_1
    end if

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)        = &
       mhd_mag_en_all(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2-4.0d0*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * mhd_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       idim)**2 * inv_rho

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    if (.not. MHD_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = sqrt(half*(cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))
       if (mhd_boris_simplification) then
          ! equation (88)
          csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = mhd_gamma_alfven(w, ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
             ixOmax3) * csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
       end if
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),dxlevel(3),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3))),&
           mhd_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))*inv_rho*kmax)
    end if

  end subroutine mhd_get_csound_prim

  !> Calculate cmax_idim for semirelativistic MHD
  subroutine mhd_get_csound_semirelati(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     csound,gamma2)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    ! here w is primitive variables
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out):: csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision :: AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), kmax
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), Alfven_speed2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), idim_Alfven_speed2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndir)

    integer :: ix1,ix2,ix3

    if(B0field) then
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,b0i)
    else
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))
    end if

    inv_rho = 1.d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

    Alfven_speed2=sum(B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,&
       dim=ndim+1)*inv_rho
    gamma2 = 1.0d0/(1.d0+Alfven_speed2*inv_squared_c)

    AvMinCs2=1.d0-gamma2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim))**2*inv_squared_c
    ! equatoin (69)
    Alfven_speed2=Alfven_speed2*AvMinCs2

    ! squared sound speed
    if(mhd_energy) then
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_gamma*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,p_)*inv_rho
    else
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_gamma*mhd_adiab*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)**gamma_1
    end if

    idim_Alfven_speed2=B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idim)**2*inv_rho

    ! Va_hat^2+a_hat^2 equation (57)
    Alfven_speed2=Alfven_speed2+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(1.d0+idim_Alfven_speed2*inv_squared_c)

    AvMinCs2=(gamma2*Alfven_speed2)**2-4.0d0*gamma2*csound(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*idim_Alfven_speed2*AvMinCs2

    where(AvMinCs2<zero)
       AvMinCs2=zero
    end where

    ! equation (68) fast magnetosonic speed
    csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = sqrt(half*(gamma2*Alfven_speed2+sqrt(AvMinCs2)))

  end subroutine mhd_get_csound_semirelati

  !> Calculate isothermal thermal pressure
  subroutine mhd_get_pthermal_iso(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=mhd_adiab*pth(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**mhd_gamma

  end subroutine mhd_get_pthermal_iso

  !> Calculate thermal pressure from internal energy
  subroutine mhd_get_pthermal_eint(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                      :: iw, ix1,ix2,ix3

    pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=gamma_1*w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)

    if(has_equi_pe0) then
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_pe0_,b0i)
    end if

    if (fix_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
            pth(ix1,ix2,ix3)=small_pressure
         end if
      enddo
      enddo
      enddo
    else if (check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix1,ix2,ix3),&
              " encountered when call mhd_get_pthermal"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix1,ix2,ix3,:)
           write(*,*) "Cell number: ", ix1,ix2,ix3
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix1,ix2,ix3,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) sqrt(pth(ix1,ix2,ix3)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      enddo
      enddo
      enddo
    end if

  end subroutine mhd_get_pthermal_eint

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine mhd_get_pthermal_origin(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                      :: iw, ix1,ix2,ix3

    pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=gamma_1*(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)- mhd_kin_en(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)- mhd_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3))

    if(has_equi_pe0) then
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_pe0_,b0i)
    end if

    if (fix_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
            pth(ix1,ix2,ix3)=small_pressure
         end if
      enddo
      enddo
      enddo
    else if (check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix1,ix2,ix3),&
              " encountered when call mhd_get_pthermal"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix1,ix2,ix3,:)
           write(*,*) "Cell number: ", ix1,ix2,ix3
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix1,ix2,ix3,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) sqrt(pth(ix1,ix2,ix3)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      enddo
      enddo
      enddo
    end if

  end subroutine mhd_get_pthermal_origin

  !> Calculate thermal pressure
  subroutine mhd_get_pthermal_semirelati(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                      :: iw, ix1,ix2,ix3

    double precision :: wprim(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       nw)

    wprim=w
    call mhd_to_primitive_semirelati(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wprim,x)
    pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=wprim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)

    if (.not.fix_small_values .and. check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix1,ix2,ix3),&
              " encountered when call mhd_get_pthermal_semirelati"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix1,ix2,ix3,:)
           write(*,*) "Cell number: ", ix1,ix2,ix3
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix1,ix2,ix3,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) sqrt(pth(ix1,ix2,ix3)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      enddo
      enddo
      enddo
    end if

  end subroutine mhd_get_pthermal_semirelati

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine mhd_get_pthermal_hde(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                      :: iw, ix1,ix2,ix3

    pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=gamma_1*(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-mhd_kin_en(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3))

    if (fix_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
            pth(ix1,ix2,ix3)=small_pressure
         end if
      enddo
      enddo
      enddo
    else if (check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(pth(ix1,ix2,ix3)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix1,ix2,ix3),&
              " encountered when call mhd_get_pthermal_hde"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix1,ix2,ix3,:)
           write(*,*) "Cell number: ", ix1,ix2,ix3
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix1,ix2,ix3,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) sqrt(pth(ix1,ix2,ix3)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      enddo
      enddo
      enddo
    end if

  end subroutine mhd_get_pthermal_hde

  !> copy temperature from stored Te variable
  subroutine mhd_get_temperature_from_Te(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3, Te_)
  end subroutine mhd_get_temperature_from_Te

  !> Calculate temperature=p/rho when in e_ the internal energy is stored
  subroutine mhd_get_temperature_from_eint(w, x, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: R(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_Rfactor(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,R)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = gamma_1 * &
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        e_)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)*R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
  end subroutine mhd_get_temperature_from_eint

  !> Calculate temperature=p/rho when in e_ the total energy is stored
  subroutine mhd_get_temperature_from_etot(w, x, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: R(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_Rfactor(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,R)
    call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=res(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_))

  end subroutine mhd_get_temperature_from_etot

  subroutine mhd_get_temperature_from_etot_with_equi(w, x, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: R(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_Rfactor(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,R)
    call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=res(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_rho0_,b0i)))

  end subroutine mhd_get_temperature_from_etot_with_equi

  subroutine mhd_get_temperature_from_eint_with_equi(w, x, ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: R(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_Rfactor(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,R)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = (gamma_1 * &
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        e_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_pe0_,b0i)) /((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_) +block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_rho0_,b0i))*R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

  end subroutine mhd_get_temperature_from_eint_with_equi

  subroutine mhd_get_temperature_equi(w,x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: R(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_Rfactor(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,R)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)= &
       block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_pe0_,b0i)/(block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,equi_rho0_,b0i)*R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

  end subroutine mhd_get_temperature_equi

  subroutine mhd_get_rho_equi(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_rho0_,b0i)
  end subroutine mhd_get_rho_equi

  subroutine mhd_get_pe_equi(w,x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, res)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out):: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_pe0_,b0i)
  end subroutine mhd_get_pe_equi

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine mhd_get_csound2(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision    :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    
    call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    if(mhd_energy) then
      call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_gamma*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    else
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_gamma*mhd_adiab*rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)**gamma_1
    end if
  end subroutine mhd_get_csound2

  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine mhd_get_p_total(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,p)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: p(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,p)

    p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = p(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) + 0.5d0 * sum(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(:))**2, dim=ndim+1)

  end subroutine mhd_get_p_total

  !> Calculate fluxes within ixO^L without any splitting
  subroutine mhd_get_flux(wC,w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nwflux)

    double precision             :: ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision             :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision             :: vHall(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:,:,:,:) :: Jambi, btot
    double precision, allocatable, dimension(:,:,:) :: tmp2, tmp3

    ! Get flux of density
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

    if(mhd_energy) then
      ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         p_)+0.5d0*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)
    else
      ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_adiab*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)**mhd_gamma+0.5d0*sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)
    end if

    if (mhd_Hall) then
      call mhd_getv_Hall(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vHall)
    end if

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))
      end if
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(mhd_energy) then
      if (mhd_internal_e) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))*(wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:)),&
           dim=ndim+1)
        if(mhd_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim) * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(:))**2,dim=ndim+1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idim)) * sum(vHall(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)),dim=ndim+1)
        end if
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              psi_)
        else
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))
        if (mhd_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir)) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir)) - vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idim)) + vHall(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)*w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))
        end if
      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(mhd_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * mhd_eta_ambi/rho**2
      allocate(Jambi(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:3))
      call mhd_get_Jambi(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Jambi)
      allocate(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3))
      btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1:3))
      allocate(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
         tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      !tmp2 = Btot^2
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)*btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
        case(2)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
        case(3)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
      endselect

      if(mhd_energy .and. .not. mhd_internal_e) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

  end subroutine mhd_get_flux

  !> Calculate fluxes within ixO^L without any splitting
  subroutine mhd_get_flux_hde(wC,w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nwflux)

    double precision             :: pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision             :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:,:,:,:) :: Jambi, btot
    double precision, allocatable, dimension(:,:,:) :: tmp2, tmp3

    ! Get flux of density
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    ! pgas is time dependent only
    if(mhd_energy) then
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
    else
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_adiab*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)**mhd_gamma
    end if

    ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+0.5d0*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))
      end if
    end do

    ! Get flux of energy
    if(mhd_energy) then
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*(wC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)+pgas(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              psi_)
        else
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))
      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(mhd_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * mhd_eta_ambi/rho**2
      allocate(Jambi(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:3))
      call mhd_get_Jambi(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Jambi)
      allocate(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3))
      btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1:3))
      allocate(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
         tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      !tmp2 = Btot^2
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)*btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
        case(2)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
        case(3)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
      endselect

      if(mhd_energy) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

  end subroutine mhd_get_flux_hde

  !> Calculate fluxes within ixO^L with possible splitting
  subroutine mhd_get_flux_split(wC,w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nwflux)

    double precision             :: pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndir)
    double precision             :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision             :: vHall(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir
    double precision, allocatable, dimension(:,:,:,:) :: Jambi, btot
    double precision, allocatable, dimension(:,:,:) :: tmp2, tmp3
    double precision :: tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)


    call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
    ! Get flux of density
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    ! pgas is time dependent only
    if(mhd_energy) then
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
    else
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_adiab*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**mhd_gamma
      if(has_equi_pe0) then
        pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe0_,b0i)
      endif
    end if

    if (mhd_Hall) then
      call mhd_getv_Hall(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vHall)
    end if

    if(B0field) then
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,idim)
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:,idim),dim=ndim+1)
    else
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))
    end if

    ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+0.5d0*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    if(B0field) then
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idim)-block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir,idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idim))
        else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)-block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir,idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idim))
        end if
      end do
    else
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idim)
        else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)
        end if
      end do
    end if

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(mhd_energy) then
      if (mhd_internal_e) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))*(wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))-B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:)),&
           dim=ndim+1)

        if (mhd_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
           if (mhd_etah>zero) then
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idim) * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  mag(:))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:),&
                 dim=ndim+1) - B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,idim) * sum(vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)),dim=ndim+1)
           end if
        end if
      end if
      if(has_equi_pe0) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=  f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim)) * block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe0_,idim) * inv_gamma_1
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              psi_)
        else
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)-B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))

        if (mhd_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (mhd_etah>zero) then
            f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(idir)) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(idir)) - vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,idir)*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,idim) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,idim)*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,idir)
          end if
        end if

      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idim))
    end if

    ! Contributions of ambipolar term in explicit scheme
    if(mhd_ambipolar_exp.and. .not.stagger_grid) then
      ! ambipolar electric field
      ! E_ambi=-eta_ambi*JxBxB=-JaxBxB=B^2*Ja-(Ja dot B)*B
      !Ja=eta_ambi*J=J * mhd_eta_ambi/rho**2
      allocate(Jambi(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:3))
      call mhd_get_Jambi(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Jambi)
      allocate(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3))
      if(B0field) then
        do idir=1,3
          btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir)) + block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir,idim)
        enddo
      else
        btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1:3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1:3))
      endif
      allocate(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
         tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      !tmp2 = Btot^2
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:3)**2,dim=ndim+1)
      !tmp3 = J_ambi dot Btot
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)*btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:),dim=ndim+1)

      select case(idim)
        case(1)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          if(B0field) tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
        case(2)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
          if(B0field) tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(3)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,3)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(3)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
        case(3)
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) *Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          if(B0field) tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(1)) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2))= f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(2)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * Jambi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1) - tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1)
      endselect

      if(mhd_energy .and. .not. mhd_internal_e) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_) +  tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *  tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif

      deallocate(Jambi,btot,tmp2,tmp3)
    endif

  end subroutine mhd_get_flux_split

  !> Calculate semirelativistic fluxes within ixO^L without any splitting
  subroutine mhd_get_flux_semirelati(wC,w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nwflux)

    double precision             :: pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision             :: SA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndir), B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)
    integer                      :: idirmin, iw, idir, jdir, kdir

    ! gas thermal pressure
    if(mhd_energy) then
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
    else
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mhd_adiab*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)**mhd_gamma
    end if

    ! Get flux of density
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))
    end do
    ! E=-uxB=Bxu
    if(B0field) then
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,idim)
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:,idim),dim=ndim+1)
    else
      B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))
    end if
    E=0.d0
    do idir=1,ndir; do jdir=1,ndir; do kdir=1,ndir
      if(lvc(idir,jdir,kdir)==1)then
        E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)+B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(kdir))
      else if(lvc(idir,jdir,kdir)==-1)then
        E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)-B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(kdir))
      end if
    end do; end do; end do

    pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)+half*(sum(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))**2,&
       dim=ndim+1)+sum(E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,&
       dim=ndim+1)*inv_squared_c)

    ! Get flux of momentum
    if(B0field) then
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))+pgas-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)-E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)*E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)*inv_squared_c-block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir,idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idim))
        else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)-E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)*E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)*inv_squared_c-block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir,idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(idim))
        end if
      end do
    else
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))+pgas-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)-E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)*E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)*inv_squared_c
        else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(idim))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)-E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)*E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)*inv_squared_c
        end if
      end do
    end if

    ! Get flux of total energy
    if(mhd_internal_e) then
      ! Get flux of internal energy
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*wC(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)
    else if(mhd_energy) then
      SA=0.d0
      do jdir=1,ndir; do kdir=1,ndir
        if(lvc(idim,jdir,kdir)==1)then
          SA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=SA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(kdir))
        else if(lvc(idim,jdir,kdir)==-1) then
          SA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=SA(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)-E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(kdir))
        end if
      end do; end do
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*(half*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)*sum(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(:))**2,&
         dim=ndim+1)+mhd_gamma*pgas*inv_gamma_1)+SA(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              psi_)
        else
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idim))*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)-B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))
      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idim))
    end if

  end subroutine mhd_get_flux_semirelati

  !> Source terms J.E in internal energy.
  !> For the ambipolar term E = ambiCoef * JxBxB=ambiCoef * B^2(-J_perpB)
  !=> the source term J.E = ambiCoef * B^2 * J_perpB^2 = ambiCoef * JxBxB^2/B^2
  !> ambiCoef is calculated as mhd_ambi_coef/rho^2,  see also the subroutine mhd_get_Jambi
  subroutine add_source_ambipolar_internal_energy(qdt,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     wCT,w,x,ie)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: jxbxb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)

    call mhd_get_jxbxb(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,jxbxb)
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sum(jxbxb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3)**2,&
       dim=ndim+1) / mhd_mag_en_all(wCT, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp,wCT,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie)+qdt * tmp

  end subroutine add_source_ambipolar_internal_energy

  subroutine mhd_get_jxbxb(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: res(:,:,:,:)

    double precision  :: btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)
    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       b2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    res=0.d0
    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    !!!here we know that current has nonzero values only for components in the range idirmin, 3
 
    if(B0field) then
      do idir=1,3
        btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir)) + block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir,b0i)
      enddo
    else
      btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1:3))
    endif

    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3)*btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3),dim=ndim+1) !J.B
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       sum(btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3)**2,&
       dim=ndim+1) !B2,!B2,!B2
    do idir=1,idirmin-1
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) = btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    enddo
    do idir=idirmin,3
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) = btot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) - current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idir) * b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    enddo
  end subroutine mhd_get_jxbxb

  !> Sets the sources for the ambipolar
  !> this is used for the STS method
  ! The sources are added directly (instead of fluxes as in the explicit)
  !> at the corresponding indices
  !>  store_flux_var is explicitly called for each of the fluxes one by one
  subroutine sts_set_source_ambipolar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
     fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,igrid,nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3) :: tmp,ff
    double precision :: fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nflux,1:ndim)
    double precision :: fE(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       sdim:3)
    double precision  :: btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3),tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: i, ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3, ie_

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;

    fluxall=zero

    call mhd_get_jxbxb(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmp)

    !set electric field in tmp: E=nuA * jxbxb, where nuA=-etaA/rho^2
    do i=1,3
      !tmp(ixA^S,i) = -(mhd_eta_ambi/w(ixA^S, rho_)**2) * tmp(ixA^S,i)
      call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,i),w,x)
    enddo

    if(mhd_energy .and. .not.mhd_internal_e) then
      !btot should be only mag. pert.
      btot(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:3)=0.d0
      !if(B0field) then
      !  do i=1,ndir
      !    btot(ixA^S, i) = w(ixA^S,mag(i)) + block%B0(ixA^S,i,0)
      !  enddo
      !else
        btot(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           1:ndir) = w(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           mag(1:ndir))
      !endif
      call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmp,btot,ff)
      call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      !- sign comes from the fact that the flux divergence is a source now
      wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    endif

    if(stagger_grid) then
      if(ndir>ndim) then
        !!!Bz
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           1) = tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,2)
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           2) = -tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1)
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,3) = 0.d0
        call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1+ndir,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndim)
        wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(ndir))=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if
      fE=0.d0
      call update_faces_ambipolar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,tmp,fE,&
         btot)
      ixAmax1=ixOmax1;ixAmax2=ixOmax2;ixAmax3=ixOmax3;
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;
      wres(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mag(1:ndim))=-btot(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         1:ndim)
    else
      !write curl(ele) as the divergence
      !m1={0,ele[[3]],-ele[[2]]}
      !m2={-ele[[3]],0,ele[[1]]}
      !m3={ele[[2]],-ele[[1]],0}

      !!!Bx
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1) = 0.d0
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         2) = tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,3)
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         3) = -tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,2)
      call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,2,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      !flux divergence is a source now
      wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1))=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      !!!By
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         1) = -tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,3)
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,2) = 0.d0
      ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         3) = tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1)
      call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
      if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,3,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(2))=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

      if(ndir==3) then
        !!!Bz
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           1) = tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,2)
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           2) = -tmp(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1)
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,3) = 0.d0
        call get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,tmp2)
        if(fix_conserve_at_step) fluxall(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1+ndir,1:ndim)=ff(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndim)
        wres(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(ndir))=-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if

    end if

    if(fix_conserve_at_step) then
      fluxall=my_dt*fluxall
      call store_flux(igrid,fluxall,1,ndim,nflux)
      if(stagger_grid) then
        call store_edge(igrid,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           my_dt*fE,1,ndim)
      end if
    end if

  end subroutine sts_set_source_ambipolar

  !> get ambipolar electric field and the integrals around cell faces
  subroutine update_faces_ambipolar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,ECC,fE,circ)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! amibipolar electric field at cell centers
    double precision, intent(in)       :: ECC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)
    double precision, intent(out)      :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)
    double precision, intent(out)      :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixAmin1,&
       ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3
    integer                            :: idim1,idim2,idir,ix1,ix2,ix3

    fE=zero
    ! calcuate ambipolar electric field on cell edges from cell centers
    do idir=sdim,3
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
      ixCmin3=ixOmin3+kr(idir,3)-1;
     do ix3=0,1
     do ix2=0,1
     do ix1=0,1
        if( ix1==1 .and. 1==idir  .or. ix2==1 .and. 2==idir  .or. ix3==1 .and. &
           3==idir ) cycle
        ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;ixAmin3=ixCmin3+ix3;
        ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;ixAmax3=ixCmax3+ix3;
        fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)+ECC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,idir)
     end do
     end do
     end do
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)*0.25d0*block%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,idir)
    end do

    ! Calculate circulation on each face to get value of line integral of
    ! electric field in the positive idir direction.
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
    ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;

    circ=zero

    do idim1=1,ndim ! Coordinate perpendicular to face
      do idim2=1,ndim
        do idir=sdim,3 ! Direction of line integral
          ! Assemble indices
          hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
          hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
          hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
          ! Add line integrals in direction idir
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
             hxCmin3:hxCmax3,idir))
        end do
      end do
      circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)/block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)
    end do

  end subroutine update_faces_ambipolar

  !> use cell-center flux to get cell-face flux
  !> and get the source term as the divergence of the flux
  subroutine get_flux_on_cell_face(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ff,src)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, dimension(:,:,:,:), intent(inout) :: ff
    double precision, intent(out) :: src(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: ffc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: dxinv(ndim)
    integer :: idims, ix1,ix2,ix3, ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,&
       ixAmax3, ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3, ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    dxinv=1.d0/dxlevel
    ! cell corner flux in ffc
    ffc=0.d0
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=ixOmin1-1
    ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;
    do ix3=0,1
    do ix2=0,1
    do ix1=0,1
      ixBmin1=ixCmin1+ix1;ixBmin2=ixCmin2+ix2;ixBmin3=ixCmin3+ix3;
      ixBmax1=ixCmax1+ix1;ixBmax2=ixCmax2+ix2;ixBmax3=ixCmax3+ix3;
      ffc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:ndim)=ffc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1:ndim)+ff(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1:ndim)
    end do
    end do
    end do
    ffc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       1:ndim)=0.5d0**ndim*ffc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       1:ndim)
    ! flux at cell face
    ff(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=0.d0
    do idims=1,ndim
      ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
      ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
      ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=ixBmin1
      ixCmin2=ixBmin2;ixCmin3=ixBmin3;
      do ix3=0,1 
      do ix2=0,1 
      do ix1=0,1 
         if( ix1==0 .and. 1==idims  .or. ix2==0 .and. 2==idims  .or. ix3==0 &
            .and. 3==idims ) then
           ixBmin1=ixCmin1-ix1;ixBmin2=ixCmin2-ix2;ixBmin3=ixCmin3-ix3;
           ixBmax1=ixCmax1-ix1;ixBmax2=ixCmax2-ix2;ixBmax3=ixCmax3-ix3;
           ff(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims)=ff(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idims)+ffc(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
              idims)
         end if
      end do
      end do
      end do
      ff(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idims)=ff(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idims)*0.5d0**(ndim-1)
    end do
    src=0.d0
    if(slab_uniform) then
      do idims=1,ndim
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idims)=dxinv(idims)*ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+ff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)-ff(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,idims)
      end do
    else
      do idims=1,ndim
        ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idims)=ff(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idims)*block%surfaceC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)
        ixBmin1=ixOmin1-kr(idims,1);ixBmin2=ixOmin2-kr(idims,2)
        ixBmin3=ixOmin3-kr(idims,3);ixBmax1=ixOmax1-kr(idims,1)
        ixBmax2=ixOmax2-kr(idims,2);ixBmax3=ixOmax3-kr(idims,3);
        src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+ff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idims)-ff(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,idims)
      end do
      src(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=src(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%dvolume(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if
  end subroutine get_flux_on_cell_face

  !> Calculates the explicit dt for the ambipokar term
  !> This function is used by both explicit scheme and STS method
  function get_ambipolar_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,&
     x)  result(dtnew)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: dtnew

    double precision              :: coef
    double precision              :: dxarr(ndim)
    double precision              :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    dxarr(1)=dx1;dxarr(2)=dx2;dxarr(3)=dx3;
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = mhd_mag_en_all(w,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp,w,x)
    coef = maxval(abs(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))
    if(coef/=0.d0) then
      coef=1.d0/coef
    else
      coef=bigdouble
    end if
    if(slab_uniform) then
      dtnew=minval(dxarr(1:ndim))**2.0d0*coef
    else
      dtnew=minval(block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndim))**2.0d0*coef
    end if

  end function get_ambipolar_dt

  !> multiply res by the ambipolar coefficient
  !> The ambipolar coefficient is calculated as -mhd_eta_ambi/rho^2
  !> The user may mask its value in the user file
  !> by implemneting usr_mask_ambipolar subroutine
  subroutine multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    tmp=0.d0
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=-&
       mhd_eta_ambi/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2
    if (associated(usr_mask_ambipolar)) then
      call usr_mask_ambipolar(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,tmp)
    end if

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
  end subroutine multiplyAmbiCoef

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine mhd_add_source(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
     wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source
    use mod_cak_force, only: cak_add_source

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    !TODO local_timestep support is only added for splitting
    ! but not for other nonideal terms such gravity, RC, viscosity,..
    ! it will also only work for divbfix  'linde', which does not require
    ! modification as it does not use dt in the update

    if (.not. qsourcesplit) then
      if(mhd_internal_e) then
        ! Source for solving internal energy
        active = .true.
        call add_source_internal_e(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
           wCTprim)
      else
        if(has_equi_pe0) then
          active = .true.
          call add_pe0_divv(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,&
             ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
             wCT,w,x)
        end if
      end if

      ! Source for B0 splitting
      if (B0field) then
        active = .true.
        call add_source_B0split(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
           w,x)
      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(mhd_eta)>smalldouble)then
        active = .true.
        call add_source_res2(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      end if

      if (mhd_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      end if

      if(mhd_hydrodynamic_e) then
        ! Source for solving hydrodynamic energy
        active = .true.
        call add_source_hydrodynamic_e(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
           w,x,wCTprim)
      else if (mhd_semirelativistic) then
        ! add sources for semirelativistic MHD
        active = .true.
        call add_source_semirelativistic(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
           w,x,wCTprim)
      end if
    end if

    
    if(source_split_divb .eqv. qsourcesplit) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_ct)
        continue ! Do nothing
      case (divb_linde)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_glm)
        active = .true.
        call add_source_glm(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
        call add_source_janhunen(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
        call add_source_powel(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
        call add_source_glm(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      case (divb_multigrid)
        continue ! Do nothing
      case (divb_none)
        ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
   

    if(mhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
         wCTprim,w,x,qsourcesplit,active, rc_fl)
    end if

    if(mhd_viscosity) then
      call viscosity_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
         mhd_energy,qsourcesplit,active)
    end if

    if(mhd_gravity) then
      call gravity_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,wCTprim,w,&
         x,gravity_energy,gravity_rhov,qsourcesplit,active)
    end if

    if (mhd_cak_force) then
      call cak_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,mhd_energy,&
         qsourcesplit,active)
    end if

    ! update temperature from new pressure, density, and old ionization degree
    if(mhd_partial_ionization) then
      if(.not.qsourcesplit) then
        active = .true.
        call mhd_update_temperature(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      end if
    end if

  end subroutine mhd_add_source

  subroutine add_pe0_divv(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    double precision                :: divv(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call mhd_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,v)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,&
           sixthorder=.true.)
      else
        call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,&
           fourthorder=.true.)
      end if
    else
     call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
        ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv)
    end if
    if(local_timestep) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-dtfactor*block%dt(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)*block%equi_vars(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe0_,0)*divv(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    else
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_)-qdt*block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_pe0_,0)*divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    endif
  end subroutine add_pe0_divv

  !> Compute the Lorentz force (JxB)
  subroutine get_Lorentz_force(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,JxB)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: JxB(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3)
    double precision                :: a(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3), b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       3)
    integer                         :: idir, idirmin
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)

    b=0.0d0
    do idir = 1, ndir
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir) = mhd_mag_i_all(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir)
    end do

    ! store J current in a
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)

    a=0.0d0
    do idir=7-2*ndir,3
      a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
    end do

    call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,JxB)
  end subroutine get_Lorentz_force

  !> Compute 1/(1+v_A^2/c^2) for semirelativistic MHD, where v_A is the Alfven
  !> velocity
  subroutine mhd_gamma2_alfven(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, gamma_A2)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(out) :: gamma_A2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision              :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    ! mhd_get_rho cannot be used as x is not a param
    if(has_equi_rho0) then
      rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,b0i)
    else
      rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    endif
    ! Compute the inverse of 1 + B^2/(rho * c^2)
    gamma_A2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1.0d0/(1.0d0+mhd_mag_en_all(w, ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*inv_squared_c)
  end subroutine mhd_gamma2_alfven

  !> Compute 1/sqrt(1+v_A^2/c^2) for semirelativisitic MHD, where v_A is the
  !> Alfven velocity
  function mhd_gamma_alfven(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(gamma_A)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: gamma_A(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    call mhd_gamma2_alfven(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, gamma_A)
    gamma_A = sqrt(gamma_A)
  end function mhd_gamma_alfven

  subroutine mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(has_equi_rho0) then
      rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho0_,b0i)
    else
      rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    endif

  end subroutine mhd_get_rho

  !> handle small or negative internal energy
  subroutine mhd_handle_small_ei(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ie, subname)
    use mod_global_parameters
    use mod_small_values
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ie
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    character(len=*), intent(in)    :: subname

    integer :: idir
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision              :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    flag=.false.
    if(has_equi_pe0) then
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         equi_pe0_,0)*inv_gamma_1<small_e)flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,ie)=.true.
    else
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)<small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)=.true.
    endif
    if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie))) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_pe0) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)=small_e - block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe0_,0)*inv_gamma_1
        else
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)=small_e
        endif
      case ("average")
        call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
            flag, ie)
      case default
        ! small values error shows primitive variables
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)*gamma_1
        call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
        do idir = 1, ndir
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mom(idir))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if

  end subroutine mhd_handle_small_ei

  subroutine mhd_update_temperature(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_ionization_degree

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: iz_H(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       iz_He(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        pth(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call ionization_degree_from_temperature(ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,Te_),iz_H,iz_He)

    call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       Te_)=(2.d0+3.d0*He_abundance)*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)*(1.d0+iz_H(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+He_abundance*(iz_He(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(iz_He(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+1.d0)+1.d0)))

  end subroutine mhd_update_temperature

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qdt, dtfactor,wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: a(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3),&
        b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3),&
        axb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3)
    integer :: idir

    a=0.d0
    b=0.d0
    ! for force-free field J0xB0 =0
    if(.not.B0field_forcefree) then
      ! store B0 magnetic field in b
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)=block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir,0)

      ! store J0 current in a
      do idir=7-2*ndir,3
        a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)
      end do
      call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,axb)
      if(local_timestep) then
        do idir=1,3
          axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)*block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*dtfactor
        enddo
      else
        axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           :)=axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*qdt
      endif
      ! add J0xB0 source term in momentum equations
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1:ndir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1:ndir))+axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)
    end if

    if(total_energy) then
      a=0.d0
      ! for free-free field -(vxB0) dot J0 =0
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)=wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))
      ! store full magnetic field B0+B1 in b
      if(.not.B0field_forcefree) b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)=b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0)
      ! store velocity in a
      call mhd_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1:ndir))
      call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,axb)
      if(local_timestep) then
        do idir=1,3
          axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)*block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*dtfactor
        enddo
      else
        axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           :)=axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*qdt
      endif
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)-axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)*block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)
      end do
      if(mhd_ambipolar) then
        !reuse axb
        call mhd_get_jxbxb(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,axb)
        ! source J0 * E
        do idir=sdim,3
          !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
          call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
             axb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idir),wCT,x)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)+axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)*block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)
        enddo
      endif
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_B0')

  end subroutine add_source_B0split

  !> Source terms for semirelativistic MHD Gombosi 2002 JCP 177, 176
  subroutine add_source_semirelativistic(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
     wCTprim)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in), optional :: wCTprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    double precision :: B(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3),&
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3),&
        E(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3),&
        divE(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: idir, idirmin

    ! store B0 magnetic field in b
    B=0.d0
    if(B0field) then
      B(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir,0)
    else
      B(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndir))
    end if
    
    
    call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,B,wCTprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,mom(1:ndir)),E)
   
    call divvector(E,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divE)
    ! curl E => B
    call curlvector(E,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,B,idirmin,1,3)
    ! E x (curl E) => v
    call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,E,B,v)
    ! add source term in momentum equations (1/c0^2-1/c^2)(E dot divE - E x curlE)
    ! equation (26) and (27)
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))+qdt*(inv_squared_c0-inv_squared_c)*(E(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)*divE(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)-v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idir))
    end do

  end subroutine add_source_semirelativistic

  !> Source terms for internal energy version of MHD
  subroutine add_source_internal_e(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,wCTprim)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision                :: divv(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           mom(:)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
           ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,sixthorder=.true.)
      else
        call divvector(wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           mom(:)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
           ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv,fourthorder=.true.)
      end if
    else
      call divvector(wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mom(:)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divv)
    end if
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-qdt*wCTprim(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)*divv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    if(mhd_ambipolar)then
      call add_source_ambipolar_internal_energy(qdt,ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3,wCT,w,x,e_)
    end if

    if(fix_small_values) then
      call mhd_handle_small_ei(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,e_,&
         'add_source_internal_e')
    end if
  end subroutine add_source_internal_e

  !> Source terms for hydrodynamic energy version of MHD
  subroutine add_source_hydrodynamic_e(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
     wCTprim)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_gravity

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in), optional :: wCTprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    double precision :: B(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3),&
        J(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3),&
        JxB(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3)
    integer :: idir, idirmin, idims, ix1,ix2,ix3
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision :: bu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:ndir), tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir), Vaoc

    
    
 !call get_current(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    ! get current in fourth order accuracy in Cartesian
    call curlvector(wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       mag(1:ndir)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,current,idirmin,1,ndir,.true.)
    ! get Lorentz force JxB
    call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,current,wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,mag(1:ndir)),JxB)
   

    if(mhd_semirelativistic) then
      ! (v . nabla) v
      do idir=1,ndir
        do idims=1,ndim
          call gradient(wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,mom(idir)),ixImin1,ixImin2,ixImin3,ixImax1,&
             ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
             idims,J(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idims))
        end do
        B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=sum(wCTprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(1:ndir))*J(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1:ndir),dim=ndim+1)
      end do
      ! nabla p
      do idir=1,ndir
        call gradient(wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           p_),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
           ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,J(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,idir))
      end do

      if(mhd_gravity) then
        gravity_field=0.d0
        call usr_gravity(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,x,&
           gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1:ndim))
        do idir=1,ndir
          B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)*(B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)-gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir))+J(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir)-JxB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir)
        end do
      else
        do idir=1,ndir
          B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_)*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)+J(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)-JxB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
        end do
      end if

      b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sum(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)>smalldouble)
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=1.d0/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      else where
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=0.d0
      end where
      ! unit vector of magnetic field
      do idir=1,ndir
        bu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do

      !b2(ixO^S)=b2(ixO^S)/w(ixO^S,rho_)*inv_squared_c
      !b2(ixO^S)=b2(ixO^S)/(1.d0+b2(ixO^S))
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
         ! Va^2/c^2
         Vaoc=b2(ix1,ix2,ix3)/w(ix1,ix2,ix3,rho_)*inv_squared_c
         ! Va^2/c^2 / (1+Va^2/c^2)
         b2(ix1,ix2,ix3)=Vaoc/(1.d0+Vaoc)
      end do
      end do
      end do
      ! bu . F
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sum(bu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir)*B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir),dim=ndim+1)
      ! Rempel 2017 ApJ 834, 10 equation (54)
      do idir=1,ndir
        J(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)=b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)-bu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      end do
      !! Rempel 2017 ApJ 834, 10 equation (29) add SR force at momentum equation
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))+qdt*J(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)
      end do
      ! Rempel 2017 ApJ 834, 10 equation (30) add work of Lorentz force and SR force
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)+qdt*sum(wCTprim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1:ndir))*(JxB(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)+J(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)),dim=ndim+1)
    else
      ! add work of Lorentz force
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)+qdt*sum(wCTprim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1:ndir))*JxB(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir),dim=ndim+1)
    end if

  end subroutine add_source_hydrodynamic_e

  !> Add resistive source to w within ixO Uses 3 point stencil (1 neighbour) in
  !> each direction, non-conservative. If the fourthorder precompiler flag is
  !> set, uses fourth order central difference for the laplacian. Then the
  !> stencil is 5 (2 neighbours).
  subroutine add_source_res1(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idir,jdir,kdir,&
       idirmin,idim,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,ix
    integer :: lxOmin1,lxOmin2,lxOmin3,lxOmax1,lxOmax2,lxOmax3, kxOmin1,&
       kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3

    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: gradeta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), Bf(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)

    ! Calculating resistive sources involve one extra layer
    if (mhd_4th_order) then
      ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmin3=ixOmin3-2;ixAmax1=ixOmax1+2
      ixAmax2=ixOmax2+2;ixAmax3=ixOmax3+2;
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
      ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    end if

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2.or.ixImin3>ixAmin3.or.ixImax3<ixAmax3) call &
       mpistop("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)

    if (mhd_eta>zero)then
       eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)=mhd_eta
       gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImin3,ixImax1,&
          ixImax2,ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,&
          idirmin,x,current,eta)
       ! assumes that eta is not function of current?
       do idim=1,ndim
          call gradient(eta,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,tmp)
          gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idim)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       end do
    end if

    if(B0field) then
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir,0)
    else
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndir))
    end if

    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (mhd_4th_order) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3)=Bf(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,idir)
         do idim=1,ndim
            lxOmin1=ixOmin1+2*kr(idim,1);lxOmin2=ixOmin2+2*kr(idim,2)
            lxOmin3=ixOmin3+2*kr(idim,3);lxOmax1=ixOmax1+2*kr(idim,1)
            lxOmax2=ixOmax2+2*kr(idim,2);lxOmax3=ixOmax3+2*kr(idim,3);
            jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
            jxOmin3=ixOmin3+kr(idim,3);jxOmax1=ixOmax1+kr(idim,1)
            jxOmax2=ixOmax2+kr(idim,2);jxOmax3=ixOmax3+kr(idim,3);
            hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
            hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
            hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
            kxOmin1=ixOmin1-2*kr(idim,1);kxOmin2=ixOmin2-2*kr(idim,2)
            kxOmin3=ixOmin3-2*kr(idim,3);kxOmax1=ixOmax1-2*kr(idim,1)
            kxOmax2=ixOmax2-2*kr(idim,2);kxOmax3=ixOmax3-2*kr(idim,3);
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+(-tmp2(lxOmin1:lxOmax1,lxOmin2:lxOmax2,&
               lxOmin3:lxOmax3)+16.0d0*tmp2(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3)-30.0d0*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+16.0d0*tmp2(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3)-tmp2(kxOmin1:kxOmax1,kxOmin2:kxOmax2,&
               kxOmin3:kxOmax3)) /(12.0d0 * dxlevel(idim)**2)
         end do
       else
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3)=Bf(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,idir)
         do idim=1,ndim
            jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
            jxOmin3=ixOmin3+kr(idim,3);jxOmax1=ixOmax1+kr(idim,1)
            jxOmax2=ixOmax2+kr(idim,2);jxOmax3=ixOmax3+kr(idim,3);
            hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
            hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
            hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+(tmp2(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
               jxOmin3:jxOmax3)-2.0d0*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+tmp2(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
               hxOmin3:hxOmax3))/dxlevel(idim)**2
         end do
       end if

       ! Multiply by eta
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)

       ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
       if (mhd_eta<zero)then
          do jdir=1,ndim; do kdir=idirmin,3
             if (lvc(idir,jdir,kdir)/=0)then
                if (lvc(idir,jdir,kdir)==1)then
                   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3)-gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3,jdir)*current(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
                else
                   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3)+gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3,jdir)*current(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir)
                end if
             end if
          end do; end do
       end if

       ! Add sources related to eta*laplB-grad(eta) x J to B and e
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mag(idir))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       if(total_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*Bf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir)
       end if
    end do ! idir

    if(mhd_energy) then
      ! de/dt+=eta*J**2
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)+qdt*eta(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)*sum(current(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,dim=ndim+1)
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_res1')

  end subroutine add_source_res1

  !> Add resistive source to w within ixO
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_res2(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),curlj(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)
    double precision :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3),tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idir,idirmin,&
       idirmin1

    ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmin3=ixOmin3-2;ixAmax1=ixOmax1+2
    ixAmax2=ixOmax2+2;ixAmax3=ixOmax3+2;

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2.or.ixImin3>ixAmin3.or.ixImax3<ixAmax3) call &
       mpistop("Error in add_source_res2: Non-conforming input limits")

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,current)

    tmpvec=zero
    if(mhd_eta>zero)then
      do idir=idirmin,3
        tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)=current(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)*mhd_eta
      end do
    else
      call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,x,&
         current,eta)
      do idir=idirmin,3
        tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)=current(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)*eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)
      end do
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    call curlvector(tmpvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,curlj,idirmin1,1,3)
    if(stagger_grid) then
      if(ndim==2.and.ndir==3) then
        ! if 2.5D
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(ndir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(ndir))-qdt*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,ndir)
      end if
    else
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))-qdt*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir)
    end if

    if(mhd_energy) then
      if(mhd_eta>zero)then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=qdt*mhd_eta*sum(current(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)**2,dim=ndim+1)
      else
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=qdt*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,:)**2,dim=ndim+1)
      end if
      if(total_energy) then
        ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
        ! de1/dt= eta J^2 - B1 dot curl(eta J)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-qdt*sum(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mag(1:ndir))*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1:ndir),dim=ndim+1)
      else
        ! add eta*J**2 source term in the internal energy equation
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_res2')
  end subroutine add_source_res2

  !> Add Hyper-resistive source to w within ixO
  !> Uses 9 point stencil (4 neighbours) in each direction.
  subroutine add_source_hyperres(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    !.. local ..
    double precision                :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision                :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3),tmpvec2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3),tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),ehyper(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)
    integer                         :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,&
       ixAmax3,idir,jdir,kdir,idirmin,idirmin1

    ixAmin1=ixOmin1-3;ixAmin2=ixOmin2-3;ixAmin3=ixOmin3-3;ixAmax1=ixOmax1+3
    ixAmax2=ixOmax2+3;ixAmax3=ixOmax3+3;
    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2.or.ixImin3>ixAmin3.or.ixImax3<ixAmax3) call &
       mpistop("Error in add_source_hyperres: Non-conforming input limits")

    call get_current(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,current)
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
          jdir)=current(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,jdir)
    end do

    ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmin3=ixOmin3-2;ixAmax1=ixOmax1+2
    ixAmax2=ixOmax2+2;ixAmax3=ixOmax3+2;
    call curlvector(tmpvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmpvec2,idirmin1,1,3)

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndir)=zero
    call curlvector(tmpvec2,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmpvec,idirmin1,1,3)
    ehyper(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
       1:ndir) = - tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
       1:ndir)*mhd_eta_hyper

    ixAmin1=ixOmin1;ixAmin2=ixOmin2;ixAmin3=ixOmin3;ixAmax1=ixOmax1
    ixAmax2=ixOmax2;ixAmax3=ixOmax3;
    tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndir)=zero
    call curlvector(ehyper,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,tmpvec2,idirmin1,1,3)

    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))-tmpvec2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)*qdt
    end do

    if(total_energy) then
      ! de/dt= +div(B x Ehyper)
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
      ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
      tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndir)=zero
      do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
        tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir) = tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
           idir)+ lvc(idir,jdir,kdir)*wCT(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,mag(jdir))*ehyper(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,kdir)
      end do; end do; end do
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
      call divvector(tmpvec2,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)+tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)*qdt
    end if

    if (fix_small_values)  call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_hyperres')

  end subroutine add_source_hyperres

  subroutine add_source_glm(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme or GLM-MHD scheme
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision:: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer          :: idim,idir
    double precision :: gradPsi(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)


    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mhd_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_) = abs(mhd_glm_alpha)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dxlevel(:)))*w(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(block%ds(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:),&
           dim=ndim+1))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           psi_)
      end if
    end if

    if(mhd_glm_extended) then
      ! gradient of Psi
      if(total_energy) then
        do idim=1,ndim
          select case(typegrad)
          case("central")
            call gradient(wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               psi_),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
               ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,gradPsi)
          case("limited")
            call gradientS(wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               psi_),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
               ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,gradPsi)
          end select
          ! e  = e  -qdt (b . grad(Psi))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idim))*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        end do
      end if

      ! We calculate now div B
      call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb,&
          mhd_divb_4thorder)

      ! m = m - qdt b div b
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))-qdt*mhd_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do
    end if

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_glm')

  end subroutine add_source_glm

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb, mhd_divb_4thorder)

    ! calculate velocity
    call mhd_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)

    if (total_energy) then
      ! e = e - qdt (v . b) * div b
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-qdt*sum(v(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)),&
         dim=ndim+1)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))-qdt*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))-qdt*mhd_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb, mhd_divb_4thorder)

    call mhd_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vel)
    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))-qdt*vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_janhunen')

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    ! Add Linde's divB related sources to wnew within ixO
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer :: idim, idir, ixpmin1,ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3, i1,&
       i2,i3, iside
    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       graddivb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    logical, dimension(-1:1,-1:1,-1:1) :: leveljump

    ! Calculate div B
    ixpmin1=ixOmin1-1;ixpmin2=ixOmin2-1;ixpmin3=ixOmin3-1;ixpmax1=ixOmax1+1
    ixpmax2=ixOmax2+1;ixpmax3=ixOmax3+1;
    call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixpmin1,&
       ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3,divb, mhd_divb_4thorder)

    ! for AMR stability, retreat one cell layer from the boarders of level jump
    do i3=-1,1
    do i2=-1,1
    do i1=-1,1
      if(i1==0.and.i2==0.and.i3==0) cycle
      if(neighbor_type(i1,i2,i3,block%igrid)==2 .or. neighbor_type(i1,i2,i3,&
         block%igrid)==4) then
        leveljump(i1,i2,i3)=.true.
      else
        leveljump(i1,i2,i3)=.false.
      end if
    end do
    end do
    end do

    ixpmin1=ixOmin1;ixpmin2=ixOmin2;ixpmin3=ixOmin3;ixpmax1=ixOmax1
    ixpmax2=ixOmax2;ixpmax3=ixOmax3;
    do idim=1,ndim
      select case(idim)
       case(1)
          do iside=1,2
            i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3)
            i3=kr(3,1)*(2*iside-3);
            if (leveljump(i1,i2,i3)) then
              if (iside==1) then
                ixpmin1=ixOmin1-i1
              else
                ixpmax1=ixOmax1-i1
              end if
            end if
          end do
       
       case(2)
          do iside=1,2
            i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3)
            i3=kr(3,2)*(2*iside-3);
            if (leveljump(i1,i2,i3)) then
              if (iside==1) then
                ixpmin2=ixOmin2-i2
              else
                ixpmax2=ixOmax2-i2
              end if
            end if
          end do
       
       case(3)
          do iside=1,2
            i1=kr(1,3)*(2*iside-3);i2=kr(2,3)*(2*iside-3)
            i3=kr(3,3)*(2*iside-3);
            if (leveljump(i1,i2,i3)) then
              if (iside==1) then
                ixpmin3=ixOmin3-i3
              else
                ixpmax3=ixOmax3-i3
              end if
            end if
          end do
       
      end select
    end do

    ! Add Linde's diffusive terms
    do idim=1,ndim
       ! Calculate grad_idim(divb)
       select case(typegrad)
       case("central")
         call gradient(divb,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixpmin1,ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3,idim,graddivb)
       case("limited")
         call gradientS(divb,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixpmin1,ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3,idim,graddivb)
       end select

       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab_uniform) then
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)=graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)*divbdiff/(1.0d0/dxlevel(1)**2+&
             1.0d0/dxlevel(2)**2+1.0d0/dxlevel(3)**2)
       else
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)=graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)*divbdiff /(1.0d0/block%ds(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
             1)**2+1.0d0/block%ds(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3,2)**2+1.0d0/block%ds(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,ixpmin3:ixpmax3,3)**2)
       end if

       w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
          mag(idim))=w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
          mag(idim))+graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3)

       if (typedivbdiff=='all' .and. total_energy) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
            e_)=w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
            e_)+wCT(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
            mag(idim))*graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
            ixpmin3:ixpmax3)
       end if
    end do

    if (fix_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_linde')

  end subroutine add_source_linde


  !> get dimensionless div B = |divB| * volume / area / |B|
  subroutine get_normalized_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)

    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                   :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), dsurface(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idims

    call get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
    invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(mhd_mag_en_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3))
    where(invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/=0.d0)
      invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    end where
    if(slab_uniform) then
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=0.5d0*abs(divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/sum(1.d0/dxlevel(:))
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;
      ixAmax1=ixOmax1-1;ixAmax2=ixOmax2-1;ixAmax3=ixOmax3-1;
      dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)= sum(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:),dim=ndim+1)
      do idims=1,ndim
        ixAmin1=ixOmin1-kr(idims,1);ixAmin2=ixOmin2-kr(idims,2)
        ixAmin3=ixOmin3-kr(idims,3);ixAmax1=ixOmax1-kr(idims,1)
        ixAmax2=ixOmax2-kr(idims,2);ixAmax3=ixOmax3-kr(idims,3);
        dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+block%surfaceC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           ixAmin3:ixAmax3,idims)
      end do
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=abs(divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*invB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    end if

  end subroutine get_normalized_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)  :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer, intent(out) :: idirmin
    integer :: idir, idirmin0

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)

    idirmin0 = 7-2*ndir

    call curlvector(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       mag(1:ndir)),ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,current,idirmin,idirmin0,ndir)

    if(B0field) current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin0:3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin0:3)+block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin0:3)
  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine mhd_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_cak_force, only: cak_get_dt

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx1,dx2,dx3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    integer                       :: idirmin,idim
    double precision              :: dxarr(ndim)
    double precision              :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    dtnew = bigdouble

    dxarr(1)=dx1;dxarr(2)=dx2;dxarr(3)=dx3;
    if (mhd_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/mhd_eta
    else if (mhd_eta<zero)then
       call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
       call usr_special_resistivity(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,x,&
          current,eta)
       dtnew=bigdouble
       do idim=1,ndim
         if(slab_uniform) then
           dtnew=min(dtnew,dtdiffpar/(smalldouble+maxval(eta(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)/dxarr(idim)**2)))
         else
           dtnew=min(dtnew,dtdiffpar/(smalldouble+maxval(eta(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)/block%ds(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)**2)))
         end if
       end do
    end if

    if(mhd_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/mhd_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1:ndim))**4/mhd_eta_hyper,dtnew)
      end if
    end if

    if(mhd_radiative_cooling) then
      call cooling_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x,&
         rc_fl)
    end if

    if(mhd_viscosity) then
      call viscosity_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    end if

    if(mhd_gravity) then
      call gravity_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    end if

    if(mhd_ambipolar_exp) then
      dtnew=min(dtdiffpar*get_ambipolar_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,&
         dx2,dx3,x),dtnew)
    endif

    if (mhd_cak_force) then
      call cak_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    end if

  end subroutine mhd_get_dt

  ! Add geometrical source terms to w
  subroutine mhd_add_source_geom(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_geometry
    use mod_rotating_frame, only: rotating_frame_add_source

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, dtfactor,x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    integer          :: iw,idir, h1xmin1,h1xmin2,h1xmin3,h1xmax1,h1xmax2,&
       h1xmax3, h2xmin1,h2xmin2,h2xmin3,h2xmax1,h2xmax2,h2xmax3
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    ! 1/rho
    invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.d0/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)
    ! include dt in invr, invr is always used with qdt
    if(local_timestep) then
      invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * dtfactor/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)
    else
      invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)
    endif  


    select case (coordinate)
    case (cylindrical)
      call mhd_get_p_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
      if(phi_>0) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)**2+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mphi_)**2*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)+invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mr_)*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,br_))
        if(.not.stagger_grid) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)+invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mr_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mphi_)) *invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        end if
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
      if(mhd_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    case (spherical)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmin3=ixOmin3-kr(1,3);h1xmax1=ixOmax1-kr(1,1)
       h1xmax2=ixOmax2-kr(1,2);h1xmax3=ixOmax3-kr(1,3)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmin3=ixOmin3-kr(2,3);h2xmax1=ixOmax1-kr(2,1)
       h2xmax2=ixOmax2-kr(2,2);h2xmax3=ixOmax3-kr(2,3);
       call mhd_get_p_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp1)
       ! m1
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)-block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,h1xmin3:h1xmax3,&
          1))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       do idir=2,ndir
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(idir))**2*invrho(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3)-wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))**2
       end do
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       ! b1
       if(mhd_glm) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))+invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*2.0d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,psi_)
       end if

       
       ! m2
       ! This will make hydrostatic p=const an exact solution
       if(local_timestep) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
       else
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = qdt * tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
       endif  
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,2)-block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          h2xmin3:h2xmax3,2)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(2))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mag(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mag(2)))
       if(ndir==3) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(3))**2*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(3))**2)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(1)))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
         if(mhd_glm) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) + dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,psi_)
         end if
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       end if
      

       if(ndir==3) then
         ! m3
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(1))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(1)))  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(3))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(3))) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(3))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
         ! b3
         if(.not.stagger_grid) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(1)))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2)) /(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,rho_)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2)))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)
         end if
       end if
    end select

    if (mhd_rotating_frame) then
       call rotating_frame_add_source(qdt,dtfactor,ixImin1,ixImin2,ixImin3,&
          ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
          ixOmax3,wCT,w,x)
    endif

  end subroutine mhd_add_source_geom

  ! Add geometrical source terms to w
  subroutine mhd_add_source_geom_split(qdt,dtfactor, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, dtfactor, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    integer          :: iw,idir, h1xmin1,h1xmin2,h1xmin3,h1xmax1,h1xmax2,&
       h1xmax3, h2xmin1,h2xmin2,h2xmin3,h2xmax1,h2xmax2,h2xmax3
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    if(has_equi_rho0) then
      invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_) + block%equi_vars(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho0_,b0i))
    else
      invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)
    end if
    ! include dt in invr, invr is always used with qdt
    if(local_timestep) then
      invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * dtfactor/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)
    else
      invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)
    endif  

    select case (coordinate)
    case (cylindrical)
      call mhd_get_p_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
      if(phi_>0) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)**2+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mphi_)**2*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)+qdt*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*(-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mr_)*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,br_))
        if(.not.stagger_grid) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)+invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mr_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mphi_)) *invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
        end if
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
      if(mhd_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    case (spherical)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmin3=ixOmin3-kr(1,3);h1xmax1=ixOmax1-kr(1,1)
       h1xmax2=ixOmax2-kr(1,2);h1xmax3=ixOmax3-kr(1,3)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmin3=ixOmin3-kr(2,3);h2xmax1=ixOmax1-kr(2,1)
       h2xmax2=ixOmax2-kr(2,2);h2xmax3=ixOmax3-kr(2,3);
       call mhd_get_p_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp1)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if(B0field) then
         tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=sum(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,:,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(:)),dim=ndim+1)
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       end if
       ! m1
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)-block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,h1xmin3:h1xmax3,&
          1))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(idir))**2*invrho(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)-wCT(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))**2
           if(B0field) tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)-2.0d0*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,idir,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(idir))
         end do
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(1))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       ! b1
       if(mhd_glm) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))+invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*2.0d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,psi_)
       end if

       
       ! m2
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if(B0field) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       end if
       if(local_timestep) then
         tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) = block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       else
         tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) = qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       endif  
       ! This will make hydrostatic p=const an exact solution
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))+tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,2)-block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          h2xmin3:h2xmax3,2)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom(2))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mag(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mag(2)))
       if (B0field) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,1,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(2)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mag(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,2,0)
       end if
       if(ndir==3) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(3))**2*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(3))**2)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)-2.0d0*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,3,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,mag(3))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2))
         end if
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(2))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(1)))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
         if(B0field) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2,0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,1,0))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)
         end if
         if(mhd_glm) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) + dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,psi_)
         end if
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       end if
      

       if(ndir==3) then
         ! m3
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(1))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(1)))  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(3))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(3))) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2))
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,1,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,mag(3)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,mag(1))*block%B0(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,&
               0)  +(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               2,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(3)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,3,0)) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,2))
         end if
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom(3))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
         ! b3
         if(.not.stagger_grid) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(1)))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2)) *invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))
           if (B0field) then
              tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mom(1))*block%B0(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,0) -wCT(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(3))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,1,0))*invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mom(3))*block%B0(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,0) -wCT(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,3,0))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,2)) *invrho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,2))
           end if
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)*invr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)
         end if
       end if
    end select
  end subroutine mhd_add_source_geom_split

  !> Compute 2 times total magnetic energy
  function mhd_mag_en_all(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if (B0field) then
      mge = sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,&
         b0i))**2, dim=ndim+1)
    else
      mge = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(:))**2,&
          dim=ndim+1)
    end if
  end function mhd_mag_en_all

  !> Compute full magnetic field by direction
  function mhd_mag_i_all(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: mgf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if (B0field) then
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mag(idir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir,b0i)
    else
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(idir))
    end if
  end function mhd_mag_i_all

  !> Compute evolving magnetic energy
  function mhd_mag_en(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    mge = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        mag(:))**2, dim=ndim+1)
  end function mhd_mag_en

  !> compute kinetic energy
  function mhd_kin_en_origin(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    if (present(inv_rho)) then
      ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(:))**2, dim=ndim+1) * inv_rho
    else
      if(has_equi_rho0) then
        ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(:))**2, dim=ndim+1) / (w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            rho_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_rho0_,0))
      else
        ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, mom(:))**2, dim=ndim+1) / w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3, rho_)
      end if
    end if
  end function mhd_kin_en_origin

  !> compute kinetic energy
  function mhd_kin_en_boris(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, inv_rho) result(ke)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    if (present(inv_rho)) then
      ke=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)*inv_rho*inv_squared_c)
      ke=0.5d0*sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(:)))**2,dim=ndim+1)*ke**2*inv_rho
    else
      ke=1.d0/(1.d0+sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))**2,dim=ndim+1)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_)*inv_squared_c)
      ke=0.5d0*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(:))**2,dim=ndim+1)*ke**2/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, rho_)
    end if
  end function mhd_kin_en_boris

  subroutine mhd_getv_Hall(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: vHall(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:3)

    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)
    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call mhd_get_rho(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3) = zero
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3) = - mhd_etah*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,idirmin:3)
    do idir = idirmin, 3
       vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir) = vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

  end subroutine mhd_getv_Hall

  subroutine mhd_get_Jambi(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, allocatable, intent(inout) :: res(:,:,:,:)


    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)

    res = 0d0

    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
 
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3)=-current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3)
    do idir = idirmin, 3
      call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,res(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,idir),w,x)
    enddo

  end subroutine mhd_get_Jambi

    ! COMMENTED  because we have that in cmax now:
!  subroutine mhd_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
!    use mod_global_parameters
!
!    integer, intent(in) :: ixI^L, ixO^L
!    double precision, intent(in)    :: dx^D
!    double precision, intent(in)    :: w(ixI^S,1:nw)
!    double precision, intent(in)    :: x(ixI^S,1:ndim)
!    double precision, intent(out)   :: dthall
!    !.. local ..
!    double precision :: dxarr(ndim)
!    double precision :: bmag(ixI^S)
!
!    dthall=bigdouble
!
!
!    ^D&dxarr(^D)=dx^D;
!
!    if (.not. B0field) then
!       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
!    else
!       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,b0i))**2))
!    end if
!
!    if(slab_uniform) then
!      dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(mhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    else
!      dthall=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2.0d0/(mhd_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_)))
!    end if
!
!  end subroutine mhd_getdt_Hall

  subroutine mhd_modify_wLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,wLC,wRC,wLp,wRp,s,&
     idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    type(state)                     :: s
    double precision                :: dB(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), dPsi(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    if(stagger_grid) then
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=s%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=s%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=s%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=s%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)
    else
      ! Solve the Riemann problem for the linear 2x2 system for normal
      ! B-field and GLM_Psi according to Dedner 2002:
      ! This implements eq. (42) in Dedner et al. 2002 JcP 175
      ! Gives the Riemann solution on the interface
      ! for the normal B component and Psi in the GLM-MHD system.
      ! 23/04/2013 Oliver Porth
      dB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)   = wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idir)) - wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idir))
      dPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_) - wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_)

      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))   = 0.5d0 * (wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idir)) + wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idir))) - 0.5d0/cmax_global * &
         dPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)       = 0.5d0 * (wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_) + wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_)) - 0.5d0*cmax_global * dB(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)

      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir)) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))
      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)

      if(total_energy) then
        wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)-half*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
        wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)-half*wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
      end if
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir)) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir)) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)
      ! modify total energy according to the change of magnetic field
      if(total_energy) then
        wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+half*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
        wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+half*wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
      end if
    end if

    if(associated(usr_set_wLR)) call usr_set_wLR(ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine mhd_modify_wLR

  subroutine mhd_boundary_adjust(igrid,psb)
    use mod_global_parameters
    integer, intent(in) :: igrid
    type(state), target :: psb(max_blocks)

    integer :: iB, idims, iside, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, i1,i2,i3

    block=>ps(igrid)
    dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
    dxlevel(3)=rnode(rpdx3_,igrid);
    do idims=1,ndim
       ! to avoid using as yet unknown corner info in more than 1D, we
       ! fill only interior mesh ranges of the ghost cell ranges at first,
       ! and progressively enlarge the ranges to include corners later
       do iside=1,2
          i1=kr(1,idims)*(2*iside-3);i2=kr(2,idims)*(2*iside-3)
          i3=kr(3,idims)*(2*iside-3);
          if (neighbor_type(i1,i2,i3,igrid)/=1) cycle
          iB=(idims-1)*2+iside
          if(.not.boundary_divbfix(iB)) cycle
          if(any(typeboundary(:,iB)==bc_special)) then
            ! MF nonlinear force-free B field extrapolation and data driven
            ! require normal B of the first ghost cell layer to be untouched by
            ! fixdivB=0 process, set boundary_divbfix_skip(iB)=1 in par file
            select case (idims)
            case (1)
               if (iside==2) then
                  ! maximal boundary
                  ixOmin1=ixGhi1+1-nghostcells+boundary_divbfix_skip(2*1)
                  ixOmin2=ixGlo2;ixOmin3=ixGlo3;
                  ixOmax1=ixGhi1;ixOmax2=ixGhi2;ixOmax3=ixGhi3;
               else
                  ! minimal boundary
                  ixOmin1=ixGlo1;ixOmin2=ixGlo2;ixOmin3=ixGlo3;
                  ixOmax1=ixGlo1-1+nghostcells-boundary_divbfix_skip(2*1-1)
                  ixOmax2=ixGhi2;ixOmax3=ixGhi3;
               end if 
            case (2)
               if (iside==2) then
                  ! maximal boundary
                  ixOmin1=ixGlo1
                  ixOmin2=ixGhi2+1-nghostcells+boundary_divbfix_skip(2*2)
                  ixOmin3=ixGlo3;
                  ixOmax1=ixGhi1;ixOmax2=ixGhi2;ixOmax3=ixGhi3;
               else
                  ! minimal boundary
                  ixOmin1=ixGlo1;ixOmin2=ixGlo2;ixOmin3=ixGlo3;
                  ixOmax1=ixGhi1
                  ixOmax2=ixGlo2-1+nghostcells-boundary_divbfix_skip(2*2-1)
                  ixOmax3=ixGhi3;
               end if 
            case (3)
               if (iside==2) then
                  ! maximal boundary
                  ixOmin1=ixGlo1;ixOmin2=ixGlo2
                  ixOmin3=ixGhi3+1-nghostcells+boundary_divbfix_skip(2*3);
                  ixOmax1=ixGhi1;ixOmax2=ixGhi2;ixOmax3=ixGhi3;
               else
                  ! minimal boundary
                  ixOmin1=ixGlo1;ixOmin2=ixGlo2;ixOmin3=ixGlo3;
                  ixOmax1=ixGhi1;ixOmax2=ixGhi2
                  ixOmax3=ixGlo3-1+nghostcells-boundary_divbfix_skip(2*3-1);
               end if 
            end select
            call fixdivB_boundary(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
               ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,psb(igrid)%w,&
               psb(igrid)%x,iB)
          end if
       end do
    end do

  end subroutine mhd_boundary_adjust

  subroutine fixdivB_boundary(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,iB)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)

    double precision :: dx1x2,dx1x3,dx2x1,dx2x3,dx3x1,dx3x2
    integer :: ix1,ix2,ix3,ixFmin1,ixFmin2,ixFmin3,ixFmax1,ixFmax2,ixFmax3

    select case(iB)
     case(1)
       ! 2nd order CD for divB=0 to set normal B component better
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
       
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=w(ix1+1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) +dx1x2*(w(ix1,&
              ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) +dx1x3*(w(ix1,&
              ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-w(ix1,&
              ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=( (w(ix1+1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+w(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,mag(1)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,1)+(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,&
              mag(2))+w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              2)-(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+w(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*block%surfaceC(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,2)+(w(ix1,ixFmin2:ixFmax2,&
              ixFmin3+1:ixFmax3+1,mag(3))+w(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,mag(3)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,3)-(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              mag(3))+w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,&
              mag(3)))*block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,&
              3) )/block%surfaceC(ix1-1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              1)-w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
      
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case(2)
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
       
       
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         dx1x3=dxlevel(1)/dxlevel(3)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=w(ix1-1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1)) -dx1x2*(w(ix1,&
              ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2))) -dx1x3*(w(ix1,&
              ixFmin2:ixFmax2,ixFmin3+1:ixFmax3+1,mag(3))-w(ix1,&
              ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))=( (w(ix1-1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))+w(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,mag(1)))*block%surfaceC(ix1-1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,1)-(w(ix1,ixFmin2+1:ixFmax2+1,ixFmin3:ixFmax3,&
              mag(2))+w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              2)+(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(2))+w(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,mag(2)))*block%surfaceC(ix1,&
              ixFmin2-1:ixFmax2-1,ixFmin3:ixFmax3,2)-(w(ix1,ixFmin2:ixFmax2,&
              ixFmin3+1:ixFmax3+1,mag(3))+w(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,mag(3)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              ixFmin3:ixFmax3,3)+(w(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,&
              mag(3))+w(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,&
              mag(3)))*block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3-1:ixFmax3-1,&
              3) )/block%surfaceC(ix1,ixFmin2:ixFmax2,ixFmin3:ixFmax3,1)-w(ix1,&
              ixFmin2:ixFmax2,ixFmin3:ixFmax3,mag(1))
         end do
       end if
      
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case(3)
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
       
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
              ix2+1,ixFmin3:ixFmax3,mag(2)) +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              ixFmin3:ixFmax3,mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,&
              ixFmin3:ixFmax3,mag(1))) +dx2x3*(w(ixFmin1:ixFmax1,ix2,&
              ixFmin3+1:ixFmax3+1,mag(3))-w(ixFmin1:ixFmax1,ix2,&
              ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,&
              mag(2))=( (w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              2)+(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,&
              mag(1))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              1)-(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,&
              1)+(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,&
              mag(3))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              3)-(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(3))+w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,&
              3) )/block%surfaceC(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,&
              2)-w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
      
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case(4)
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
       
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         dx2x3=dxlevel(2)/dxlevel(3)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,mag(2))=w(ixFmin1:ixFmax1,&
              ix2-1,ixFmin3:ixFmax3,mag(2)) -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              ixFmin3:ixFmax3,mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,&
              ixFmin3:ixFmax3,mag(1))) -dx2x3*(w(ixFmin1:ixFmax1,ix2,&
              ixFmin3+1:ixFmax3+1,mag(3))-w(ixFmin1:ixFmax1,ix2,&
              ixFmin3-1:ixFmax3-1,mag(3)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,ixFmin3:ixFmax3,&
              mag(2))=( (w(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2-1,ixFmin3:ixFmax3,&
              2)-(w(ixFmin1+1:ixFmax1+1,ix2,ixFmin3:ixFmax3,&
              mag(1))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              1)+(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,ixFmin3:ixFmax3,&
              1)-(w(ixFmin1:ixFmax1,ix2,ixFmin3+1:ixFmax3+1,&
              mag(3))+w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              3)+(w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              mag(3))+w(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3-1:ixFmax3-1,&
              3) )/block%surfaceC(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,&
              2)-w(ixFmin1:ixFmax1,ix2,ixFmin3:ixFmax3,mag(2))
         end do
       end if
      
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     
     case(5)
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3+1
       ixFmax3=ixOmax3+1
       if(slab_uniform) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmax3,ixFmin3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,mag(3))=w(ixFmin1:ixFmax1,&
              ixFmin2:ixFmax2,ix3+1,mag(3)) +dx3x1*(w(ixFmin1+1:ixFmax1+1,&
              ixFmin2:ixFmax2,ix3,mag(1))-w(ixFmin1-1:ixFmax1-1,&
              ixFmin2:ixFmax2,ix3,mag(1))) +dx3x2*(w(ixFmin1:ixFmax1,&
              ixFmin2+1:ixFmax2+1,ix3,mag(2))-w(ixFmin1:ixFmax1,&
              ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmax3,ixFmin3,-1
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,&
              mag(3))=( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,&
              mag(3))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              3)+(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,&
              mag(1))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              1)-(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(1))+w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,&
              1)+(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,&
              mag(2))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              2)-(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(2))+w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,&
              2) )/block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,&
              3)-w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
     case(6)
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       ixFmin3=ixOmin3-1
       ixFmax3=ixOmax3-1
       if(slab_uniform) then
         dx3x1=dxlevel(3)/dxlevel(1)
         dx3x2=dxlevel(3)/dxlevel(2)
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,mag(3))=w(ixFmin1:ixFmax1,&
              ixFmin2:ixFmax2,ix3-1,mag(3)) -dx3x1*(w(ixFmin1+1:ixFmax1+1,&
              ixFmin2:ixFmax2,ix3,mag(1))-w(ixFmin1-1:ixFmax1-1,&
              ixFmin2:ixFmax2,ix3,mag(1))) -dx3x2*(w(ixFmin1:ixFmax1,&
              ixFmin2+1:ixFmax2+1,ix3,mag(2))-w(ixFmin1:ixFmax1,&
              ixFmin2-1:ixFmax2-1,ix3,mag(2)))
         end do
       else
         do ix3=ixFmin3,ixFmax3
           w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3+1,&
              mag(3))=( (w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,&
              mag(3))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(3)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3-1,&
              3)-(w(ixFmin1+1:ixFmax1+1,ixFmin2:ixFmax2,ix3,&
              mag(1))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              1)+(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(1))+w(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ixFmin2:ixFmax2,ix3,&
              1)-(w(ixFmin1:ixFmax1,ixFmin2+1:ixFmax2+1,ix3,&
              mag(2))+w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              2)+(w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              mag(2))+w(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ixFmin2-1:ixFmax2-1,ix3,&
              2) )/block%surfaceC(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,&
              3)-w(ixFmin1:ixFmax1,ixFmin2:ixFmax2,ix3,mag(3))
         end do
       end if
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
          ixGmax2,ixGmax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  
  subroutine mhd_clean_divb_multigrid(qdt, qt, active)
    use mod_forest
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_geometry

    double precision, intent(in) :: qdt    !< Current time step
    double precision, intent(in) :: qt     !< Current time
    logical, intent(inout)       :: active !< Output if the source is active
    integer                      :: iigrid, igrid, id
    integer                      :: n, nc, lvl, ixmin1,ixmin2,ixmin3,ixmax1,&
       ixmax2,ixmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idim
    type(tree_node), pointer     :: pnode
    double precision             :: tmp(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
       ixGlo3:ixGhi3), grad(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3, ndim)
    double precision             :: res
    double precision, parameter  :: max_residual = 1d-3
    double precision, parameter  :: residual_reduction = 1d-10
    integer, parameter           :: max_its      = 50
    double precision             :: residual_it(max_its), max_divb

    mg%operator_type = mg_laplacian

    ! Set boundary conditions
    do n = 1, 2*ndim
       idim = (n+1)/2
       select case (typeboundary(mag(idim), n))
       case (bc_symm)
          ! d/dx B = 0, take phi = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_asymm)
          ! B = 0, so grad(phi) = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_cont)
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_special)
          ! Assume Dirichlet boundary conditions, derivative zero
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case (bc_periodic)
          ! Nothing to do here
       case default
          write(*,*) "mhd_clean_divb_multigrid warning: unknown boundary type"
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       end select
    end do

    ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmin3=ixMlo3-1;ixmax1=ixMhi1+1
    ixmax2=ixMhi2+1;ixmax3=ixMhi3+1;
    max_divb = 0.0d0

    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);

       call get_divb(ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
           1:nw), ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, ixMlo1,ixMlo2,&
          ixMlo3,ixMhi1,ixMhi2,ixMhi3, tmp, mhd_divb_4thorder)
       mg%boxes(id)%cc(1:nc,1:nc,1:nc, mg_irhs) = tmp(ixMlo1:ixMhi1,&
          ixMlo2:ixMhi2,ixMlo3:ixMhi3)
       max_divb = max(max_divb, maxval(abs(tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
          ixMlo3:ixMhi3))))
    end do

    ! Solve laplacian(phi) = divB
    if(stagger_grid) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, max_divb, 1, MPI_DOUBLE_PRECISION,&
          MPI_MAX, icomm, ierrmpi)

      if (mype == 0) print *, "Performing multigrid divB cleaning"
      if (mype == 0) print *, "iteration vs residual"
      ! Solve laplacian(phi) = divB
      do n = 1, max_its
         call mg_fas_fmg(mg, n>1, max_res=residual_it(n))
         if (mype == 0) write(*, "(I4,E11.3)") n, residual_it(n)
         if (residual_it(n) < residual_reduction * max_divb) exit
      end do
      if (mype == 0 .and. n > max_its) then
         print *, "divb_multigrid warning: not fully converged"
         print *, "current amplitude of divb: ", residual_it(max_its)
         print *, "multigrid smallest grid: ", &
              mg%domain_size_lvl(:, mg%lowest_lvl)
         print *, "note: smallest grid ideally has <= 8 cells"
         print *, "multigrid dx/dy/dz ratio: ", mg%dr(:, 1)/mg%dr(1, 1)
         print *, "note: dx/dy/dz should be similar"
      end if
    else
      do n = 1, max_its
         call mg_fas_vcycle(mg, max_res=res)
         if (res < max_residual) exit
      end do
      if (res > max_residual) call mpistop("divb_multigrid: no convergence")
    end if


    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);

       ! Compute the gradient of phi
       tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = mg%boxes(id)%cc(:,:,:,&
           mg_iphi)

       if(stagger_grid) then
         do idim =1, ndim
           ixCmin1=ixMlo1-kr(idim,1);ixCmin2=ixMlo2-kr(idim,2)
           ixCmin3=ixMlo3-kr(idim,3);
           ixCmax1=ixMhi1;ixCmax2=ixMhi2;ixCmax3=ixMhi3;
           call gradientx(tmp,ps(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
              ixGhi3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim,&
              grad(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,idim),.false.)
           ! Apply the correction B* = B - gradient(phi)
           ps(igrid)%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
              idim)=ps(igrid)%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3,idim)-grad(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              ixCmin3:ixCmax3,idim)
         end do
         ! store cell-center magnetic energy
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3) = sum(ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3, mag(1:ndim))**2, dim=ndim+1)
         ! change cell-center magnetic field
         call mhd_face_to_center(ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
            ps(igrid))
       else
         do idim = 1, ndim
            call gradient(tmp,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
               ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,idim,grad(ixGlo1:ixGhi1,&
               ixGlo2:ixGhi2,ixGlo3:ixGhi3, idim))
         end do
         ! store cell-center magnetic energy
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3) = sum(ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3, mag(1:ndim))**2, dim=ndim+1)
         ! Apply the correction B* = B - gradient(phi)
         ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
             mag(1:ndim)) = ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3, mag(1:ndim)) - grad(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3, :)
       end if

       if(total_energy) then
         ! Determine magnetic energy difference
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3) = 0.5_dp * (sum(ps(igrid)%w(ixMlo1:ixMhi1,&
            ixMlo2:ixMhi2,ixMlo3:ixMhi3, mag(1:ndim))**2,&
             dim=ndim+1) - tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))
         ! Keep thermal pressure the same
         ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
             e_) = ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
             e_) + tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)
       end if
    end do

    active = .true.

  end subroutine mhd_clean_divb_multigrid
 

  subroutine mhd_update_faces(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,wprim,fC,fE,sCT,s,&
     vcts)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)

    select case(type_ct)
    case('average')
      call update_faces_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fC,fE,&
         sCT,s)
    case('uct_contact')
      call update_faces_contact(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,wprim,&
         fC,fE,sCT,s,vcts)
    case('uct_hll')
      call update_faces_hll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fE,sCT,s,vcts)
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine mhd_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fC,fE,sCT,&
     s)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)

    integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,&
       jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3,ixCmmin1,ixCmmin2,ixCmmin3,&
       ixCmmax1,ixCmmax2,ixCmmax3
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2
    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3) :: E_resi, E_ambi

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(mhd_ambipolar_exp) call get_ambipolar_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT%w,x,E_ambi)

    do idim1=1,ndim
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=sdim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
            ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
            ixCmin3=ixOmin3+kr(idir,3)-1;
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
            jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
            jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmin3=ixCmin3+kr(idim2,3);hxCmax1=ixCmax1+kr(idim2,1)
            hxCmax2=ixCmax2+kr(idim2,2);hxCmax3=ixCmax3+kr(idim2,3);
            ! Interpolate to edges
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=quarter*(fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,iwdim1,idim2)+fC(jxCmin1:jxCmax1,&
               jxCmin2:jxCmax2,jxCmin3:jxCmax3,iwdim1,&
               idim2)-fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               iwdim2,idim1)-fC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
               hxCmin3:hxCmax3,iwdim2,idim1))

            ! add resistive electric field at cell edges E=-vxB+eta J
            if(mhd_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_resi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)
            ! add ambipolar electric field
            if(mhd_ambipolar_exp) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_ambi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)

            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=qdt*s%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)*fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)

            if (.not.slab) then
              where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                 r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
      do idim2=1,ndim
        do idir=sdim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
            hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
            hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
            ! Add line integrals in direction idir
            circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) > 1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)=bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
    end do

    end associate

  end subroutine update_faces_average

  !> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
  subroutine update_faces_contact(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,wp,fC,fE,&
     sCT,s,vcts)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)

    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! electric field at cell centers
    double precision                   :: ECC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)
    double precision                   :: Ein(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)
    ! gradient of E at left and right side of a cell face
    double precision                   :: EL(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),ER(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    ! gradient of E at left and right side of a cell corner
    double precision                   :: ELC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),ERC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3) :: E_resi, E_ambi
    ! total magnetic field at cell centers
    double precision                   :: Btot(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! current on cell edges
    double precision :: jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       sdim:3)
    ! location at cell faces
    double precision :: xs(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3,&
       1:ndim)
    double precision :: gradi(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3)
    integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,&
       jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,&
       ixAmax2,ixAmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2,ix1,&
       ix2,ix3

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm,wCTs=>sCT%ws)

    if(B0field) then
      Btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)=wp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndim))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim,0)
    else
      Btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndim)=wp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         mag(1:ndim))
    end if
    ECC=0.d0
    ! Calculate electric field at cell centers
    do idim1=1,ndim; do idim2=1,ndim; do idir=sdim,3
      if(lvc(idim1,idim2,idir)==1)then
         ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)=ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)+Btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idim1)*wp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            mom(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)=ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)-Btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idim1)*wp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            mom(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(mhd_ambipolar_exp) call get_ambipolar_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT%w,x,E_ambi)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ! evaluate electric field along cell edges according to equation (41)
    do idim1=1,ndim
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=sdim,3 ! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
            ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
            ixCmin3=ixOmin3+kr(idir,3)-1;
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
            jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
            jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmin3=ixCmin3+kr(idim2,3);hxCmax1=ixCmax1+kr(idim2,1)
            hxCmax2=ixCmax2+kr(idim2,2);hxCmax3=ixCmax3+kr(idim2,3);
            ! average cell-face electric field to cell edges
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=quarter*(fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,iwdim1,idim2)+fC(jxCmin1:jxCmax1,&
               jxCmin2:jxCmax2,jxCmin3:jxCmax3,iwdim1,&
               idim2)-fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               iwdim2,idim1)-fC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
               hxCmin3:hxCmax3,iwdim2,idim1))
            if(partial_energy) Ein(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)

            ! add slope in idim2 direction from equation (50)
            ixAmin1=ixCmin1;ixAmin2=ixCmin2;ixAmin3=ixCmin3;
            ixAmax1=ixCmax1+kr(idim1,1);ixAmax2=ixCmax2+kr(idim1,2)
            ixAmax3=ixCmax3+kr(idim1,3);
            EL(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=fC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,iwdim1,idim2)-ECC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,ixAmin3:ixAmax3,idir)
            hxCmin1=ixAmin1+kr(idim2,1);hxCmin2=ixAmin2+kr(idim2,2)
            hxCmin3=ixAmin3+kr(idim2,3);hxCmax1=ixAmax1+kr(idim2,1)
            hxCmax2=ixAmax2+kr(idim2,2);hxCmax3=ixAmax3+kr(idim2,3);
            ER(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=fC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,iwdim1,idim2)-ECC(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir)
            where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)>0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=EL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)
            else where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)<0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3)
            else where
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=0.5d0*(EL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)+EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3))
            end where
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmin3=ixCmin3+kr(idim2,3);hxCmax1=ixCmax1+kr(idim2,1)
            hxCmax2=ixCmax2+kr(idim2,2);hxCmax3=ixCmax3+kr(idim2,3);
            where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
               idim1)>0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=ER(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)
            else where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
               idim1)<0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3)
            else where
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=0.5d0*(ER(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)+ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3))
            end where
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)+0.25d0*(ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)+ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3))

            ! add slope in idim1 direction from equation (50)
            jxCmin1=ixCmin1+kr(idim2,1);jxCmin2=ixCmin2+kr(idim2,2)
            jxCmin3=ixCmin3+kr(idim2,3);jxCmax1=ixCmax1+kr(idim2,1)
            jxCmax2=ixCmax2+kr(idim2,2);jxCmax3=ixCmax3+kr(idim2,3);
            ixAmin1=ixCmin1;ixAmin2=ixCmin2;ixAmin3=ixCmin3;
            ixAmax1=ixCmax1+kr(idim2,1);ixAmax2=ixCmax2+kr(idim2,2)
            ixAmax3=ixCmax3+kr(idim2,3);
            EL(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=-fC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,iwdim2,idim1)-ECC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,ixAmin3:ixAmax3,idir)
            hxCmin1=ixAmin1+kr(idim1,1);hxCmin2=ixAmin2+kr(idim1,2)
            hxCmin3=ixAmin3+kr(idim1,3);hxCmax1=ixAmax1+kr(idim1,1)
            hxCmax2=ixAmax2+kr(idim1,2);hxCmax3=ixAmax3+kr(idim1,3);
            ER(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3)=-fC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
               ixAmin3:ixAmax3,iwdim2,idim1)-ECC(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir)
            where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim2)>0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=EL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)
            else where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim2)<0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3)
            else where
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=0.5d0*(EL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)+EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3))
            end where
            hxCmin1=ixCmin1+kr(idim1,1);hxCmin2=ixCmin2+kr(idim1,2)
            hxCmin3=ixCmin3+kr(idim1,3);hxCmax1=ixCmax1+kr(idim1,1)
            hxCmax2=ixCmax2+kr(idim1,2);hxCmax3=ixCmax3+kr(idim1,3);
            where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
               idim2)>0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=ER(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)
            else where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,hxCmin3:hxCmax3,&
               idim2)<0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3)
            else where
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)=0.5d0*(ER(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3)+ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
                 jxCmin3:jxCmax3))
            end where
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)+0.25d0*(ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)+ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3))
            ! difference between average and upwind interpolated E
            if(partial_energy) Ein(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)-Ein(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)
            ! add resistive electric field at cell edges E=-vxB+eta J
            if(mhd_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_resi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)
            ! add ambipolar electric field
            if(mhd_ambipolar_exp) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_ambi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)

            ! times time step and edge length
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)*qdt*s%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)
          end if
        end do
      end do
    end do

    if(partial_energy) then
      ! add upwind diffused magnetic energy back to energy
      ! calculate current density at cell edges
      jce=0.d0
      do idim1=1,ndim
        do idim2=1,ndim
          do idir=sdim,3
            if (lvc(idim1,idim2,idir)==0) cycle
            ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
            ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
            ixCmin3=ixOmin3+kr(idir,3)-1;
            ixBmax1=ixCmax1-kr(idir,1)+1;ixBmax2=ixCmax2-kr(idir,2)+1
            ixBmax3=ixCmax3-kr(idir,3)+1;
            ixBmin1=ixCmin1;ixBmin2=ixCmin2;ixBmin3=ixCmin3;
            ! current at transverse faces
            xs(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
               :)=x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,:)
            xs(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
               idim2)=x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
               idim2)+half*s%dx(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
               ixBmin3:ixBmax3,idim2)
            call gradientx(wCTs(ixGslo1:ixGshi1,ixGslo2:ixGshi2,&
               ixGslo3:ixGshi3,idim2),xs,ixGslo1,ixGslo2,ixGslo3,ixGshi1,&
               ixGshi2,ixGshi3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
               idim1,gradi,.false.)
            if (lvc(idim1,idim2,idir)==1) then
              jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                 idir)=jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                 idir)+gradi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
            else
              jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                 idir)=jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                 idir)-gradi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
            end if
          end do
        end do
      end do
      if(nwextra>0) block%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         nw)=0.d0
      do idir=sdim,3
        ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
        ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
        ixCmin3=ixOmin3+kr(idir,3)-1;
        ! E dot J on cell edges
        Ein(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)=Ein(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)*jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)
        ! average from cell edge to cell center
        jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idir)=0.d0
        do ix3=-1,0
        do ix2=-1,0
        do ix1=-1,0
           if( ix1==-1 .and. 1==idir  .or. ix2==-1 .and. 2==idir  .or. ix3==-1 &
              .and. 3==idir ) cycle
           ixAmin1=ixOmin1+ix1;ixAmin2=ixOmin2+ix2;ixAmin3=ixOmin3+ix3;
           ixAmax1=ixOmax1+ix1;ixAmax2=ixOmax2+ix2;ixAmax3=ixOmax3+ix3;
           jce(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              idir)=jce(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              idir)+Ein(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,idir)
        end do
        end do
        end do
        where(jce(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)<0.d0)
          jce(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)=0.d0
        end where
        ! save additional numerical resistive heating to an extra variable
        if(nwextra>0) block%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           nw)=block%w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           nw)+0.25d0*jce(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+qdt*0.25d0*jce(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)
      end do
    end if

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
      do idim2=1,ndim
        do idir=sdim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
            hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
            hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
            ! Add line integrals in direction idir
            circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) > smalldouble)
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)=bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
    end do

    end associate

  end subroutine update_faces_contact

  !> update faces
  subroutine update_faces_hll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,fE,sCT,s,vcts)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: qt, qdt
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)
    type(state)                        :: sCT, s
    type(ct_velocity)                  :: vcts

    double precision                   :: vtilL(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,2)
    double precision                   :: vtilR(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,2)
    double precision                   :: bfacetot(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    double precision                   :: btilL(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    double precision                   :: btilR(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndim)
    double precision                   :: cp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,2)
    double precision                   :: cm(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,2)
    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    ! non-ideal electric field on cell edges
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3) :: E_resi, E_ambi
    integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
       ixCpmin1,ixCpmin2,ixCpmin3,ixCpmax1,ixCpmax2,ixCpmax3,jxCmin1,jxCmin2,&
       jxCmin3,jxCmax1,jxCmax2,jxCmax3,ixCmmin1,ixCmmin2,ixCmmin3,ixCmmax1,&
       ixCmmax2,ixCmmax3
    integer                            :: idim1,idim2,idir

    associate(bfaces=>s%ws,bfacesCT=>sCT%ws,x=>s%x,vbarC=>vcts%vbarC,&
       cbarmin=>vcts%cbarmin,cbarmax=>vcts%cbarmax)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! Loop over components of electric field

    ! idir: electric field component we need to calculate
    ! idim1: directions in which we already performed the reconstruction
    ! idim2: directions in which we perform the reconstruction

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT,s,E_resi)

    ! if there is ambipolar diffusion, get E_ambi
    if(mhd_ambipolar_exp) call get_ambipolar_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT%w,x,E_ambi)

    do idir=sdim,3
      ! Indices
      ! idir: electric field component
      ! idim1: one surface
      ! idim2: the other surface
      ! cyclic permutation: idim1,idim2,idir=1,2,3
      ! Velocity components on the surface
      ! follow cyclic premutations:
      ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-1+kr(idir,1);ixCmin2=ixOmin2-1+kr(idir,2)
      ixCmin3=ixOmin3-1+kr(idir,3);

      ! Set indices and directions
      idim1=mod(idir,3)+1
      idim2=mod(idir+1,3)+1

      jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
      jxCmin3=ixCmin3+kr(idim1,3);jxCmax1=ixCmax1+kr(idim1,1)
      jxCmax2=ixCmax2+kr(idim1,2);jxCmax3=ixCmax3+kr(idim1,3);
      ixCpmin1=ixCmin1+kr(idim2,1);ixCpmin2=ixCmin2+kr(idim2,2)
      ixCpmin3=ixCmin3+kr(idim2,3);ixCpmax1=ixCmax1+kr(idim2,1)
      ixCpmax2=ixCmax2+kr(idim2,2);ixCpmax3=ixCmax3+kr(idim2,3);

      ! Reconstruct transverse transport velocities
      call reconstruct(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,vbarC(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,idim1,1),vtilL(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2),vtilR(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2))

      call reconstruct(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,vbarC(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,idim2,2),vtilL(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1),vtilR(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1))

      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      if(B0field) then
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1)+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1,idim1)
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2)+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2,idim2)
      else
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim1)
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idim2)
      end if
      call reconstruct(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim2,&
         bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
         btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1),&
         btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim1))

      call reconstruct(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
         ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,&
         bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
         btilL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2),&
         btilR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idim2))

      ! Take the maximum characteristic

      cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(cbarmin(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,&
         idim1),cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1))
      cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(cbarmax(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,&
         idim1),cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1))

      cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)=max(cbarmin(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
         cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2))
      cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)=max(cbarmax(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,idim2),&
         cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim2))
     

      ! Calculate eletric field
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=-(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim2) + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim2) - cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)*(btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim2)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim2)))/(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)+cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)) +(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) - cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)*(btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)))/(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)+cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2))

      ! add resistive electric field at cell edges E=-vxB+eta J
      if(mhd_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)+E_resi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)
      ! add ambipolar electric field
      if(mhd_ambipolar_exp) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)+E_ambi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)

      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=qdt*s%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)*fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)

      if (.not.slab) then
        where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           r_)+half*dxlevel(r_)).lt.1.0d-9)
          fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=zero
        end where
      end if

    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)=zero

    ! Calculate circulation on each face: interal(fE dot dl)
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
      do idim2=1,ndim
        do idir=sdim,3 ! Direction of line integral
          ! Assemble indices
          if(lvc(idim1,idim2,idir)/=0) then
            hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
            hxCmin3=ixCmin3-kr(idim2,3);hxCmax1=ixCmax1-kr(idim2,1)
            hxCmax2=ixCmax2-kr(idim2,2);hxCmax3=ixCmax3-kr(idim2,3);
            ! Add line integrals in direction idir
            circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)-fE(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,hxCmin3:hxCmax3,idir))
          end if
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1) > 1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3))
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)=circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)=bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim1)
    end do

    end associate
  end subroutine update_faces_hll

  !> calculate eta J at cell edges
  subroutine get_resistive_electric_field(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,sCT,s,&
     jce)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    type(state), intent(in)            :: sCT, s
    ! current on cell edges
    double precision :: jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       sdim:3)

    ! current on cell centers
    double precision :: jcc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       7-2*ndir:3)
    ! location at cell faces
    double precision :: xs(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3,&
       1:ndim)
    ! resistivity
    double precision :: eta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: gradi(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3)
    integer :: ix1,ix2,ix3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,ixBmin1,ixBmin2,ixBmin3,&
       ixBmax1,ixBmax2,ixBmax3,idir,idirmin,idim1,idim2

    associate(x=>s%x,dx=>s%dx,w=>s%w,wCT=>sCT%w,wCTs=>sCT%ws)
    ! calculate current density at cell edges
    jce=0.d0
    do idim1=1,ndim
      do idim2=1,ndim
        do idir=sdim,3
          if (lvc(idim1,idim2,idir)==0) cycle
          ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
          ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
          ixCmin3=ixOmin3+kr(idir,3)-1;
          ixBmax1=ixCmax1-kr(idir,1)+1;ixBmax2=ixCmax2-kr(idir,2)+1
          ixBmax3=ixCmax3-kr(idir,3)+1;
          ixBmin1=ixCmin1;ixBmin2=ixCmin2;ixBmin3=ixCmin3;
          ! current at transverse faces
          xs(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
             :)=x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,:)
          xs(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
             idim2)=x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
             idim2)+half*dx(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
             idim2)
          call gradientx(wCTs(ixGslo1:ixGshi1,ixGslo2:ixGshi2,ixGslo3:ixGshi3,&
             idim2),xs,ixGslo1,ixGslo2,ixGslo3,ixGshi1,ixGshi2,ixGshi3,ixCmin1,&
             ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim1,gradi,.true.)
          if (lvc(idim1,idim2,idir)==1) then
            jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)+gradi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
          else
            jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)-gradi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
          end if
        end do
      end do
    end do
    ! get resistivity
    if(mhd_eta>zero)then
      jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         :)=jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:)*mhd_eta
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
      ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
      call get_current(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,jcc)
      call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,idirmin,x,jcc,&
         eta)
      ! calcuate eta on cell edges
      do idir=sdim,3
        ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
        ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
        ixCmin3=ixOmin3+kr(idir,3)-1;
        jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=0.d0
       do ix3=0,1
       do ix2=0,1
       do ix1=0,1
          if( ix1==1 .and. 1==idir  .or. ix2==1 .and. 2==idir  .or. ix3==1 &
             .and. 3==idir ) cycle
          ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;ixAmin3=ixCmin3+ix3;
          ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;ixAmax3=ixCmax3+ix3;
          jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idir)=jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
             idir)+eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)
       end do
       end do
       end do
        jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)=jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)*0.25d0
        jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)=jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)*jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)
      end do
    end if

    end associate
  end subroutine get_resistive_electric_field

  !> get ambipolar electric field on cell edges
  subroutine get_ambipolar_electric_field(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,fE)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)      :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,sdim:3)

    double precision :: jxbxb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:3)
    integer :: idir,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ix1,ix2,ix3

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    call mhd_get_jxbxb(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,jxbxb)
    ! calcuate electric field on cell edges from cell centers
    do idir=sdim,3
      !set electric field in jxbxb: E=nuA * jxbxb, where nuA=-etaA/rho^2
      !jxbxb(ixA^S,i) = -(mhd_eta_ambi/w(ixA^S, rho_)**2) * jxbxb(ixA^S,i)
      call multiplyAmbiCoef(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,jxbxb(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,idir),w,x)
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1
      ixCmin3=ixOmin3+kr(idir,3)-1;
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)=0.d0
     do ix3=0,1
     do ix2=0,1
     do ix1=0,1
        if( ix1==1 .and. 1==idir  .or. ix2==1 .and. 2==idir  .or. ix3==1 .and. &
           3==idir ) cycle
        ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;ixAmin3=ixCmin3+ix3;
        ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;ixAmax3=ixCmax3+ix3;
        fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
           idir)+jxbxb(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,idir)
     end do
     end do
     end do
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)*0.25d0
    end do

  end subroutine get_ambipolar_electric_field

  !> calculate cell-center values from face-center values
  subroutine mhd_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    type(state)                        :: s

    integer                            :: fxOmin1,fxOmin2,fxOmin3,fxOmax1,&
       fxOmax2,fxOmax3, gxOmin1,gxOmin2,gxOmin3,gxOmax1,gxOmax2,gxOmax3,&
        hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, jxOmin1,jxOmin2,&
       jxOmin3,jxOmax1,jxOmax2,jxOmax3, kxOmin1,kxOmin2,kxOmin3,kxOmax1,&
       kxOmax2,kxOmax3, idim

    associate(w=>s%w, ws=>s%ws)

    ! calculate cell-center values from face-center values in 2nd order
    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the w arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
      hxOmin3=ixOmin3-kr(idim,3);hxOmax1=ixOmax1-kr(idim,1)
      hxOmax2=ixOmax2-kr(idim,2);hxOmax3=ixOmax3-kr(idim,3);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idim))=half/s%surface(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)*(ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)*s%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)+ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         hxOmin3:hxOmax3,idim)*s%surfaceC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         hxOmin3:hxOmax3,idim))
    end do

    ! calculate cell-center values from face-center values in 4th order
    !do idim=1,ndim
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);

    !  ! Interpolate to cell barycentre using fourth order central formula
    !  w(ixO^S,mag(idim))=(0.0625d0/s%surface(ixO^S,idim))*&
    !         ( -ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !     +9.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !     +9.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !           -ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) )
    !end do

    ! calculate cell-center values from face-center values in 6th order
    !do idim=1,ndim
    !  fxO^L=ixO^L-3*kr(idim,^D);
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);
    !  kxO^L=ixO^L+2*kr(idim,^D);

    !  ! Interpolate to cell barycentre using sixth order central formula
    !  w(ixO^S,mag(idim))=(0.00390625d0/s%surface(ixO^S,idim))* &
    !     (  +3.0d0*ws(fxO^S,idim)*s%surfaceC(fxO^S,idim) &
    !       -25.0d0*ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !      +150.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !      +150.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !       -25.0d0*ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) &
    !        +3.0d0*ws(kxO^S,idim)*s%surfaceC(kxO^S,idim) )
    !end do

    end associate

  end subroutine mhd_face_to_center

  !> calculate magnetic field from vector potential
  subroutine b_from_vector_potential(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
     ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x)
    use mod_global_parameters
    use mod_constrained_transport

    integer, intent(in)                :: ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,&
       ixIsmax2,ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout)    :: ws(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:nws)
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    double precision                   :: Adummy(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,ixIsmin3:ixIsmax3,1:3)

    call b_from_vector_potentialA(ixIsmin1,ixIsmin2,ixIsmin3,ixIsmax1,ixIsmax2,&
       ixIsmax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ws, x, Adummy)

  end subroutine b_from_vector_potential

  subroutine Rfactor_from_temperature_ionization(w,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     Rfactor)
    use mod_global_parameters
    use mod_ionization_degree
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: Rfactor(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: iz_H(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       iz_He(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    call ionization_degree_from_temperature(ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,Te_),iz_H,iz_He)
    ! assume the first and second ionization of Helium have the same degree
    Rfactor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=(1.d0+iz_H(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+0.1d0*(1.d0+iz_He(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(1.d0+iz_He(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))))/(2.d0+3.d0*He_abundance)

  end subroutine Rfactor_from_temperature_ionization

  subroutine Rfactor_from_constant_ionization(w,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     Rfactor)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: Rfactor(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    Rfactor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=RR

  end subroutine Rfactor_from_constant_ionization

end module mod_mhd_phys
