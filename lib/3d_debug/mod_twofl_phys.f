!> Magneto-hydrodynamics module
module mod_twofl_phys

#include "amrvac.h"

  use mod_physics
  use mod_global_parameters, only: std_len
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  use mod_functions_bfield, only: get_divb,mag
  use mod_comm_lib, only: mpistop

  implicit none
  private
  !! E_c = E_kin + E_mag + E_int
  !! E_n = E_kin + E_int
  integer, public, parameter              :: EQ_ENERGY_TOT=2
  !! E_c = E_int
  !! E_n = E_int
  integer, public, parameter              :: EQ_ENERGY_INT=1
  !! E_n, E_c are calculated from density as c_adiab rho^gamma
  !! No energy equation => no variable assigned for it
  integer, public, parameter              :: EQ_ENERGY_NONE=0
  !! E_c = E_kin + E_int
  !! E_n = E_kin + E_int
  integer, public, parameter              :: EQ_ENERGY_KI=3

  integer, public, protected              :: twofl_eq_energy = EQ_ENERGY_TOT

  !> Whether hyperdiffusivity is used
  logical, public, protected              :: twofl_hyperdiffusivity = .false.
  logical, public, protected              :: twofl_dump_hyperdiffusivity_coef &
     = .false.
  double precision, public, protected, allocatable :: c_shk(:)
  double precision, public, protected, allocatable :: c_hyp(:)

  !> Whether thermal conduction is used
  logical, public, protected              :: twofl_thermal_conduction_c = &
     .false.
  !> type of TC used: 1: adapted module (mhd implementation), 2: adapted module (hd implementation)
  integer, parameter, private             :: MHD_TC =1
  integer, parameter, private             :: HD_TC =2
  integer, protected                      :: use_twofl_tc_c = MHD_TC

  !> Whether radiative cooling is added
  logical, public, protected              :: twofl_radiative_cooling_c = &
     .false.
  type(rc_fluid), public, allocatable :: rc_fl_c

  !> Whether viscosity is added
  logical, public, protected              :: twofl_viscosity = .false.

  !> Whether gravity is added: common flag for charges and neutrals
  logical, public, protected              :: twofl_gravity = .false.

  !> whether dump full variables (when splitting is used) in a separate dat file
  logical, public, protected              :: twofl_dump_full_vars = .false.

  !> Whether Hall-MHD is used
  logical, public, protected              :: twofl_Hall = .false.

  type(tc_fluid), public, allocatable :: tc_fl_c
  type(te_fluid), public, allocatable :: te_fl_c

  type(tc_fluid), allocatable :: tc_fl_n
  logical, public, protected              :: twofl_thermal_conduction_n = &
     .false.
  logical, public, protected              :: twofl_radiative_cooling_n = &
     .false.
  type(rc_fluid), allocatable :: rc_fl_n

  !> Whether TRAC method is used
  logical, public, protected              :: twofl_trac = .false.

  !> Whether GLM-MHD is used
  logical, public, protected              :: twofl_glm = .false.

  !> Which TRAC method is used            
  integer, public, protected              :: twofl_trac_type=1

  !> Height of the mask used in the TRAC method
  double precision, public, protected              :: twofl_trac_mask = 0.d0

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: twofl_glm_alpha = 0.5d0

  !> MHD fourth order
  logical, public, protected              :: twofl_4th_order = .false.

  !> Index of the density (in the w array)
  integer, public             :: rho_c_

  !> Indices of the momentum density
  integer, allocatable, public :: mom_c(:)

  !> Index of the energy density (-1 if not present)
  integer, public             :: e_c_=-1

  !> Index of the cutoff temperature for the TRAC method
  integer, public              :: Tcoff_c_
  integer, public              :: Tweight_c_

  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> equi vars flags
  logical, public :: has_equi_rho_c0 = .false.  
  logical, public :: has_equi_pe_c0 = .false.  

  !> equi vars indices in the state%equi_vars array
  integer, public :: equi_rho_c0_ = -1
  integer, public :: equi_pe_c0_ = -1
  logical, public                         :: twofl_equi_thermal_c = .false.

  logical, public                         :: twofl_equi_thermal = .false.
  !neutrals:

  integer, public              :: rho_n_
  integer, allocatable, public :: mom_n(:)
  integer, public              :: e_n_
  integer, public              :: Tcoff_n_
  integer, public              :: Tweight_n_
  logical, public :: has_equi_rho_n0 = .false. 
  logical, public :: has_equi_pe_n0 = .false.  
  integer, public :: equi_rho_n0_ = -1
  integer, public :: equi_pe_n0_ = -1

  ! related to collisions:
  !> collisional alpha
  double precision, public                :: twofl_alpha_coll = 0d0
  logical, public                         :: twofl_alpha_coll_constant = &
     .true.
  !> whether include thermal exchange collisional terms
  logical, public                         :: twofl_coll_inc_te = .true.
  !> whether include ionization/recombination inelastic collisional terms
  logical, public                         :: twofl_coll_inc_ionrec = .false.
  logical, public                         :: twofl_equi_thermal_n = .false.
  double precision, public                :: dtcollpar = -1d0 !negative value does not impose restriction on the timestep
  !> whether dump collisional terms in a separte dat file
  logical, public, protected              :: twofl_dump_coll_terms = .false.

  ! TODO Helium abundance not used, radiative cooling init uses it
  ! not in parameters list anymore 
  double precision, public, protected  :: He_abundance = 0d0
  ! two fluid is only H plasma
  double precision, public, protected  :: Rc = 2d0
  double precision, public, protected  :: Rn = 1d0

  !> The adiabatic index
  double precision, public                :: twofl_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: twofl_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public                :: twofl_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public                :: twofl_eta_hyper = 0.0d0

  !> The MHD Hall coefficient
  double precision, public                :: twofl_etah = 0.0d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'linde'

  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'uct_contact'

  !> Whether divB is computed with a fourth order approximation
  logical, public, protected :: twofl_divb_4thorder = .false.

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'

  !> clean initial divB
  logical, public :: clean_initial_divb     = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*3)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*3)=0

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.

  logical :: twofl_cbounds_species = .true.

  !> added from modules: gravity
  !> source split or not
  logical :: grav_split= .false.

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

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

  ! Public methods
  public :: twofl_phys_init
  public :: twofl_to_conserved
  public :: twofl_to_primitive
  public :: get_divb
  public :: get_rhoc_tot
  public :: twofl_get_v_c_idim
  ! TODO needed for the roe, see if can be used for n
  public :: twofl_get_csound2_c_from_conserved
  public :: get_rhon_tot
  public :: get_alpha_coll
  public :: get_gamma_ion_rec
  public :: twofl_get_v_n_idim
  public :: get_current
  public :: twofl_get_pthermal_c
  public :: twofl_get_pthermal_n
  public :: twofl_face_to_center
  public :: get_normalized_divb
  public :: b_from_vector_potential
  public :: usr_mask_gamma_ion_rec
  public :: usr_mask_alpha

  
  public :: twofl_clean_divb_multigrid
 

  abstract interface

    subroutine implicit_mult_factor_subroutine(ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
        step_dt, JJ, res)
      integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in) :: step_dt
      double precision, intent(in) :: JJ(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
      double precision, intent(out) :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)

    end subroutine implicit_mult_factor_subroutine

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

    subroutine mask_subroutine2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,res1, res2)
      use mod_global_parameters
      integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:nw)
      double precision, intent(inout) :: res1(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3),res2(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
    end subroutine mask_subroutine2

  end interface

   procedure (implicit_mult_factor_subroutine),&
       pointer :: calc_mult_factor => null()
   integer, protected ::  twofl_implicit_calc_mult_method = 1
  procedure(mask_subroutine), pointer  :: usr_mask_alpha => null()
  procedure(mask_subroutine2), pointer  :: usr_mask_gamma_ion_rec => null()

contains

  !> Read this module"s parameters from a file
  subroutine twofl_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /twofl_list/ twofl_eq_energy, twofl_gamma, twofl_adiab,twofl_eta,&
        twofl_eta_hyper, twofl_etah, twofl_glm_alpha,&
       twofl_thermal_conduction_c, use_twofl_tc_c, twofl_radiative_cooling_c,&
        twofl_Hall, twofl_gravity,twofl_viscosity, twofl_4th_order,&
        typedivbfix, source_split_divb, divbdiff,typedivbdiff, type_ct,&
        divbwave, SI_unit, B0field,B0field_forcefree, Bdip, Bquad, Boct, Busr,&
       twofl_equi_thermal_c,twofl_equi_thermal,twofl_dump_full_vars,&
        has_equi_rho_c0, has_equi_pe_c0, twofl_hyperdiffusivity,&
       twofl_dump_hyperdiffusivity_coef,has_equi_pe_n0, has_equi_rho_n0,&
        twofl_thermal_conduction_n, twofl_radiative_cooling_n,&
         twofl_alpha_coll,twofl_alpha_coll_constant,twofl_coll_inc_te,&
        twofl_coll_inc_ionrec,twofl_equi_thermal_n,dtcollpar,&
       twofl_dump_coll_terms,twofl_implicit_calc_mult_method,boundary_divbfix,&
        boundary_divbfix_skip, twofl_divb_4thorder, clean_initial_divb,&
         twofl_trac, twofl_trac_type, twofl_trac_mask,twofl_cbounds_species 

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, twofl_list, end=111)
111    close(unitpar)
    end do
 
  end subroutine twofl_read_params

  subroutine twofl_init_hyper(files)
    use mod_global_parameters
    use mod_hyperdiffusivity, only: hyperdiffusivity_init
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hyperdiffusivity_list/ c_shk, c_hyp

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hyperdiffusivity_list, end=113)
113    close(unitpar)
    end do
    
    call hyperdiffusivity_init()

    !!DEBUG
    if(mype .eq. 0) then
      print*, "Using Hyperdiffusivity"
      print*, "C_SHK ", c_shk(:)
      print*, "C_HYP ", c_hyp(:)
    endif

  end subroutine twofl_init_hyper

  !> Write this module's parameters to a snapsoht
  subroutine twofl_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = twofl_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine twofl_write_info

  subroutine twofl_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    !use mod_gravity, only: gravity_init
    use mod_supertimestepping, only: sts_init, add_sts_method,&
       set_conversion_methods_to_head, set_error_handling_to_head
    
    use mod_multigrid_coupling
   
    integer :: itr, idir

    call twofl_read_params(par_files)
    physics_type = "twofl"
    if (twofl_cbounds_species) then
      number_species = 2
    endif
    phys_energy=.true.
  !> Solve total energy equation or not
  ! for the two fluid the true value means 
  ! E_charges = E_mag + E_kin_charges + E_int_charges
  ! E_neutrals =  E_kin_neutrals + E_int_neutrals
    phys_total_energy=.false.

  !> Solve internal energy instead of total energy
  ! for the two fluid the true value means 
  ! E_charges = E_int_charges
  ! E_neutrals = E_int_neutrals
    phys_internal_e=.false.

  ! For the two fluid phys_energy=.true. and phys_internal_e=.false. and phys_total_energy = .false. means
  ! E_charges = E_kin_charges + E_int_charges
  ! E_neutrals =  E_kin_neutrals + E_int_neutrals
    phys_gamma = twofl_gamma

    if(twofl_eq_energy == EQ_ENERGY_INT) then
      phys_internal_e = .true.
    elseif(twofl_eq_energy == EQ_ENERGY_TOT) then
      phys_total_energy = .true.
    elseif(twofl_eq_energy == EQ_ENERGY_NONE) then
      phys_energy = .false.
    endif

    phys_trac=twofl_trac
    phys_trac_type=twofl_trac_type

    if(.not. phys_energy) then
      if(twofl_thermal_conduction_n) then
        twofl_thermal_conduction_n=.false.
        if(mype==0) write(*,*)&
            'WARNING: set twofl_thermal_conduction_n=F when twofl_energy=F'
      end if
      if(twofl_radiative_cooling_n) then
        twofl_radiative_cooling_n=.false.
        if(mype==0) write(*,*)&
            'WARNING: set twofl_radiative_cooling_n=F when twofl_energy=F'
      end if
      if(twofl_thermal_conduction_c) then
        twofl_thermal_conduction_c=.false.
        if(mype==0) write(*,*)&
            'WARNING: set twofl_thermal_conduction_c=F when twofl_energy=F'
      end if
      if(twofl_radiative_cooling_c) then
        twofl_radiative_cooling_c=.false.
        if(mype==0) write(*,*)&
            'WARNING: set twofl_radiative_cooling_c=F when twofl_energy=F'
      end if
      if(twofl_trac) then
        twofl_trac=.false.
        if(mype==0) write(*,*) 'WARNING: set twofl_trac=F when twofl_energy=F'
      end if
    end if
    
    if(twofl_trac .and. twofl_trac_type .le. 3) then
      twofl_trac_mask=bigdouble
      if(mype==0) write(*,*)&
          'WARNING: set twofl_trac_mask==bigdouble for global TRAC method'
    end if
    phys_trac_mask=twofl_trac_mask

    ! set default gamma for polytropic/isothermal process
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    
    case ('multigrid')
       type_divb = divb_multigrid
       use_multigrid = .true.
       mg%operator_type = mg_laplacian
       phys_global_source_after => twofl_clean_divb_multigrid
   
    case ('glm')
      twofl_glm          = .true.
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
      twofl_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_lindeglm
    case ('ct')
      type_divb = divb_ct
      stagger_grid = .true.
    case default
      call mpistop('Unknown divB fix')
    end select

    allocate(start_indices(number_species))
    allocate(stop_indices(number_species))
    start_indices(1)=1
    !allocate charges first and the same order as in mhd module
    rho_c_ = var_set_fluxvar("rho_c", "rho_c")
    !set variables from mod_variables to point to charges vars
    iw_rho = rho_c_

    allocate(mom_c(ndir))
    do idir=1,ndir
      mom_c(idir) = var_set_fluxvar("m_c","v_c",idir)
    enddo

    allocate(iw_mom(ndir))
    iw_mom(1:ndir) = mom_c(1:ndir)

    ! Set index of energy variable
    if (phys_energy) then
      e_c_ = var_set_fluxvar("e_c", "p_c")
      iw_e = e_c_
    else
      e_c_ = -1
    end if

  ! ambipolar sts assumes mag and energy charges are continuous
    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)

    if (twofl_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    ! set cutoff temperature when using the TRAC method, as well as an auxiliary weight
    Tweight_c_ = -1
    if(twofl_trac) then
      Tcoff_c_ = var_set_wextra()
      iw_tcoff = Tcoff_c_
      if(twofl_trac_type > 2) then
        Tweight_c_ = var_set_wextra()
      endif
    else
      Tcoff_c_ = -1
    end if

  !now allocate neutrals

    ! TODO so far number_species is only used to treat them differently
    ! in the solvers (different cbounds)
    if (twofl_cbounds_species) then
      stop_indices(1)=nwflux
      start_indices(2)=nwflux+1
    endif

    ! Determine flux variables
    rho_n_ = var_set_fluxvar("rho_n", "rho_n")
    allocate(mom_n(ndir))
    do idir=1,ndir
      mom_n(idir) = var_set_fluxvar("m_n","v_n",idir)
    enddo
    if (phys_energy) then
      e_n_ = var_set_fluxvar("e_n", "p_n")
    else
      e_n_     = -1
    end if

    Tweight_n_ = -1
    if(twofl_trac) then
      Tcoff_n_ = var_set_wextra()
      if(twofl_trac_type > 2) then
        Tweight_n_ = var_set_wextra()
      endif
    else
      Tcoff_n_ = -1
    end if

    stop_indices(number_species)=nwflux

    ! set indices of equi vars and update number_equi_vars
    number_equi_vars = 0
    if(has_equi_rho_n0) then
      number_equi_vars = number_equi_vars + 1
      equi_rho_n0_ = number_equi_vars
    endif  
    if(has_equi_pe_n0) then
      number_equi_vars = number_equi_vars + 1
      equi_pe_n0_ = number_equi_vars
      phys_equi_pe=.true.
    endif  
    if(has_equi_rho_c0) then
      number_equi_vars = number_equi_vars + 1
      equi_rho_c0_ = number_equi_vars
      iw_equi_rho = equi_rho_c0_
    endif  
    if(has_equi_pe_c0) then
      number_equi_vars = number_equi_vars + 1
      equi_pe_c0_ = number_equi_vars
      iw_equi_p = equi_pe_c0_
      phys_equi_pe=.true.
    endif  

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! determine number of stagger variables
    nws=ndim

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    if(ndim>1) then
      if(twofl_glm) then
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

    phys_get_dt              => twofl_get_dt
    phys_get_cmax            => twofl_get_cmax
    phys_get_a2max           => twofl_get_a2max
    !phys_get_tcutoff         => twofl_get_tcutoff_c
    if(twofl_cbounds_species) then
      if (mype .eq. 0) print*,&
          "Using different cbounds for each species nspecies = ",&
          number_species
      phys_get_cbounds         => twofl_get_cbounds_species
      phys_get_H_speed         => twofl_get_H_speed_species
    else
      if (mype .eq. 0) print*, "Using same cbounds for all species"
      phys_get_cbounds         => twofl_get_cbounds_one
      phys_get_H_speed         => twofl_get_H_speed_one
    endif
    phys_get_flux            => twofl_get_flux
    phys_add_source_geom     => twofl_add_source_geom
    phys_add_source          => twofl_add_source
    phys_to_conserved        => twofl_to_conserved
    phys_to_primitive        => twofl_to_primitive
    phys_check_params        => twofl_check_params
    phys_check_w             => twofl_check_w
    phys_write_info          => twofl_write_info
    phys_handle_small_values => twofl_handle_small_values
    !set equilibrium variables for the new grid
    if(number_equi_vars>0) then
      phys_set_equi_vars => set_equi_vars_grid
    endif
    ! convert_type is not known here, so associate the corresp. subroutine in check_params
    if(type_divb==divb_glm) then
      phys_modify_wLR => twofl_modify_wLR
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_get_ct_velocity => twofl_get_ct_velocity
      phys_update_faces => twofl_update_faces
      phys_face_to_center => twofl_face_to_center
      phys_modify_wLR => twofl_modify_wLR
    else if(ndim>1) then
      phys_boundary_adjust => twofl_boundary_adjust
    end if

    
    ! clean initial divb
    if(clean_initial_divb) phys_clean_divb => twofl_clean_divb_multigrid
   

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call twofl_physical_units()

    if(.not. phys_energy .and. (twofl_thermal_conduction_c.or. &
       twofl_thermal_conduction_n)) then
      call mpistop("thermal conduction needs twofl_energy=T")
    end if

    ! initialize thermal conduction module
    if (twofl_thermal_conduction_c .or. twofl_thermal_conduction_n) then
      phys_req_diagonal = .true.
      call sts_init()
      call tc_init_params(twofl_gamma)
    endif
    if (twofl_thermal_conduction_c) then
      allocate(tc_fl_c)
      if(has_equi_pe_c0 .and. has_equi_rho_c0) then
        tc_fl_c%get_temperature_from_eint => &
           twofl_get_temperature_from_eint_c_with_equi
        if(phys_internal_e) then
          tc_fl_c%get_temperature_from_conserved => &
             twofl_get_temperature_from_eint_c_with_equi
        else
            if(twofl_eq_energy == EQ_ENERGY_KI) then
              tc_fl_c%get_temperature_from_conserved => &
                 twofl_get_temperature_from_eki_c_with_equi
            else
              tc_fl_c%get_temperature_from_conserved => &
                 twofl_get_temperature_from_etot_c_with_equi
            endif
        endif
        if(twofl_equi_thermal_c) then
          tc_fl_c%has_equi = .true.
          tc_fl_c%get_temperature_equi => twofl_get_temperature_c_equi
          tc_fl_c%get_rho_equi => twofl_get_rho_c_equi
        else  
          tc_fl_c%has_equi = .false.
        endif
      else
        if(phys_internal_e) then
          tc_fl_c%get_temperature_from_conserved => &
             twofl_get_temperature_from_eint_c
         else
            if(twofl_eq_energy == EQ_ENERGY_KI) then
              tc_fl_c%get_temperature_from_conserved => &
                 twofl_get_temperature_from_eki_c
            else
              tc_fl_c%get_temperature_from_conserved => &
                 twofl_get_temperature_from_etot_c
          endif  
         endif  
        tc_fl_c%get_temperature_from_eint => twofl_get_temperature_from_eint_c
      endif
      if(use_twofl_tc_c .eq. MHD_TC) then
        call tc_get_mhd_params(tc_fl_c,tc_c_params_read_mhd)
        call add_sts_method(twofl_get_tc_dt_mhd_c,&
           twofl_sts_set_source_tc_c_mhd,e_c_,1,e_c_,1,.false.)
      else if(use_twofl_tc_c .eq. HD_TC) then
        call tc_get_hd_params(tc_fl_c,tc_c_params_read_hd)
        call add_sts_method(twofl_get_tc_dt_hd_c,twofl_sts_set_source_tc_c_hd,&
           e_c_,1,e_c_,1,.false.)
      endif
      if(.not. phys_internal_e) then
        call set_conversion_methods_to_head(twofl_e_to_ei_c, twofl_ei_to_e_c)
      endif
      call set_error_handling_to_head(twofl_tc_handle_small_e_c)
      tc_fl_c%get_rho => get_rhoc_tot
      tc_fl_c%e_ = e_c_
      tc_fl_c%Tcoff_ = Tcoff_c_
    end if
    if (twofl_thermal_conduction_n) then
      allocate(tc_fl_n)
      call tc_get_hd_params(tc_fl_n,tc_n_params_read_hd)
      if(has_equi_pe_n0 .and. has_equi_rho_n0) then
        tc_fl_n%get_temperature_from_eint => &
           twofl_get_temperature_from_eint_n_with_equi
        if(twofl_equi_thermal_n) then
          tc_fl_n%has_equi = .true.
          tc_fl_n%get_temperature_equi => twofl_get_temperature_n_equi
          tc_fl_n%get_rho_equi => twofl_get_rho_n_equi
        else  
          tc_fl_n%has_equi = .false.
        endif
      else
        tc_fl_n%get_temperature_from_eint => twofl_get_temperature_from_eint_n
      endif
      if(phys_internal_e) then
        if(has_equi_pe_n0 .and. has_equi_rho_n0) then
          tc_fl_n%get_temperature_from_conserved => &
             twofl_get_temperature_from_eint_n_with_equi
        else
          tc_fl_n%get_temperature_from_conserved => &
             twofl_get_temperature_from_eint_n
        endif 
        call add_sts_method(twofl_get_tc_dt_hd_n,twofl_sts_set_source_tc_n_hd,&
           e_n_,1,e_n_,1,.false.)
      else
        if(has_equi_pe_n0 .and. has_equi_rho_n0) then
          tc_fl_n%get_temperature_from_conserved => &
             twofl_get_temperature_from_etot_n_with_equi
        else
          tc_fl_n%get_temperature_from_conserved => &
             twofl_get_temperature_from_etot_n
        endif 
        call add_sts_method(twofl_get_tc_dt_hd_n,twofl_sts_set_source_tc_n_hd,&
           e_n_,1,e_n_,1,.false.)
        call set_conversion_methods_to_head(twofl_e_to_ei_n, twofl_ei_to_e_n)
      endif
      call set_error_handling_to_head(twofl_tc_handle_small_e_n)
      tc_fl_n%get_rho => get_rhon_tot
      tc_fl_n%e_ = e_n_
      tc_fl_n%Tcoff_ = Tcoff_n_
    end if


    if(.not. phys_energy .and. (twofl_radiative_cooling_c.or. &
       twofl_radiative_cooling_n)) then
      call mpistop("radiative cooling needs twofl_energy=T")
    end if

    if(twofl_equi_thermal .and. (.not. has_equi_pe_c0 .or. .not.  &
       has_equi_pe_n0)) then
      call mpistop("twofl_equi_thermal=T has_equi_pe_n0 and has _equi_pe_c0=T")

    endif

    ! initialize thermal conduction module
    if (twofl_radiative_cooling_c .or. twofl_radiative_cooling_n) then
    ! Initialize radiative cooling module
      call radiative_cooling_init_params(twofl_gamma,He_abundance)
      if(twofl_radiative_cooling_c) then
        allocate(rc_fl_c)
        call radiative_cooling_init(rc_fl_c,rc_params_read_c)
        rc_fl_c%get_rho => get_rhoc_tot
        rc_fl_c%get_pthermal => twofl_get_pthermal_c
        rc_fl_c%get_var_Rfactor => Rfactor_c
        rc_fl_c%e_ = e_c_
        rc_fl_c%Tcoff_ = Tcoff_c_
        if(has_equi_pe_c0 .and. has_equi_rho_c0 .and. twofl_equi_thermal_c) &
           then
          rc_fl_c%has_equi = .true.
          rc_fl_c%get_rho_equi => twofl_get_rho_c_equi
          rc_fl_c%get_pthermal_equi => twofl_get_pe_c_equi
        else
          rc_fl_c%has_equi = .false.
        end if
      end if
    end if
    allocate(te_fl_c)
    te_fl_c%get_rho=> get_rhoc_tot
    te_fl_c%get_pthermal=> twofl_get_pthermal_c
    te_fl_c%get_var_Rfactor => Rfactor_c

    phys_te_images => twofl_te_images


    ! Initialize viscosity module
    !!TODO
    !if (twofl_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if(twofl_gravity) then
    !  call gravity_init()
       call grav_params_read(par_files)
    end if

    ! Initialize particles module
    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in getflux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if (twofl_hall) then
       phys_req_diagonal = .true.
       if (twofl_4th_order) then
          phys_wider_stencil = 2
       else
          phys_wider_stencil = 1
       end if
    end if

    if(twofl_hyperdiffusivity) then
      allocate(c_shk(1:nwflux))
      allocate(c_hyp(1:nwflux))
      call twofl_init_hyper(par_files)
    end if

  end subroutine twofl_phys_init


  subroutine twofl_te_images
    use mod_global_parameters
    use mod_thermal_emission

    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_c)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_c)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_c)
      case('WIvtiCCmpi','WIvtuCCmpi')
        call get_whitelight_image(unitconvert,te_fl_c)
      case default
        call mpistop("Error in synthesize emission: Unknown convert_type")
      end select
  end subroutine twofl_te_images
  

  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  twofl_sts_set_source_tc_c_mhd(ixImin1,ixImin2,ixImin3,ixImax1,&
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
       fix_conserve_at_step,my_dt,igrid,nflux,tc_fl_c)
  end subroutine twofl_sts_set_source_tc_c_mhd

  subroutine  twofl_sts_set_source_tc_c_hd(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
     fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_hd
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, igrid, nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_hd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
       fix_conserve_at_step,my_dt,igrid,nflux,tc_fl_c)
  end subroutine twofl_sts_set_source_tc_c_hd

  function twofl_get_tc_dt_mhd_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,&
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
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,x,tc_fl_c) 
  end function twofl_get_tc_dt_mhd_c

  function twofl_get_tc_dt_hd_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,&
     x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_hd
 
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_hd(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,x,tc_fl_c) 
  end function twofl_get_tc_dt_hd_c

  subroutine twofl_tc_handle_small_e_c(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, step)
    use mod_global_parameters
    use mod_small_values

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    integer, intent(in)    :: step

    character(len=140) :: error_msg

    write(error_msg,"(a,i3)") "Charges thermal conduction step ", step
    call twofl_handle_small_ei_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,e_c_,error_msg)
  end subroutine twofl_tc_handle_small_e_c

  subroutine  twofl_sts_set_source_tc_n_hd(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
     fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_hd
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, igrid, nflux
    double precision, intent(in) ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) ::  wres(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_hd(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,wres,&
       fix_conserve_at_step,my_dt,igrid,nflux,tc_fl_n)
  end subroutine twofl_sts_set_source_tc_n_hd

  subroutine twofl_tc_handle_small_e_n(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
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

    write(error_msg,"(a,i3)") "Neutral thermal conduction step ", step
    call twofl_handle_small_ei_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,e_n_,error_msg)
  end subroutine twofl_tc_handle_small_e_n

  function twofl_get_tc_dt_hd_n(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,&
     x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_hd
 
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_hd(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dx1,dx2,dx3,x,tc_fl_n) 
  end function twofl_get_tc_dt_hd_n

  subroutine tc_n_params_read_hd(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_global_parameters, only: unitpar
    type(tc_fluid), intent(inout) :: fl
    integer                      :: n
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0

    namelist /tc_n_list/ tc_saturate, tc_k_para

    do n = 1, size(par_files)
       open(unitpar, file=trim(par_files(n)), status="old")
       read(unitpar, tc_n_list, end=111)
111      close(unitpar)
    end do
    fl%tc_saturate = tc_saturate
    fl%tc_k_para = tc_k_para

  end subroutine tc_n_params_read_hd

  subroutine rc_params_read_n(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_constants, only: bigdouble
    type(rc_fluid), intent(inout) :: fl
    integer                      :: n
    ! list parameters
    integer :: ncool = 4000
    double precision :: cfrac=0.1d0
  
    !> Name of cooling curve
    character(len=std_len)  :: coolcurve='JCorona'
  
    !> Name of cooling method
    character(len=std_len)  :: coolmethod='exact'
  
    !> Fixed temperature not lower than tlow
    logical    :: Tfix=.false.
  
    !> Lower limit of temperature
    double precision   :: tlow=bigdouble
  
    !> Add cooling source in a split way (.true.) or un-split way (.false.)
    logical    :: rc_split=.false.

    namelist /rc_list_n/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix,&
        rc_split

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, rc_list_n, end=111)
111     close(unitpar)
    end do

    fl%ncool=ncool
    fl%coolcurve=coolcurve
    fl%coolmethod=coolmethod
    fl%tlow=tlow
    fl%Tfix=Tfix
    fl%rc_split=rc_split
    fl%cfrac=cfrac
  end subroutine rc_params_read_n

  !end wrappers

  ! fill in tc_fluid fields from namelist
  subroutine tc_c_params_read_mhd(fl)
    use mod_global_parameters, only: unitpar,par_files
    type(tc_fluid), intent(inout) :: fl

    integer                      :: n

    ! list parameters
    logical :: tc_perpendicular=.false.
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0
    double precision :: tc_k_perp=0d0
    character(len=std_len)  :: tc_slope_limiter="MC"

    namelist /tc_c_list/ tc_perpendicular, tc_saturate, tc_slope_limiter,&
        tc_k_para, tc_k_perp
    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, tc_c_list, end=111)
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
  end subroutine tc_c_params_read_mhd

  subroutine tc_c_params_read_hd(fl)
    use mod_global_parameters, only: unitpar,par_files
    use mod_global_parameters, only: unitpar
    type(tc_fluid), intent(inout) :: fl
    integer                      :: n
    logical :: tc_saturate=.false.
    double precision :: tc_k_para=0d0

    namelist /tc_c_list/ tc_saturate, tc_k_para

    do n = 1, size(par_files)
       open(unitpar, file=trim(par_files(n)), status="old")
       read(unitpar, tc_c_list, end=111)
111      close(unitpar)
    end do
    fl%tc_saturate = tc_saturate
    fl%tc_k_para = tc_k_para

  end subroutine tc_c_params_read_hd

!! end th cond

!!rad cool
  subroutine rc_params_read_c(fl)
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


    namelist /rc_list_c/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix,&
        rc_split

    do n = 1, size(par_files)
      open(unitpar, file=trim(par_files(n)), status="old")
      read(unitpar, rc_list_c, end=111)
111     close(unitpar)
    end do

    fl%ncool=ncool
    fl%coolcurve=coolcurve
    fl%coolmethod=coolmethod
    fl%tlow=tlow
    fl%Tfix=Tfix
    fl%rc_split=rc_split
    fl%cfrac=cfrac
  end subroutine rc_params_read_c

!! end rad cool

  !> sets the equilibrium variables
  subroutine set_equi_vars_grid_faces(igrid,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    use mod_global_parameters
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

  ! w, wnew conserved
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

    call  get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3))
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_n_) = rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom_n(:)) =  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom_n(:))
    call  get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3))
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_c_) = rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom_c(:)) =  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom_c(:))

    if (B0field) then
      ! add background magnetic field B0 to B
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0)
    else 
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(:))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))
    end if

    if(phys_energy) then
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_n_)
      if(has_equi_pe_n0) then
        wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_) = wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_n0_,0)* inv_gamma_1  
      endif
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_c_)
      if(has_equi_pe_c0) then
        wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_) = wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_c0_,0)* inv_gamma_1  
      endif
      if(B0field .and. phys_total_energy) then
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)+0.5d0*sum(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,:,0)**2,dim=ndim+1) + sum(w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))*block%B0(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0),dim=ndim+1)
      endif
    endif

  end function convert_vars_splitting

  !> copied from mod_gravity
  subroutine grav_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /grav_list/ grav_split

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, grav_list, end=111)
111    close(unitpar)
    end do

  end subroutine grav_params_read

  subroutine associate_dump_hyper()
    use mod_global_parameters
    use mod_convert, only: add_convert_method
    integer :: ii
    do ii = 1,ndim 
      if(ii==1) then
        call add_convert_method(dump_hyperdiffusivity_coef_x, nw,&
            cons_wnames(1:nw), "hyper_x")
      elseif(ii==2) then
        call add_convert_method(dump_hyperdiffusivity_coef_y, nw,&
            cons_wnames(1:nw), "hyper_y")
      else
        call add_convert_method(dump_hyperdiffusivity_coef_z, nw,&
            cons_wnames(1:nw), "hyper_z")
      endif
    enddo
  end subroutine associate_dump_hyper

  subroutine twofl_check_params
    use mod_global_parameters
    use mod_usr_methods
    use mod_convert, only: add_convert_method

    ! after user parameter setting
    gamma_1=twofl_gamma-1.d0
    if (.not. phys_energy) then
       if (twofl_gamma <= 0.0d0) call mpistop ("Error: twofl_gamma <= 0")
       if (twofl_adiab < 0.0d0) call mpistop ("Error: twofl_adiab < 0")
       small_pressure = twofl_adiab*small_density**twofl_gamma
    else
       if (twofl_gamma <= 0.0d0 .or. twofl_gamma == 1.0d0) call mpistop &
          ("Error: twofl_gamma <= 0 or twofl_gamma == 1")
       inv_gamma_1=1.d0/gamma_1
       small_e = small_pressure * inv_gamma_1
    end if

    ! this has to be done here as use_imex_scheme is not set in init subroutine, 
    ! but here it is
    if(use_imex_scheme) then
      if(has_collisions()) then
        ! implicit collisional terms update
        phys_implicit_update => twofl_implicit_coll_terms_update
        phys_evaluate_implicit => twofl_evaluate_implicit
        if(mype .eq. 1) then
            print*, "IMPLICIT UPDATE with calc_mult_factor", twofl_implicit_calc_mult_method
        endif
        if(twofl_implicit_calc_mult_method == 1) then
          calc_mult_factor => calc_mult_factor1
        else
          calc_mult_factor => calc_mult_factor2
        endif
      endif
    else
      ! check dtcoll par for explicit implementation of the coll. terms
      if(dtcollpar .le. 0d0 .or. dtcollpar .ge. 1d0) then
        if (mype .eq. 0) print*,&
"Explicit update of coll terms requires 0<dtcollpar<1, dtcollpar set to 0.8."
        dtcollpar = 0.8
      endif 
        
    endif
!    if(H_ion_fr == 0d0 .and. He_ion_fr == 0d0) then
!      call mpistop("H_ion_fr or He_ion_fr must be > 0 or use hd module")
!    endif
!    if(H_ion_fr == 1d0 .and. He_ion_fr == 1d0) then
!      call mpistop("H_ion_fr or He_ion_fr must be < 1 or use mhd module")
!    endif
    if (number_equi_vars > 0 .and. .not. associated(usr_set_equi_vars)) then
      call mpistop("usr_set_equi_vars has to be implemented in the user file")
    endif
    if(convert .or. autoconvert) then
      if(convert_type .eq. 'dat_generic_mpi') then
        if(twofl_dump_full_vars) then
          if(mype .eq. 0) print*, " add conversion method: split -> full "
          call add_convert_method(convert_vars_splitting, nw, cons_wnames,&
              "new")
        endif
        if(twofl_dump_coll_terms) then
          if(mype .eq. 0) print*, " add conversion method: dump coll terms "
          call add_convert_method(dump_coll_terms, 3, (/"alpha    ",&
              "gamma_rec", "gamma_ion"/), "_coll")
        endif
        if(twofl_hyperdiffusivity .and. twofl_dump_hyperdiffusivity_coef) then
          if(mype .eq. 0) print*,&
              " add conversion method: dump hyperdiffusivity coeff. "
          call associate_dump_hyper()
        endif
      endif
    endif
  end subroutine twofl_check_params

  subroutine twofl_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0,c_lightspeed
    !double precision :: a,b,c,d
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
      miu0=4.d0*dpi
      c_lightspeed=const_c
    end if

    a=1d0  
    b=1d0
    Rc=2d0
    Rn=1d0  

    !now the unit choice:
    !unit 1 from number density or density -> mH
    !unit 2 from 

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
  end subroutine twofl_physical_units

  subroutine twofl_check_w(primitive,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
        
    if(has_equi_rho_n0) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_n_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho_n0_,0)
    else  
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_n_) 
    endif
    where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_n_) = .true.
    if(has_equi_rho_c0) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_c_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho_c0_,0)
    else  
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_c_) 
    endif
    where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_c_) = .true.
    if(phys_energy) then
      if(primitive) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_n_)
        if(has_equi_pe_n0) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe_n0_,0)
        endif
        where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) < small_pressure) flag(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_n_) = .true.
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_c_)
        if(has_equi_pe_c0) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe_c0_,0)
        endif
        where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) < small_pressure) flag(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_c_) = .true.
      else
        if(phys_internal_e) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_n_)
          if(has_equi_pe_n0) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_n0_,0)*inv_gamma_1
          endif
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_n_) = .true.
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_c_)
          if(has_equi_pe_c0) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_c0_,0)*inv_gamma_1
          endif
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_c_) = .true.
        else
          !neutrals
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_n_)-twofl_kin_en_n(w,ixImin1,ixImin2,ixImin3,&
             ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
             ixOmax3)
          if(has_equi_pe_n0) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_n0_,0)*inv_gamma_1
          endif
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_n_) = .true.
          if(phys_total_energy) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,e_c_)-twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,&
               ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
               ixOmax3)-twofl_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
          else
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,e_c_)-twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,&
               ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
               ixOmax3)
          end if
          if(has_equi_pe_c0) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)+block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_c0_,0)*inv_gamma_1
          endif
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,e_c_) = .true.
        end if
      endif
    end if

  end subroutine twofl_check_w

  !> Transform primitive variables into conservative ones
  subroutine twofl_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    integer                         :: idir
    double precision                :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    !if (fix_small_values) then
    !  call twofl_handle_small_values(.true., w, x, ixI^L, ixO^L, 'twofl_to_conserved')
    !end if

    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(phys_energy) then
      if(phys_internal_e) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)*inv_gamma_1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)*inv_gamma_1
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_n(:))**2,dim=ndim+1)*rhon(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        if(phys_total_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mom_c(:))**2,dim=ndim+1)*rhoc(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)+twofl_mag_en(w, ixImin1,ixImin2,&
             ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
             ixOmax2,ixOmax3)
        else
          ! kinetic energy + internal energy is evolved   
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mom_c(:))**2,dim=ndim+1)*rhoc(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end if
      end if
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idir)) = rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idir))
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir)) = rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir))
    end do
  end subroutine twofl_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine twofl_to_primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    integer                         :: idir
    double precision                :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if (fix_small_values) then
      call twofl_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3, 'twofl_to_primitive')
    end if

    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)

    if(phys_energy) then
      if(phys_internal_e) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)*gamma_1
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)*gamma_1
      else
        ! neutrals evolved energy = ke + e_int 
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)-twofl_kin_en_n(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3))
        ! charges
        if(phys_total_energy) then
         ! evolved energy = ke + e_int + e_mag 
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)-twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
             ixOmax3)-twofl_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3))
        else
         ! evolved energy = ke + e_int 
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)-twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3))
        end if
      end if
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir))/rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idir))/rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

  end subroutine twofl_to_primitive

!!USED IN TC
  !> Transform internal energy to total energy
  subroutine twofl_ei_to_e_c(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
 
    ! Calculate total energy from internal, kinetic and magnetic energy
    if(twofl_eq_energy == EQ_ENERGY_KI) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)+twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    else
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)+twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3)+twofl_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    endif
  end subroutine twofl_ei_to_e_c

  !> Transform total energy to internal energy
  subroutine twofl_e_to_ei_c(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    if(twofl_eq_energy == EQ_ENERGY_KI) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)-twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    else
    ! Calculate ei = e - ek - eb
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)-twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3)-twofl_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    endif
  end subroutine twofl_e_to_ei_c

   !Neutrals 
  subroutine twofl_ei_to_e_n(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
 
    ! Calculate total energy from internal and kinetic  energy

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_n_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_n_)+twofl_kin_en_n(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3)

  end subroutine twofl_ei_to_e_n

  !> Transform total energy to internal energy
  subroutine twofl_e_to_ei_n(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    ! Calculate ei = e - ek 
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_n_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_n_)-twofl_kin_en_n(w,ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3)

    call twofl_handle_small_ei_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,e_n_,&
       "e_to_ei_n")
  end subroutine twofl_e_to_ei_n

  subroutine twofl_handle_small_values(primitive, w, x, ixImin1,ixImin2,&
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
    double precision :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call twofl_check_w(primitive, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, flag)

    if(any(flag)) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_rho_c0) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_c_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_c_) = small_density-block%equi_vars(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_c0_,0)
        else
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_c_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_c_) = small_density
        end if
        if(has_equi_rho_n0) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_n_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_n_) = small_density-block%equi_vars(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_n0_,0)
        else
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_n_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             rho_n_) = small_density
        end if
        do idir = 1, ndir
          if(small_values_fix_iw(mom_n(idir))) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               rho_n_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom_n(idir)) = 0.0d0
          end if
          if(small_values_fix_iw(mom_c(idir))) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               rho_c_)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom_c(idir)) = 0.0d0
          end if
        end do

        if(phys_energy) then
          if(primitive) then
            if(has_equi_pe_n0) then 
             tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = small_pressure - &
                block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,equi_pe_n0_,0)
            else
             tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = small_pressure
            end if  
            if(has_equi_pe_c0) then 
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = small_e - block%equi_vars(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_c0_,0)
            else
             tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = small_pressure
            end if
          else
            ! conserved  
            if(has_equi_pe_n0) then 
              tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = small_e - block%equi_vars(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_n0_,0)*inv_gamma_1 
            else
              tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = small_e
            end if  
            if(has_equi_pe_c0) then 
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = small_e - block%equi_vars(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_c0_,0)*inv_gamma_1 
            else
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = small_e
            end if  
            if(phys_internal_e) then
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_n_))
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_n_)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
              end where
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_c_))
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_c_)=tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
              end where
            else
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_n_))
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_n_) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)+twofl_kin_en_n(w,ixImin1,ixImin2,ixImin3,&
                   ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                   ixOmax2,ixOmax3)
              end where
              if(phys_total_energy) then
                where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_c_))
                  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     e_c_) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)+twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,&
                     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                     ixOmax2,ixOmax3)+twofl_mag_en(w,ixImin1,ixImin2,ixImin3,&
                     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                     ixOmax2,ixOmax3)
                end where
              else
                where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_c_))
                  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     e_c_) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     ixOmin3:ixOmax3)+twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,&
                     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
                     ixOmax2,ixOmax3)
                end where
              end if
            end if
          end if
        end if
      case ("average")
        call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
            flag)
      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek-eb)
          if(phys_energy) then
            if(phys_internal_e) then
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_c_)*gamma_1
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_n_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_n_)*gamma_1
            else
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_n_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,e_n_)-twofl_kin_en_n(w,ixImin1,ixImin2,&
                 ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
                 ixOmax1,ixOmax2,ixOmax3))
              if(phys_total_energy) then
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   e_c_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,e_c_)-twofl_kin_en_c(w,ixImin1,ixImin2,&
                   ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
                   ixOmax1,ixOmax2,ixOmax3)-twofl_mag_en(w,ixImin1,ixImin2,&
                   ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
                   ixOmax1,ixOmax2,ixOmax3))
               else
                 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                    e_c_)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    ixOmin3:ixOmax3,e_c_)-twofl_kin_en_c(w,ixImin1,ixImin2,&
                    ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
                    ixOmax1,ixOmax2,ixOmax3))

               end if  
            end if
          end if
          ! Convert momentum to velocity
          if(has_equi_rho_n0) then
            tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_n_) + block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_n0_,0)
          else  
            tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_n_) 
          end if
    
          if(has_equi_rho_c0) then
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_c_) + block%equi_vars(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_c0_,0)
          else  
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,rho_c_) 
          end if
          do idir = 1, ndir
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom_n(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mom_n(idir))/tmp1(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom_c(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mom_c(idir))/tmp2(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end do
        end if
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if
  end subroutine twofl_handle_small_values

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine twofl_get_cmax(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: vc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: cmax2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: vn(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call twofl_get_csound_c_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,cmax)
    call twofl_get_v_c_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,vc)
    call twofl_get_v_n_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,vn)
    call twofl_get_csound_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,cmax2)
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=max(abs(vn(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))+cmax2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),abs(vc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))+cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
        
  end subroutine twofl_get_cmax

  subroutine twofl_get_a2max(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a2max)
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
  end subroutine twofl_get_a2max

  ! COPIED from hd/moh_hd_phys
  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine twofl_get_tcutoff_n(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,tco_local,&
     Tmax_local)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim),w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(out) :: tco_local, Tmax_local

    double precision, parameter :: delta=0.25d0
    double precision :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       Te(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),lts(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,hxOmin1,hxOmin2,&
       hxOmin3,hxOmax1,hxOmax2,hxOmax3
    logical :: lrlt(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    
  end subroutine twofl_get_tcutoff_n

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine twofl_get_tcutoff_c(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,Tco_local,&
     Tmax_local)
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
    double precision :: ltr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       ltrc,ltrp,altr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: idims,jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3,hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,ixA1,ixA2,ixA3,ixB1,ixB2,ixB3
    integer :: jxPmin1,jxPmin2,jxPmin3,jxPmax1,jxPmax2,jxPmax3,hxPmin1,hxPmin2,&
       hxPmin3,hxPmax1,hxPmax2,hxPmax3,ixPmin1,ixPmin2,ixPmin3,ixPmax1,ixPmax2,&
       ixPmax3
    logical :: lrlt(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    ! reuse lts as rhoc
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,lts)
    if(phys_internal_e) then
      tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,e_c_)
    else
      tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,e_c_)-0.5d0*(sum(w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,mom_c(:))**2,&
         dim=ndim+1)/lts(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)+sum(w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,mag(:))**2,dim=ndim+1))
    end if
    Te(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=tmp1(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)/lts(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)*(twofl_gamma-1.d0)
    Tmax_local=maxval(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

    
    
    select case(twofl_trac_type)
    case(0)
      !> test case, fixed cutoff temperature
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         Tcoff_c_)=2.5d5/unit_temperature
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
      if(twofl_trac_type .gt. 1) then
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
      ixPmin1=ixOmin1-1;ixPmin2=ixOmin2-1;ixPmin3=ixOmin3-1;ixPmax1=ixOmax1+1
      ixPmax2=ixOmax2+1;ixPmax3=ixOmax3+1;
      ! temperature gradient at cell centers
      do idims=1,ndim
        call gradient(Te,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixPmin1,ixPmin2,ixPmin3,ixPmax1,ixPmax2,ixPmax3,idims,tmp1)
        gradT(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
           idims)=tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)
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
         ixPmin3:ixPmax3)=dsqrt(sum(bunitvec(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3,:)**2,dim=ndim+1))
      where(tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)/=0.d0)
        tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)=1.d0/tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
           ixPmin3:ixPmax3)
      elsewhere
        tmp1(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)=bigdouble
      end where
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
      ltr(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3)=max(one,&
          (exp(lts(ixPmin1:ixPmax1,ixPmin2:ixPmax2,&
         ixPmin3:ixPmax3))/ltrc)**ltrp)
  
      altr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=zero  
      do idims=1,ndim 
        hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
        hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
        hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
        jxOmin1=ixOmin1+kr(idims,1);jxOmin2=ixOmin2+kr(idims,2)
        jxOmin3=ixOmin3+kr(idims,3);jxOmax1=ixOmax1+kr(idims,1)
        jxOmax2=ixOmax2+kr(idims,2);jxOmax3=ixOmax3+kr(idims,3);
        altr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=altr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) +0.25*(ltr(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3)+two*ltr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+ltr(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3))*bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idims)**2
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           Tcoff_c_)=Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*altr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)**(0.4*ltrp)
      end do
    case(3,5)
      !> do nothing here
    case default
      call mpistop("unknown twofl_trac_type")
    end select 
   
  end subroutine twofl_get_tcutoff_c

  !> get H speed for H-correction to fix the carbuncle problem at grid-aligned shock front
  subroutine twofl_get_H_speed_one(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     Hspeed) 
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
      call twofl_get_csound_prim(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,id,&
         tmp)
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
       1)=0.5d0*abs(0.5d0 * (wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
       jxCmin3:jxCmax3,mom_c(idim))+ wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
       jxCmin3:jxCmax3,mom_n(idim))) +csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
       jxCmin3:jxCmax3,idim)- 0.5d0 * (wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
       ixCmin3:ixCmax3,mom_c(idim)) + wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
       ixCmin3:ixCmax3,mom_n(idim)))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
       ixCmin3:ixCmax3,idim))

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=ixCmax1+kr(id,1);ixAmax2=ixCmax2+kr(id,2)
      ixAmax3=ixCmax3+kr(id,3);
      ixAmin1=ixCmin1+kr(id,1);ixAmin2=ixCmin2+kr(id,2)
      ixAmin3=ixCmin3+kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(0.5d0 * (wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,mom_c(id)) + wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,mom_n(id)))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,id)-0.5d0 * (wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,mom_c(id)) + wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,mom_n(id)))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,id)))


      ixAmax1=ixCmax1-kr(id,1);ixAmax2=ixCmax2-kr(id,2)
      ixAmax3=ixCmax3-kr(id,3);
      ixAmin1=ixCmin1-kr(id,1);ixAmin2=ixCmin2-kr(id,2)
      ixAmin3=ixCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(0.5d0 * (wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,mom_c(id)) + wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,mom_n(id)))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ixCmin3:ixCmax3,id)-0.5d0 * (wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,mom_c(id)) + wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,mom_n(id)))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,id)))

    end do

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=jxCmax1+kr(id,1);ixAmax2=jxCmax2+kr(id,2)
      ixAmax3=jxCmax3+kr(id,3);
      ixAmin1=jxCmin1+kr(id,1);ixAmin2=jxCmin2+kr(id,2)
      ixAmin3=jxCmin3+kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(0.5d0 * (wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,mom_c(id)) + wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,mom_n(id)))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,id)-0.5d0 * (wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3,mom_c(id)) + wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3,mom_n(id)))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3,id)))
      ixAmax1=jxCmax1-kr(id,1);ixAmax2=jxCmax2-kr(id,2)
      ixAmax3=jxCmax3-kr(id,3);
      ixAmin1=jxCmin1-kr(id,1);ixAmin2=jxCmin2-kr(id,2)
      ixAmin3=jxCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(0.5d0 * (wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3,mom_c(id)) + wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3,mom_n(id)))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
         jxCmin3:jxCmax3,id)-0.5d0 * (wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,mom_c(id)) + wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,mom_n(id)))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
         ixAmin3:ixAmax3,id)))
    end do

  end subroutine twofl_get_H_speed_one

  !> get H speed for H-correction to fix the carbuncle problem at grid-aligned shock front
  subroutine twofl_get_H_speed_species(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     Hspeed) 
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
    ! charges
    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    do id=1,ndim
      call twofl_get_csound_prim_c(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,id,&
         tmp)
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
       mom_c(idim))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
       idim)-wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       mom_c(idim))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
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
         mom_c(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)-wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom_c(id))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         id)))
      ixAmax1=ixCmax1-kr(id,1);ixAmax2=ixCmax2-kr(id,2)
      ixAmax3=ixCmax3-kr(id,3);
      ixAmin1=ixCmin1-kr(id,1);ixAmin2=ixCmin2-kr(id,2)
      ixAmin3=ixCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom_c(id))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         id)-wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom_c(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)))
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
         mom_c(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)-wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         mom_c(id))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         id)))
      ixAmax1=jxCmax1-kr(id,1);ixAmax2=jxCmax2-kr(id,2)
      ixAmax3=jxCmax3-kr(id,3);
      ixAmin1=jxCmin1-kr(id,1);ixAmin2=jxCmin2-kr(id,2)
      ixAmin3=jxCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         1)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,1),&
         0.5d0*abs(wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         mom_c(id))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         id)-wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom_c(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)))
    end do

    ! neutrals
    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmin3=ixOmin3-1;ixAmax1=ixOmax1+1
    ixAmax2=ixOmax2+1;ixAmax3=ixOmax3+1;
    do id=1,ndim
      call twofl_get_csound_prim_n(wprim,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,id,&
         tmp)
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
       2)=0.5d0*abs(wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
       mom_n(idim))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
       idim)-wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       mom_n(idim))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
       idim))

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=ixCmax1+kr(id,1);ixAmax2=ixCmax2+kr(id,2)
      ixAmax3=ixCmax3+kr(id,3);
      ixAmin1=ixCmin1+kr(id,1);ixAmin2=ixCmin2+kr(id,2)
      ixAmin3=ixCmin3+kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2),&
         0.5d0*abs(wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom_n(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)-wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom_n(id))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         id)))
      ixAmax1=ixCmax1-kr(id,1);ixAmax2=ixCmax2-kr(id,2)
      ixAmax3=ixCmax3-kr(id,3);
      ixAmin1=ixCmin1-kr(id,1);ixAmin2=ixCmin2-kr(id,2)
      ixAmin3=ixCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2),&
         0.5d0*abs(wprim(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         mom_n(id))+csound(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         id)-wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom_n(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)))
    end do

    do id=1,ndim
      if(id==idim) cycle
      ixAmax1=jxCmax1+kr(id,1);ixAmax2=jxCmax2+kr(id,2)
      ixAmax3=jxCmax3+kr(id,3);
      ixAmin1=jxCmin1+kr(id,1);ixAmin2=jxCmin2+kr(id,2)
      ixAmin3=jxCmin3+kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2),&
         0.5d0*abs(wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom_n(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)-wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         mom_n(id))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         id)))
      ixAmax1=jxCmax1-kr(id,1);ixAmax2=jxCmax2-kr(id,2)
      ixAmax3=jxCmax3-kr(id,3);
      ixAmin1=jxCmin1-kr(id,1);ixAmin2=jxCmin2-kr(id,2)
      ixAmin3=jxCmin3-kr(id,3);
      Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         2)=max(Hspeed(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,2),&
         0.5d0*abs(wprim(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         mom_n(id))+csound(jxCmin1:jxCmax1,jxCmin2:jxCmax2,jxCmin3:jxCmax3,&
         id)-wprim(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         mom_n(id))+csound(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
         id)))
    end do

  end subroutine twofl_get_H_speed_species

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine twofl_get_cbounds_one(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    use mod_constrained_transport

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
       ixImin3:ixImax3,number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    double precision :: wmean(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       nw)
    double precision :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix1,ix2,ix3

    select case (boundspeed) 
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      call get_rhoc_tot(wLP,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
      call get_rhon_tot(wLP,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(abs(rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)  +rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      call get_rhoc_tot(wRP,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
      call get_rhon_tot(wRP,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(abs(rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) +rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(0.5*(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_n(idim))+wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_c(idim)))*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + 0.5*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_n(idim))+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_c(idim)))*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call twofl_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundL)
      call twofl_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
         ixOmin3:ixOmax3)**2*(0.5*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_n(idim))+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_c(idim)))- 0.5*(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_n(idim))+wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_c(idim))))**2
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
    ! typeboundspeed=='cmaxmean'
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux))
      call get_rhon_tot(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_n(idim))/rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call get_rhoc_tot(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_c(idim))/rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call twofl_get_csound(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(max(abs(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)), abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)) ) +csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(min(abs(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)),  abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)) ) -csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
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
           1)= max(abs(tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)),&
           abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call twofl_get_csound(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundL)
      call twofl_get_csound(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,csoundR)
      csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=max(csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(0.5*(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim))+ wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim))),0.5*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_c(idim))+ wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_n(idim))))-csoundL(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(0.5*(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim))+ wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim))),0.5*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_c(idim))+ wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_n(idim))))+csoundL(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)
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
           1)=max(0.5*(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim))+ wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim))),0.5*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_c(idim))+ wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_n(idim))))+csoundL(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if
    end select

  end subroutine twofl_get_cbounds_one

  !> Calculate fast magnetosonic wave speed
  subroutine twofl_get_csound_prim_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     csound)
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
    double precision :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    integer :: ix1,ix2


    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.d0/rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if(phys_energy) then
      call twofl_get_pthermal_c_primitive(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         csound)
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_gamma*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      call twofl_get_csound2_adiab_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         csound)
    endif

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)        = &
       twofl_mag_en_all(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2-4.0d0*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * twofl_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       idim)**2 * inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    if (.not. twofl_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = sqrt(half*(cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),dxlevel(3),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3))),&
           twofl_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))*inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*kmax)
    end if

  end subroutine twofl_get_csound_prim_c

  !> Calculate fast magnetosonic wave speed
  subroutine twofl_get_csound_prim_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    if(phys_energy) then
      call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
      call twofl_get_pthermal_n_primitive(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         csound)
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_gamma*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      call twofl_get_csound2_adiab_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         csound)
    endif
    csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = sqrt(csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

  end subroutine twofl_get_csound_prim_n

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine twofl_get_cbounds_species(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,&
     ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3,idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_variables

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
    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3
    integer :: ix1,ix2,ix3

    select case (boundspeed) 
    case (1)
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      ! charges
      call get_rhoc_tot(wLP,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(abs(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      call get_rhoc_tot(wRP,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(abs(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idim))*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idim))*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call twofl_get_csound_prim_c(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundL)
      call twofl_get_csound_prim_c(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundR)


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
         ixOmin3:ixOmax3,mom_c(idim)) - wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_c(idim)))**2
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
    
      ! neutrals

      call get_rhon_tot(wLP,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(abs(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      call get_rhon_tot(wRP,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(abs(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idim))*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idim))*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call twofl_get_csound_prim_n(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundL)
      call twofl_get_csound_prim_n(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundR)


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
         ixOmin3:ixOmax3,mom_n(idim)) - wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_n(idim)))**2
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)=umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,2)=sign(one,cmin(ix1,ix2,ix3,2))*max(abs(cmin(ix1,&
               ix2,ix3,2)),Hspeed(ix1,ix2,ix3,2))
            cmax(ix1,ix2,ix3,2)=sign(one,cmax(ix1,ix2,ix3,2))*max(abs(cmax(ix1,&
               ix2,ix3,2)),Hspeed(ix1,ix2,ix3,2))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)=abs(umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if

    case (2)
    ! typeboundspeed=='cmaxmean'
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux))
     ! charges 

      call get_rhoc_tot(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_c(idim))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call twofl_get_csound_c_idim(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
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
      !neutrals
      
      call get_rhon_tot(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_n(idim))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      call twofl_get_csound_n(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)=max(abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)=min(abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),zero)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,2)=sign(one,cmin(ix1,ix2,ix3,2))*max(abs(cmin(ix1,&
               ix2,ix3,2)),Hspeed(ix1,ix2,ix3,2))
            cmax(ix1,ix2,ix3,2)=sign(one,cmax(ix1,ix2,ix3,2))*max(abs(cmax(ix1,&
               ix2,ix3,2)),Hspeed(ix1,ix2,ix3,2))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)= abs(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
    case (3)
      ! Miyoshi 2005 JCP 208, 315 equation (67)
      call twofl_get_csound_c_idim(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundL)
      call twofl_get_csound_c_idim(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
         csoundR)
      csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=max(csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim)))-csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
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
           mom_c(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if
      call twofl_get_csound_n(wLp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csoundL)
      call twofl_get_csound_n(wRp,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csoundR)
      csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=max(csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3),csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)=min(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim)))-csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        if(H_correction) then
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
            cmin(ix1,ix2,ix3,2)=sign(one,cmin(ix1,ix2,ix3,2))*max(abs(cmin(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,2))
            cmax(ix1,ix2,ix3,2)=sign(one,cmax(ix1,ix2,ix3,2))*max(abs(cmax(ix1,&
               ix2,ix3,1)),Hspeed(ix1,ix2,ix3,2))
          end do
          end do
          end do
        end if
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)=max(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim)),wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim)))+csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      end if

    end select

  end subroutine twofl_get_cbounds_species

  !> prepare velocities for ct methods
  subroutine twofl_get_ct_velocity(vcts,wLp,wRp,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     idim,cmax,cmin)
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
         mom_c(idim))+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idim)))
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
         1)=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_c(idimN))
      vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         1)=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_c(idimN))
      vcts%vbarC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         1)=(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,1) +vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim,1))/(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)+vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim))

      vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         2)=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_c(idimE))
      vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,&
         2)=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_c(idimE))
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

  end subroutine twofl_get_ct_velocity

  subroutine twofl_get_csound_c_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
     csound)
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
#if (!defined(ONE_FLUID) || ONE_FLUID==0) && (defined(A_TOT) && A_TOT == 1)
    double precision :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
#endif
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
#if (!defined(ONE_FLUID) || ONE_FLUID==0) && (defined(A_TOT) && A_TOT == 1)
    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/(rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)) 
#else
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.d0/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
#endif

    call twofl_get_csound2_c_from_conserved(w,x,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       csound)

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = twofl_mag_en_all(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)

    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2-4.0d0*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * twofl_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       idim)**2 * inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    if (.not. twofl_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = sqrt(half*(cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),dxlevel(3),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3))),&
           twofl_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))*inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*kmax)
    end if

  end subroutine twofl_get_csound_c_idim

  !> Calculate fast magnetosonic wave speed when cbounds_species=false
  subroutine twofl_get_csound_prim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
    double precision :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
#if (defined(A_TOT) && A_TOT == 1)
    double precision :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
#endif
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
#if  (defined(A_TOT) && A_TOT == 1)
    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/(rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) 
#else
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.d0/rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
#endif

    call twofl_get_csound2_primitive(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound)

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)        = &
       twofl_mag_en_all(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2-4.0d0*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * twofl_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       idim)**2 * inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    if (.not. twofl_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = sqrt(half*(cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),dxlevel(3),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3))),&
           twofl_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))*inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*kmax)
    end if

    contains
    !TODO copy it inside
    subroutine twofl_get_csound2_primitive(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       csound2)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
      double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,nw)
      double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1:ndim)
      double precision, intent(out)   :: csound2(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3)
      double precision  :: pth_c(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
      double precision  :: pth_n(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
  
      if(phys_energy) then
        call twofl_get_pthermal_c_primitive(w,x,ixImin1,ixImin2,ixImin3,&
           ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3,pth_c)
        call twofl_get_pthermal_n_primitive(w,x,ixImin1,ixImin2,ixImin3,&
           ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3,pth_n)
        call twofl_get_csound2_from_pthermal(w,x,ixImin1,ixImin2,ixImin3,&
           ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3,pth_c,pth_n,csound2)
      else
        call twofl_get_csound2_adiab(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
           csound2)
      endif
    end subroutine twofl_get_csound2_primitive

  end subroutine twofl_get_csound_prim

  subroutine twofl_get_csound2(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
    double precision  :: pth_c(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision  :: pth_n(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(phys_energy) then
      call twofl_get_pthermal_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth_c)
      call twofl_get_pthermal_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth_n)
      call twofl_get_csound2_from_pthermal(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth_c,&
         pth_n,csound2)
    else
      call twofl_get_csound2_adiab(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
    endif
  end subroutine twofl_get_csound2

  subroutine twofl_get_csound2_adiab(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision  :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision  :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=twofl_gamma*twofl_adiab*max((rhoc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)**twofl_gamma + rhon(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)**twofl_gamma)/(rhoc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)+ rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)),rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**gamma_1,rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**gamma_1) 
  end subroutine twofl_get_csound2_adiab

  subroutine twofl_get_csound(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
    double precision :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
#if (defined(A_TOT) && A_TOT == 1)
    double precision :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
#endif
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
#if (defined(A_TOT) && A_TOT == 1)
    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = 1d0/(rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) 
#else
    inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.d0/rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
#endif

    call twofl_get_csound2(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound)

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = twofl_mag_en_all(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)

    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   = b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2-4.0d0*csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * twofl_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       idim)**2 * inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    if (.not. twofl_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = sqrt(half*(cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),dxlevel(3),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)+AvMinCs2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3))),&
           twofl_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))*inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*kmax)
    end if

  end subroutine twofl_get_csound

  subroutine twofl_get_csound2_from_pthermal(w,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     pth_c,pth_n,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: pth_c(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(in)    :: pth_n(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision  :: csound1(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
#if !defined(C_TOT) || C_TOT == 0
    csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=twofl_gamma*max((pth_c(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + pth_n(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))/(rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)),pth_n(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        pth_c(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
#else
    csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=twofl_gamma*(csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + csound1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))/(rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

#endif
  end subroutine twofl_get_csound2_from_pthermal

! end cbounds_species=false

  subroutine twofl_get_csound_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: pe_n1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    call twofl_get_csound2_n_from_conserved(w,x,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       csound)
    csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = sqrt(csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
  end subroutine twofl_get_csound_n

  !> separate routines so that it is faster
  !> Calculate temperature=p/rho when in e_ the internal energy is stored
  subroutine twofl_get_temperature_from_eint_n(w, x, ixImin1,ixImin2,ixImin3,&
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

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 1d0/Rn * gamma_1 * &
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        e_n_) /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_n_)

  end subroutine twofl_get_temperature_from_eint_n

  subroutine twofl_get_temperature_from_eint_n_with_equi(w, x, ixImin1,ixImin2,&
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

      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/Rn * (gamma_1 * w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          e_n_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_pe_n0_,b0i)) /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_n_) +block%equi_vars(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_n0_,b0i))
  end subroutine twofl_get_temperature_from_eint_n_with_equi

!  subroutine twofl_get_temperature_n_pert_from_tot(Te, ixI^L, ixO^L, res)
!    use mod_global_parameters
!    integer, intent(in)          :: ixI^L, ixO^L
!    double precision, intent(in) :: Te(ixI^S)
!    double precision, intent(out):: res(ixI^S)
!      res(ixO^S) = Te(ixO^S) -1d0/Rn * &
!                block%equi_vars(ixO^S,equi_pe_n0_,0)/block%equi_vars(ixO^S,equi_rho_n0_,0)
!  end subroutine twofl_get_temperature_n_pert_from_tot

  subroutine twofl_get_temperature_n_equi(w,x, ixImin1,ixImin2,ixImin3,ixImax1,&
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
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/Rn * block%equi_vars(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_n0_,&
         b0i)/block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         equi_rho_n0_,b0i)
  end subroutine twofl_get_temperature_n_equi

  subroutine twofl_get_rho_n_equi(w, x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho_n0_,b0i)
  end subroutine twofl_get_rho_n_equi

  subroutine twofl_get_pe_n_equi(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_pe_n0_,b0i)
  end subroutine twofl_get_pe_n_equi

  !> Calculate temperature=p/rho when in e_ the total energy is stored
  !> this does not check the values of twofl_energy and twofl_internal_e, 
  !>  twofl_energy = .true. and twofl_internal_e = .false.
  !> also check small_values is avoided
  subroutine twofl_get_temperature_from_etot_n(w, x, ixImin1,ixImin2,ixImin3,&
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
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1d0/Rn * &
       (gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       e_n_)- twofl_kin_en_n(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)))/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_n_)
  end subroutine twofl_get_temperature_from_etot_n

  subroutine twofl_get_temperature_from_etot_n_with_equi(w, x, ixImin1,ixImin2,&
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
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1d0/Rn * &
       (gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       e_n_)- twofl_kin_en_n(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3)) +  block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,equi_pe_n0_,b0i))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_n_) +block%equi_vars(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_n0_,b0i))
            
  end subroutine twofl_get_temperature_from_etot_n_with_equi

  !> separate routines so that it is faster
  !> Calculate temperature=p/rho when in e_ the internal energy is stored
  subroutine twofl_get_temperature_from_eint_c(w, x, ixImin1,ixImin2,ixImin3,&
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

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 1d0/Rc * gamma_1 * &
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        e_c_) /w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_c_)

  end subroutine twofl_get_temperature_from_eint_c

  subroutine twofl_get_temperature_from_eint_c_with_equi(w, x, ixImin1,ixImin2,&
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
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/Rc * (gamma_1 * w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          e_c_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_pe_c0_,b0i)) /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_c_) +block%equi_vars(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_c0_,b0i))
  end subroutine twofl_get_temperature_from_eint_c_with_equi

!  subroutine twofl_get_temperature_c_pert_from_tot(Te, ixI^L, ixO^L, res)
!    use mod_global_parameters
!    integer, intent(in)          :: ixI^L, ixO^L
!    double precision, intent(in) :: Te(ixI^S)
!    double precision, intent(out):: res(ixI^S)
!      res(ixO^S) = Te(ixO^S) -1d0/Rc * &
!                block%equi_vars(ixO^S,equi_pe_c0_,0)/block%equi_vars(ixO^S,equi_rho_c0_,0)
!  end subroutine twofl_get_temperature_c_pert_from_tot

  subroutine twofl_get_temperature_c_equi(w,x, ixImin1,ixImin2,ixImin3,ixImax1,&
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
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0/Rc * block%equi_vars(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_c0_,&
         b0i)/block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         equi_rho_c0_,b0i)
  end subroutine twofl_get_temperature_c_equi

  subroutine twofl_get_rho_c_equi(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
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
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho_c0_,b0i)
  end subroutine twofl_get_rho_c_equi

  subroutine twofl_get_pe_c_equi(w,x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_pe_c0_,b0i)
  end subroutine twofl_get_pe_c_equi

  !> Calculate temperature=p/rho when in e_ the total energy is stored
  !> this does not check the values of twofl_energy and twofl_internal_e, 
  !>  twofl_energy = .true. and twofl_internal_e = .false.
  !> also check small_values is avoided
  subroutine twofl_get_temperature_from_etot_c(w, x, ixImin1,ixImin2,ixImin3,&
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
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1d0/Rc * &
       (gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       e_c_)- twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)- twofl_mag_en(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_c_)
  end subroutine twofl_get_temperature_from_etot_c
  subroutine twofl_get_temperature_from_eki_c(w, x, ixImin1,ixImin2,ixImin3,&
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
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1d0/Rc * &
       (gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       e_c_)- twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)))/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_c_)
  end subroutine twofl_get_temperature_from_eki_c

  subroutine twofl_get_temperature_from_etot_c_with_equi(w, x, ixImin1,ixImin2,&
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
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1d0/Rc * &
       (gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       e_c_)- twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)- twofl_mag_en(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)) +  block%equi_vars(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_c0_,b0i))/(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_c_) &
       +block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       equi_rho_c0_,b0i))
            
  end subroutine twofl_get_temperature_from_etot_c_with_equi

  subroutine twofl_get_temperature_from_eki_c_with_equi(w, x, ixImin1,ixImin2,&
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
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1d0/Rc * &
       (gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       e_c_)- twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3)) +  block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,equi_pe_c0_,b0i))/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_c_) +block%equi_vars(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_c0_,b0i))
            
  end subroutine twofl_get_temperature_from_eki_c_with_equi

  subroutine twofl_get_csound2_adiab_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision  :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=twofl_gamma*twofl_adiab*rhon(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)**gamma_1

  end subroutine twofl_get_csound2_adiab_n

  subroutine twofl_get_csound2_n_from_conserved(w,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(phys_energy) then
      call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
      call twofl_get_pthermal_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_gamma*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      call twofl_get_csound2_adiab_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         csound2)
    endif
  end subroutine twofl_get_csound2_n_from_conserved

  !! TO DELETE
  subroutine twofl_get_csound2_n_from_primitive(w,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(phys_energy) then
      call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
      call twofl_get_pthermal_n_primitive(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         csound2)
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_gamma*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      call twofl_get_csound2_adiab_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         csound2)
    endif
  end subroutine twofl_get_csound2_n_from_primitive

  subroutine twofl_get_csound2_adiab_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision  :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
    csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=twofl_gamma*twofl_adiab* rhoc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)**gamma_1

  end subroutine twofl_get_csound2_adiab_c

  subroutine twofl_get_csound2_c_from_conserved(w,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(phys_energy) then
      call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
      call twofl_get_pthermal_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_gamma*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      call twofl_get_csound2_adiab_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         csound2)
    endif
  end subroutine twofl_get_csound2_c_from_conserved

  !> Calculate fluxes within ixO^L.
  subroutine twofl_get_flux(wC,w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
       ixOmin3:ixOmax3),tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, allocatable:: vHall(:,:,:,:)
    integer                      :: idirmin, iw, idir, jdir, kdir

    ! value at the interfaces, idim =  block%iw0 --> b0i 
    ! reuse tmp, used afterwards
    ! value at the interface so we can't put momentum
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
    ! Get flux of density
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom_c(idim))*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    ! pgas is time dependent only
    if(phys_energy) then
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_c_)
    else
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_adiab*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**twofl_gamma
      if(has_equi_pe_c0) then
        pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_c0_,b0i)
      end if
    end if

    if (twofl_Hall) then
      allocate(vHall(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir))
      call twofl_getv_Hall(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vHall)
    end if

    if(B0field) tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sum(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,:,idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mag(:)),dim=ndim+1)

    ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + 0.5d0*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3, mag(:))**2, dim=ndim+1)

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    do idir=1,ndir
      if(idim==idir) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir))=ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))
        if(B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir))=f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir))= -w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))
      end if
      if (B0field) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir))=f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idim,idim)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir,idim)
      end if
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir))=f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir))+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir))
    end do

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(phys_energy) then
      if (phys_internal_e) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom_c(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_c_)
      else if(twofl_eq_energy == EQ_ENERGY_KI) then

        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim))*(wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)+pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim))*(wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_c(:)),&
           dim=ndim+1)

        if (B0field) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              e_c_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              e_c_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mom_c(idim)) * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) - sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom_c(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(:)),dim=ndim+1) * block%B0(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,idim)
        end if

        if (twofl_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
           if (twofl_etah>zero) then
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_c_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 e_c_) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 idim) * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  mag(:))**2,dim=ndim+1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mag(idim)) * sum(vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)),dim=ndim+1)
              if (B0field) then
                 f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                    e_c_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                    e_c_) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    ixOmin3:ixOmax3,idim) * tmp(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2,ixOmin3:ixOmax3) - &
                    sum(vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                    :)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                    mag(:)),dim=ndim+1) * block%B0(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,idim)
              end if
           end if
        end if
      end if !total_energy
      ! add flux of equilibrium internal energy corresponding to pe_c0
      if(has_equi_pe_c0) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=  f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim)) * block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_c0_,idim) * inv_gamma_1
#else
        if(phys_internal_e) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=  f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom_c(idim)) * block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe_c0_,idim) * inv_gamma_1
        else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=  f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom_c(idim)) * block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe_c0_,idim) * twofl_gamma * inv_gamma_1
        end if
#endif
      end if
    end if !phys_energy

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (twofl_glm) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              psi_)
        else
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idir))

        if (B0field) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))=f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(idir))+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom_c(idim))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir,idim)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mom_c(idir))*block%B0(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim,idim)
        end if

        if (twofl_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (twofl_etah>zero) then
            if (B0field) then
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mag(idir)) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mag(idir)) - vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)*(w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mag(idim))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,idim,idim)) + vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)*(w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mag(idir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,idir,idim))
            else
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mag(idir)) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mag(idir)) - vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir)*w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mag(idim)) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mag(idir))
            end if
          end if
        end if

      end if
    end do

    if (twofl_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)  = cmax_global**2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(idim))
    end if

    if (twofl_Hall) then
      deallocate(vHall)
    end if

    !!neutrals
    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_n_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom_n(idim))*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    if(phys_energy) then
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          e_n_)
    else
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_adiab*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**twofl_gamma
      if(has_equi_pe_n0) then
        pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_n0_,b0i)
      end if
    end if
    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
        !if(idim==idir) then
        !  f(ixO^S,mom_c(idir)) = pgas(ixO^S)
        !else
        !  f(ixO^S,mom_c(idir)) = 0.0d0
        !end if
        !f(ixO^S,mom_c(idir))=f(ixO^S,mom_c(idir))+w(ixO^S,mom_c(idim))*wC(ixO^S,mom_c(idir))
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom_n(idim)) * wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idir))
    end do

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        mom_n(idim)) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        mom_n(idim)) + pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    if(phys_energy) then
      !reuse pgas for storing a in the term: div (u_n * a) and make multiplication at the end
      pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_)
      if(.not. phys_internal_e) then
        ! add pressure perturbation
        pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_n_)
      end if
      ! add flux of equilibrium internal energy corresponding to pe_n0
      if(has_equi_pe_n0) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
        pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_n0_,idim) * inv_gamma_1
#else
        pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_n0_,idim) * twofl_gamma * inv_gamma_1
#endif
      end if
      ! add u_n * a in the flux
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          e_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idim)) * pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

      ! Viscosity fluxes - viscInDiv
      !if (hd_viscosity) then
      !  call visc_get_flux_prim(w, x, ixI^L, ixO^L, idim, f, phys_energy)
      !endif
    end if

  end subroutine twofl_get_flux

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine twofl_add_source(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
     wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    !use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    if (.not. qsourcesplit) then
      ! Source for solving internal energy
      if(phys_internal_e) then
        active = .true.
        call internal_energy_add_source_n(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
           w,x)
        call internal_energy_add_source_c(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
           w,x,e_c_)
      else 
#if !defined(E_RM_W0) || E_RM_W0==1
        ! add -p0 div v source terms when equi are present
        if(has_equi_pe_n0) then
          active = .true.
          call add_pe_n0_divv(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
        endif
        if(has_equi_pe_c0) then
          active = .true.
          call add_pe_c0_divv(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
        endif
#endif
        if(twofl_eq_energy == EQ_ENERGY_KI) then
          active = .true.
          call add_source_lorentz_work(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
             ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,&
             wCT,x)
        endif
      endif

      ! Source for B0 splitting
      if (B0field) then
        active = .true.
        call add_source_B0split(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(twofl_eta)>smalldouble)then
        active = .true.
        call add_source_res2(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      end if

      if (twofl_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
      end if
      !it is not added in a split manner
      if(.not. use_imex_scheme .and. has_collisions()) then
        active = .true.
        call  twofl_explicit_coll_terms_update(qdt,ixImin1,ixImin2,ixImin3,&
           ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3,w,wCT,x)
      endif
    
      if(twofl_hyperdiffusivity) then
        active = .true.
        call  add_source_hyperdiffusive(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,&
           wCT,x)
      endif  

    end if

    
    if(source_split_divb .eqv. qsourcesplit) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
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
      case (divb_linde)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
      case (divb_ct)
        continue ! Do nothing
      case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
   

    if(twofl_radiative_cooling_c) then
      call radiative_cooling_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
         wCTprim,w,x,qsourcesplit,active,rc_fl_c)
    end if
    if(twofl_radiative_cooling_n) then
      call radiative_cooling_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
         wCTprim,w,x,qsourcesplit,active,rc_fl_n)
    end if
!
!    if(twofl_viscosity) then
!      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,&
!           w,x,phys_energy,qsourcesplit,active)
!    end if
!
    if(twofl_gravity) then
      call gravity_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
         twofl_eq_energy .eq. EQ_ENERGY_KI .or. phys_total_energy,qsourcesplit,&
         active)
    end if

  end subroutine twofl_add_source

  subroutine add_pe_n0_divv(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
    double precision                :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)

    call twofl_get_v_n(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,v)
    call add_geom_PdivV(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v,&
       -block%equi_vars(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       equi_pe_n0_,0),w,x,e_n_)

  end subroutine add_pe_n0_divv

  subroutine add_pe_c0_divv(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
    double precision                :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)

    call twofl_get_v_c(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,v)
    call add_geom_PdivV(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v,&
       -block%equi_vars(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       equi_pe_c0_,0),w,x,e_c_)

  end subroutine add_pe_c0_divv

  subroutine add_geom_PdivV(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v,p,w,x,ind)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ind
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: p(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: divv(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

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
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ind)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,ind)+qdt*p(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
  end subroutine add_geom_PdivV

  !> Compute the Lorentz force (JxB)
  subroutine get_lorentz(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
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
       3), tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,3)
    integer                         :: idir, idirmin
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3)

    b=0.0d0
    do idir = 1, ndir
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir) = twofl_mag_i_all(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
  end subroutine get_lorentz

  subroutine  add_source_lorentz_work(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,wCT,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: a(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,3), b(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    
    call get_lorentz(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,a)
    call twofl_get_v_c(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,b)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_c_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_c_)+qdt*sum(a(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)*b(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir),dim=ndim+1)

  end subroutine  add_source_lorentz_work

  !> Calculate v_n vector
  subroutine twofl_get_v_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndir)
    double precision              :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer :: idir

    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)

    do idir=1,ndir
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom_n(idir)) / rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

  end subroutine twofl_get_v_n

  subroutine get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    if(has_equi_rho_n0) then
      rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_n_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho_n0_,b0i)
    else  
      rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_n_) 
    endif

  end subroutine get_rhon_tot

  subroutine twofl_get_pthermal_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    integer :: ix1,ix2,ix3, iw

    if(phys_energy) then
      if(phys_internal_e) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=gamma_1*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_n_)
      else
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_n_)- twofl_kin_en_n(w,ixImin1,ixImin2,ixImin3,&
           ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3))
      end if
      if(has_equi_pe_n0) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_n0_,b0i)
      endif
    else
      call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_adiab*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**twofl_gamma
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
              " encountered when call twofl_get_pthermal_n"
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

  end subroutine twofl_get_pthermal_n

  subroutine twofl_get_pthermal_n_primitive(w,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     pth)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(phys_energy) then
      if(has_equi_pe_n0) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_n_) + block%equi_vars(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_n0_,b0i)
      else
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_n_) 
      endif
    else
      call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_adiab*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**twofl_gamma
    end if
  end subroutine twofl_get_pthermal_n_primitive

  !> Calculate v component
  subroutine twofl_get_v_n_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision              :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom_n(idim)) / rhon(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine twofl_get_v_n_idim

  subroutine internal_energy_add_source_n(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
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
    double precision                :: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir),divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call twofl_get_pthermal_n(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    call twofl_get_v_n(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,v)
    call add_geom_PdivV(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v,-pth,w,x,e_n_)

    if(fix_small_values .and. .not. has_equi_pe_n0) then
      call twofl_handle_small_ei_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,e_n_,&
         'internal_energy_add_source')
    end if
  end subroutine internal_energy_add_source_n

  !> Calculate v_c vector
  subroutine twofl_get_v_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndir)
    double precision              :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer :: idir

    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
    do idir=1,ndir
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom_c(idir)) / rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

  end subroutine twofl_get_v_c

  subroutine get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    if(has_equi_rho_c0) then
      rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_c_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho_c0_,b0i)
    else  
      rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_c_) 
    endif

  end subroutine get_rhoc_tot

  subroutine twofl_get_pthermal_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer :: ix1,ix2,ix3, iw

    if(phys_energy) then
      if(phys_internal_e) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=gamma_1*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_c_)
      elseif(phys_total_energy) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_c_)- twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,&
           ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3)- twofl_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3))
      else
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=gamma_1*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_c_)- twofl_kin_en_c(w,ixImin1,ixImin2,ixImin3,&
           ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
           ixOmax3))
      end if
      if(has_equi_pe_c0) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_c0_,b0i)
      endif
    else
      call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_adiab*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**twofl_gamma
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
              " encountered when call twofl_get_pe_c1"
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

  end subroutine twofl_get_pthermal_c

  subroutine twofl_get_pthermal_c_primitive(w,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     pth)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    if(phys_energy) then
      if(has_equi_pe_c0) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_c_) + block%equi_vars(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_pe_c0_,b0i)
      else
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_c_) 
      endif
    else
      call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=twofl_adiab*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**twofl_gamma
    end if
  end subroutine twofl_get_pthermal_c_primitive

  !> Calculate v_c component
  subroutine twofl_get_v_c_idim(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision              :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom_c(idim)) / rhoc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine twofl_get_v_c_idim

  subroutine internal_energy_add_source_c(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,&
     ie)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                :: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir),divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call twofl_get_pthermal_c(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pth)
    call twofl_get_v_c(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,v)
    call add_geom_PdivV(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v,-pth,w,x,ie)
    if(fix_small_values .and. .not. has_equi_pe_c0) then
      call twofl_handle_small_ei_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,ie,&
         'internal_energy_add_source')
    end if
  end subroutine internal_energy_add_source_c

  !> handle small or negative internal energy
  subroutine twofl_handle_small_ei_c(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ie,&
      subname)
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
    double precision              :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision              :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    flag=.false.
    if(has_equi_pe_c0) then
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         equi_pe_c0_,0)*inv_gamma_1<small_e)flag(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie)=.true.
    else
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)<small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)=.true.
    endif  
    if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie))) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_pe_c0) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)=small_e - block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe_c0_,0)*inv_gamma_1
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
        ! to_primitive subroutine cannot be used as this error handling
        ! is also used in TC where e_to_ei is explicitly called
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)*gamma_1
        call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)*gamma_1
        call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
        do idir = 1, ndir
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mom_n(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3, mom_n(idir))/rhon(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mom_c(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3, mom_c(idir))/rhoc(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if

  end subroutine twofl_handle_small_ei_c

  !> handle small or negative internal energy
  subroutine twofl_handle_small_ei_n(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ie,&
      subname)
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
    double precision              :: rhoc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision              :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    flag=.false.
    if(has_equi_pe_n0) then
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)+block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         equi_pe_n0_,0)*inv_gamma_1<small_e)flag(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie)=.true.
    else
      where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)<small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ie)=.true.
    endif
    if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ie))) then
      select case (small_values_method)
      case ("replace")
        if(has_equi_pe_n0) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ie)=small_e - block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,equi_pe_n0_,0)*inv_gamma_1
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
           e_n_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_)*gamma_1
        call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)*gamma_1
        call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
        do idir = 1, ndir
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mom_n(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3, mom_n(idir))/rhon(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mom_c(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3, mom_c(idir))/rhoc(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if

  end subroutine twofl_handle_small_ei_n

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
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
      axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)=axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*qdt
      ! add J0xB0 source term in momentum equations
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(1:ndir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(1:ndir))+axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:ndir)
    end if

    if(phys_total_energy) then
      a=0.d0
      ! for free-free field -(vxB0) dot J0 =0
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)=wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))
      ! store full magnetic field B0+B1 in b
      if(.not.B0field_forcefree) b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)=b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,0)
      ! store velocity in a
      do idir=1,ndir
        call twofl_get_v_c_idim(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,&
           a(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idir))
      end do
      call cross_product(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,a,b,axb)
      axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)=axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:)*qdt
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)-axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)*block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           idir)
      end do
    end if

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_B0')

  end subroutine add_source_B0split

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
    if (twofl_4th_order) then
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

    if (twofl_eta>zero)then
       eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)=twofl_eta
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
       if (twofl_4th_order) then
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
       if (twofl_eta<zero)then
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
       if (phys_total_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*Bf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idir)
       end if
    end do ! idir

    if (phys_energy) then
       ! de/dt+=eta*J**2
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)+qdt*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:)**2,dim=ndim+1)
    end if

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixImin1,&
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

    if (twofl_eta>zero)then
       eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)=twofl_eta
    else
       call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImin3,ixImax1,&
          ixImax2,ixImax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,ixAmax2,ixAmax3,&
          idirmin,x,current,eta)
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,1:ndir)=zero
    do idir=idirmin,3
       tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
          idir)=current(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3,&
          idir)*eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2,ixAmin3:ixAmax3)
    end do
    curlj=0.d0
    call curlvector(tmpvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,curlj,idirmin1,1,3)
    if(stagger_grid.and.ndim==2.and.ndir==3) then
      ! if 2.5D
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(ndir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(ndir))-qdt*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ndir)
    else
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(1:ndir))-qdt*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir)
    end if

    if(phys_energy) then
      if(phys_total_energy) then
        ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
        ! de1/dt= eta J^2 - B1 dot curl(eta J)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)+qdt*(eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,:)**2,dim=ndim+1)-sum(wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1:ndir))*curlj(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir),dim=ndim+1))
      else
        ! add eta*J**2 source term in the internal energy equation
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)+qdt*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,:)**2,dim=ndim+1)
      end if

    end if

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixImin1,&
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
       1:ndir)*twofl_eta_hyper

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

    if (phys_energy) then
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
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)*qdt
    end if

    if (fix_small_values)  call twofl_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3,'add_source_hyperres')

  end subroutine add_source_hyperres

  subroutine add_source_glm(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme
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

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb, twofl_divb_4thorder)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (twofl_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_) = abs(twofl_glm_alpha)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           psi_) = dexp(-qdt*cmax_global*twofl_glm_alpha/minval(dxlevel(:)))*w(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           psi_) = dexp(-qdt*cmax_global*twofl_glm_alpha/minval(block%ds(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:),&
           dim=ndim+1))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           psi_)
      end if
    end if

    ! gradient of Psi
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
       if (phys_total_energy) then
       ! e  = e  -qdt (b . grad(Psi))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            e_c_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(idim))*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
       end if
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir))-qdt*twofl_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixImin1,&
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
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb, twofl_divb_4thorder)

    ! calculate velocity
    call twofl_get_v_c(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)

    if (phys_total_energy) then
      ! e = e - qdt (v . b) * div b
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)-qdt*sum(v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         :)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:)),&
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
         mom_c(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir))-qdt*twofl_mag_i_all(w,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixImin1,&
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
       ixImin3:ixImax3),vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb, twofl_divb_4thorder)

    ! b = b - qdt v * div b
    do idir=1,ndir
      call twofl_get_v_c_idim(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,vel)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(idir))-qdt*vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixImin1,&
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
       ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3,divb, twofl_divb_4thorder)

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

       if (typedivbdiff=='all' .and. phys_total_energy) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
            e_c_)=w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
            e_c_)+wCT(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
            mag(idim))*graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
            ixpmin3:ixpmax3)
       end if
    end do

    if (fix_small_values) call twofl_handle_small_values(.false.,w,x,ixImin1,&
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
       ixOmin3:ixOmax3)=sqrt(twofl_mag_en_all(w,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3))
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
       ixImin3:ixImax3,7-2*ndir:3),bvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)

    idirmin0 = 7-2*ndir

    bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(1:ndir))

    call curlvector(bvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,current,idirmin,&
       idirmin0,ndir)

    if(B0field) current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin0:3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin0:3)+block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin0:3)

  end subroutine get_current

  ! copied from gravity
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine gravity_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,energy,&
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
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
    double precision       :: vel(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                         :: idim

    double precision :: gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndim)

    if(qsourcesplit .eqv. grav_split) then
      active = .true.

      if (.not. associated(usr_gravity)) then
        write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
        write(*,*) "like the phys_gravity in mod_usr_methods.t"
        call mpistop("gravity_add_source: usr_gravity not defined")
      else
        call usr_gravity(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,x,&
           gravity_field)
      end if
  
      do idim = 1, ndim
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(idim)) + qdt * gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idim) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,rho_n_)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_c(idim)) + qdt * gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idim) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,rho_c_)
        if(energy) then
#if !defined(E_RM_W0) || E_RM_W0 == 1
          call twofl_get_v_n_idim(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,&
             ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
             idim,vel)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_n_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_n_) + qdt * gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idim) * vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,rho_n_)
          call twofl_get_v_c_idim(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,&
             ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
             idim,vel)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_) + qdt * gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idim) * vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,rho_c_)
#else
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_n_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_n_) + qdt * gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idim) *  wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mom_n(idim))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             e_c_) + qdt * gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idim) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,mom_c(idim))
#endif


        end if
      end do
    end if

  end subroutine gravity_add_source

  subroutine gravity_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim), w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: dxinv(1:ndim), max_grav
    integer                         :: idim

    double precision :: gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndim)

    dxinv(1)=one/dx1;dxinv(2)=one/dx2;dxinv(3)=one/dx3;

    if(.not. associated(usr_gravity)) then
      write(*,*) "mod_usr.t: please point usr_gravity to a subroutine"
      write(*,*) "like the phys_gravity in mod_usr_methods.t"
      call mpistop("gravity_get_dt: usr_gravity not defined")
    else
      call usr_gravity(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,gravity_field)
    end if

    do idim = 1, ndim
      max_grav = maxval(abs(gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idim)))
      max_grav = max(max_grav, epsilon(1.0d0))
      dtnew = min(dtnew, 1.0d0 / sqrt(max_grav * dxinv(idim)))
    end do

  end subroutine gravity_get_dt

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine twofl_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    !use mod_viscosity, only: viscosity_get_dt
    !use mod_gravity, only: gravity_get_dt

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
    if (twofl_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/twofl_eta
    else if (twofl_eta<zero)then
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

    if(twofl_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/twofl_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1:ndim))**4/twofl_eta_hyper,dtnew)
      end if
    end if

    ! the timestep related to coll terms: 1/(rho_n rho_c alpha)
    if(dtcollpar>0d0 .and. has_collisions()) then
        call coll_get_dt(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew)
    endif

    if(twofl_radiative_cooling_c) then
      call cooling_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x,&
         rc_fl_c)
    end if
    if(twofl_radiative_cooling_n) then
      call cooling_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x,&
         rc_fl_n)
    end if
!
!    if(twofl_viscosity) then
!      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
!    end if
!
    if(twofl_gravity) then
      call gravity_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)
    end if
    if(twofl_hyperdiffusivity) then
      call hyperdiffusivity_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,&
         dx3,x)
    end if


  end subroutine twofl_get_dt

  pure function has_collisions() result(res)
    logical :: res
    res = .not. twofl_alpha_coll_constant .or. twofl_alpha_coll >0d0
  end function has_collisions

  subroutine coll_get_dt(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: dtnew

    double precision :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        alpha(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, allocatable :: gamma_rec(:,:,:), gamma_ion(:,:,:)
    double precision :: max_coll_rate

    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)

    call get_alpha_coll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, alpha)
    max_coll_rate = maxval(alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * max(rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))

    if(twofl_coll_inc_ionrec) then
       allocate(gamma_ion(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
           gamma_rec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)) 
       call get_gamma_ion_rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, gamma_rec,&
           gamma_ion)
       max_coll_rate=max(max_coll_rate, maxval(gamma_ion(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)), maxval(gamma_rec(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3))) 
       deallocate(gamma_ion, gamma_rec) 
    endif
    dtnew = min(dtcollpar/max_coll_rate, dtnew)

  end subroutine coll_get_dt

  ! Add geometrical source terms to w
  subroutine twofl_add_source_geom(qdt,dtfactor,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     wCT,w,x)
    use mod_global_parameters
    use mod_geometry

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
       rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    ! charges

    mr_=mom_c(1); mphi_=mom_c(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_
    call get_rhoc_tot(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)

    select case (coordinate)
    case (cylindrical)
      call twofl_get_p_c_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)

      if(phi_>0) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)*(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)**2+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mphi_)**2/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)*(-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,br_))
        if(.not.stagger_grid) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             1)*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mr_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mphi_)) /rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end if
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if
      if(twofl_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         br_)+qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         psi_)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
    case (spherical)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmin3=ixOmin3-kr(1,3);h1xmax1=ixOmax1-kr(1,1)
       h1xmax2=ixOmax2-kr(1,2);h1xmax3=ixOmax3-kr(1,3)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmin3=ixOmin3-kr(2,3);h2xmax1=ixOmax1-kr(2,1)
       h2xmax2=ixOmax2-kr(2,2);h2xmax3=ixOmax3-kr(2,3);
       call twofl_get_p_c_total(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
              ixOmin3:ixOmax3,mom_c(idir))**2/rho(ixOmin1:ixOmax1,&
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
          mom_c(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom_c(1))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)
       ! b1
       if(twofl_glm) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(1))+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1)*2.0d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,psi_)
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
       ! This will make hydrostatic p=const an exact solution
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom_c(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom_c(2))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,2)-block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          h2xmin3:h2xmax3,2)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom_c(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,mom_c(2))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
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
            ixOmin3:ixOmax3,mom_c(3))**2/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
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
          mom_c(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom_c(2))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom_c(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom_c(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(1)))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
         if(B0field) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom_c(1))*block%B0(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,0) -wCT(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mom_c(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,1,0))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)
         end if
         if(twofl_glm) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) + dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,psi_)
         end if
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mag(2))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1)
       end if
      

       if(ndir==3) then
         ! m3
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=-(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom_c(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom_c(1))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mag(1)))  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom_c(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom_c(3))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
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
            mom_c(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            mom_c(3))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            1)
         ! b3
         if(.not.stagger_grid) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom_c(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom_c(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(1)))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom_c(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom_c(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mag(3)))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2)) /(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2)))
           if (B0field) then
              tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mom_c(1))*block%B0(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,3,0) -wCT(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom_c(3))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,1,0))/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,mom_c(3))*block%B0(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,0) -wCT(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom_c(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,3,0))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,2)) /(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,2)))
           end if
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              mag(3))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,1)
         end if
       end if
    end select

    ! neutrals 
    !TODO no dust: see and implement them from hd/mod_hd_phys !
    !uncomment cartesian expansion
    call get_rhon_tot(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    call twofl_get_pthermal_n(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, tmp1)

    select case (coordinate)
!    case(Cartesian_expansion)
!      !the user provides the functions of exp_factor and del_exp_factor
!      if(associated(usr_set_surface)) call usr_set_surface(ixI^L,x,block%dx,exp_factor,del_exp_factor,exp_factor_primitive)
!      tmp(ixO^S) = tmp1(ixO^S)*del_exp_factor(ixO^S)/exp_factor(ixO^S)
!      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*tmp(ixO^S)

    case (cylindrical)
      mr_   = mom_n(r_)
      if (phi_ > 0) then
         where (rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) > 0d0)
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, mphi_)**2 / rho(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mr_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mr_) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, r_)
         end where
         ! s[mphi]=(-mphi*mr/rho)/radius
         where (rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) > 0d0)
             tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mphi_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mr_) / rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mphi_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mphi_) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, r_)
         end where
      else
         ! s[mr]=2pthermal/radius
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mr_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mr_) + qdt * tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, r_)
      end if
    case (spherical)
       if(phi_>0) mphi_ = mom_n(phi_)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmin3=ixOmin3-kr(1,3);h1xmax1=ixOmax1-kr(1,1)
       h1xmax2=ixOmax2-kr(1,2);h1xmax3=ixOmax3-kr(1,3)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmin3=ixOmin3-kr(2,3);h2xmax1=ixOmax1-kr(2,1)
       h2xmax2=ixOmax2-kr(2,2);h2xmax3=ixOmax3-kr(2,3);
       ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1) - block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,h1xmin3:h1xmax3,&
           1)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       if (ndir > 1) then
         do idir = 2, ndir
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3, mom_n(idir))**2 / rho(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,ixOmin3:ixOmax3)
         end do
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)

       
       ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*(mphi**2/rho+p)/r
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1) * (block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, 2) - block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          h2xmin3:h2xmax3, 2)) / block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (ndir == 3) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, mom_n(3))**2 / rho(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)) / tan(x(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, 2))
       end if
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) - (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, mom_n(2)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, mr_)) / rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(2)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom_n(2)) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)

       if (ndir == 3) then
         ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) = -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, mom_n(3)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, mr_)) / rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)- (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, mom_n(2)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, mom_n(3))) / rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) / tan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, 2))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom_n(3)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom_n(3)) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, 1)
       end if
      
    end select

    contains
      subroutine twofl_get_p_c_total(w,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,p)
        use mod_global_parameters
    
        integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
        double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,nw)
        double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndim)
        double precision, intent(out)   :: p(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)
    
        call twofl_get_pthermal_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,p)
    
        p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = p(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3) + 0.5d0 * sum(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(:))**2, dim=ndim+1)
    
      end subroutine twofl_get_p_c_total

  end subroutine twofl_add_source_geom

  subroutine twofl_get_temp_c_pert_from_etot(w, x, ixImin1,ixImin2,ixImin3,&
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

    ! store pe1 in res
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(gamma_1*(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_c_)- twofl_kin_en_c(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)- twofl_mag_en(w,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3)))
    if(has_equi_pe_c0) then
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_pe_c0_,b0i)
      if(has_equi_rho_c0) then
        res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/(Rc * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,rho_c_)+ block%equi_vars(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_c0_,&
           b0i))) - block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_c0_,b0i)/(Rc * &
           block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           equi_rho_c0_,b0i))
      else
        ! infinite equi temperature with p0 and 0 density
        res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
      endif
    else
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/(Rc * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_c_))
    endif

  end subroutine twofl_get_temp_c_pert_from_etot

  !> Compute 2 times total magnetic energy
  function twofl_mag_en_all(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if (B0field) then
      mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:,b0i))**2, dim=ndim+1)
    else
      mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mag(:))**2, dim=ndim+1)
    end if
  end function twofl_mag_en_all

  !> Compute full magnetic field by direction
  function twofl_mag_i_all(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: mgf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if (B0field) then
      mgf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(idir))+block%B0(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,b0i)
    else
      mgf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(idir))
    end if
  end function twofl_mag_i_all

  !> Compute evolving magnetic energy
  function twofl_mag_en(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.5d0 * &
       sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mag(:))**2,&
        dim=ndim+1)
  end function twofl_mag_en

  !> compute kinetic energy of neutrals
  function twofl_kin_en_n(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if(has_equi_rho_n0) then
      ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mom_n(:))**2, dim=ndim+1) / (w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          rho_n_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho_n0_,0))
    else
      ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mom_n(:))**2, dim=ndim+1) / w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, rho_n_)
    endif

  end function twofl_kin_en_n

  subroutine twofl_get_temp_n_pert_from_etot(w, x, ixImin1,ixImin2,ixImin3,&
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

    ! store pe1 in res
    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(gamma_1*(w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_c_)- twofl_kin_en_c(w,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3)))
    if(has_equi_pe_n0) then
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_pe_n0_,b0i)
      if(has_equi_rho_n0) then
        res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/(Rn * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,rho_n_)+ block%equi_vars(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,equi_rho_n0_,&
           b0i))) - block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_n0_,b0i)/(Rn * &
           block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           equi_rho_n0_,b0i))
      else
        ! infinite equi temperature with p0 and 0 density
        res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
      endif
    else
      res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/(Rn * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,rho_n_))
    endif

  end subroutine twofl_get_temp_n_pert_from_etot

  !> compute kinetic energy of charges
  !> w are conserved variables
  function twofl_kin_en_c(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) result(ke)
    use mod_global_parameters, only: nw, ndim,block
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    if(has_equi_rho_c0) then
      ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mom_c(:))**2, dim=ndim+1) / (w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          rho_c_) + block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,equi_rho_c0_,0))
    else
      ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mom_c(:))**2, dim=ndim+1) / w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, rho_c_)
    endif
  end function twofl_kin_en_c

  subroutine twofl_getv_Hall(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vHall)
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

    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:3) = zero
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idirmin:3) = - twofl_etah*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,idirmin:3)
    do idir = idirmin, 3
       vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir) = vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          idir)/rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do

  end subroutine twofl_getv_Hall

! the following not used
!  subroutine twofl_getdt_Hall(w,x,ixI^L,ixO^L,dx^D,dthall)
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
!    ! because we have that in cmax now:
!    return
!
!    ^D&dxarr(^D)=dx^D;
!
!    if (.not. B0field) then
!       bmag(ixO^S)=sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1))
!       bmag(ixO^S)=sqrt(sum((w(ixO^S,mag(:)) + block%B0(ixO^S,1:ndir,b0i))**2))
!    end if
!
!    if(slab_uniform) then
!      dthall=dtdiffpar*minval(dxarr(1:ndim))**2.0d0/(twofl_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_c_)))
!    else
!      dthall=dtdiffpar*minval(block%ds(ixO^S,1:ndim))**2.0d0/(twofl_etah*maxval(bmag(ixO^S)/w(ixO^S,rho_c_)))
!    end if
!
!  end subroutine twofl_getdt_Hall

  subroutine twofl_modify_wLR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
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

      if(phys_total_energy) then
        wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)-half*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
        wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)-half*wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
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
      if(phys_total_energy) then
        wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)+half*wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
        wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)=wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_)+half*wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(idir))**2
      end if
    end if

    if(associated(usr_set_wLR)) call usr_set_wLR(ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine twofl_modify_wLR

  subroutine twofl_boundary_adjust(igrid,psb)
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

  end subroutine twofl_boundary_adjust

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
      
     case(2)
       
       
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
      
     case(3)
       
       
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
      
     case(4)
       
       
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
      
     
     case(5)
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
     case(6)
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
    
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  
  subroutine twofl_clean_divb_multigrid(qdt, qt, active)
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
          print *, "divb_multigrid warning: unknown b.c.: ", &
               typeboundary(mag(idim), n)
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
          ixMlo3,ixMhi1,ixMhi2,ixMhi3, tmp, twofl_divb_4thorder)
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
         call twofl_face_to_center(ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
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

       if(phys_total_energy) then
         ! Determine magnetic energy difference
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
            ixMlo3:ixMhi3) = 0.5_dp * (sum(ps(igrid)%w(ixMlo1:ixMhi1,&
            ixMlo2:ixMhi2,ixMlo3:ixMhi3, mag(1:ndim))**2,&
             dim=ndim+1) - tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))
         ! Keep thermal pressure the same
         ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
             e_c_) = ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3,&
             e_c_) + tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2,ixMlo3:ixMhi3)
       end if
    end do

    active = .true.

  end subroutine twofl_clean_divb_multigrid
 

  subroutine twofl_update_faces(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,qdt,wprim,fC,&
     fE,sCT,s,vcts)
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

  end subroutine twofl_update_faces

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
       ixImin3:ixImax3,sdim:3) :: E_resi

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
    ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1;

    ! if there is resistivity, get eta J
    if(twofl_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT,s,E_resi)

    fE=zero

    do idim1=1,ndim 
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=sdim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
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
            if(twofl_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_resi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
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
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
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
      ! Time update
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
    integer                            :: hxCmin1,hxCmin2,hxCmin3,hxCmax1,&
       hxCmax2,hxCmax3,ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,jxCmin1,&
       jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3,ixAmin1,ixAmin2,ixAmin3,ixAmax1,&
       ixAmax2,ixAmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm)

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
            mom_c(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)=ECC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idir)-Btot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            idim1)*wp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
            mom_c(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(twofl_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT,s,E_resi)
    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    fE=zero
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

            ! add current component of electric field at cell edges E=-vxB+eta J
            if(twofl_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)+E_resi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idir)
            ! times time step and edge length 
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)*qdt*s%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
               idir)
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
      ! Divide by the area of the face to get dB/dt
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
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
    double precision                   :: btilL(s%ixGsmin1:s%ixGsmax1,&
       s%ixGsmin2:s%ixGsmax2,s%ixGsmin3:s%ixGsmax3,ndim)
    double precision                   :: btilR(s%ixGsmin1:s%ixGsmax1,&
       s%ixGsmin2:s%ixGsmax2,s%ixGsmin3:s%ixGsmax3,ndim)
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
    if(twofl_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3,sCT,s,E_resi)
    fE=zero

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

      ! add current component of electric field at cell edges E=-vxB+eta J
      if(twofl_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
         idir)+E_resi(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir)
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
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2)
      ixCmin3=ixOmin3-kr(idim1,3);
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
      ! Time update
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
    if(twofl_eta>zero)then
      jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         :)=jce(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:)*twofl_eta
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
      enddo
    end if

    end associate
  end subroutine get_resistive_electric_field

  !> calculate cell-center values from face-center values
  subroutine twofl_face_to_center(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
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

  end subroutine twofl_face_to_center

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

  subroutine  hyperdiffusivity_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,&
     dx2,dx3,x)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: dx1,dx2,dx3
    double precision, intent(inout) :: dtnew

    double precision :: nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       temp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    double precision :: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),csound_dim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision              :: dxarr(ndim)
    double precision :: maxCoef
    integer :: ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, hxbmin1,&
       hxbmin2,hxbmin3,hxbmax1,hxbmax2,hxbmax3, hxmin1,hxmin2,hxmin3,hxmax1,&
       hxmax2,hxmax3, ii, jj


    dxarr(1)=dx1;dxarr(2)=dx2;dxarr(3)=dx3;
    maxCoef = smalldouble

    ! charges
    call twofl_get_v_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,vel)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,rho)
    call twofl_get_csound2_c_from_conserved(w,x,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       csound)
    csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = sqrt(csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)) + sqrt(twofl_mag_en_all(w,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3) /rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
    csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) + sqrt(sum(vel(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)**2 ,dim=ndim+1))
    do ii=1,ndim
      call div_vel_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, vel, ii,&
          divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii))
      hxmin1=ixImin1+1;hxmin2=ixImin2+1;hxmin3=ixImin3+1;
      hxmax1=ixImax1-1;hxmax2=ixImax2-1;hxmax3=ixImax3-1;
      hxbmin1=hxmin1-kr(ii,1);hxbmin2=hxmin2-kr(ii,2);hxbmin3=hxmin3-kr(ii,3)
      hxbmax1=hxmax1-kr(ii,1);hxbmax2=hxmax2-kr(ii,2);hxbmax3=hxmax3-kr(ii,3);
      csound_dim(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,&
         ii) = (csound(hxbmin1:hxbmax1,hxbmin2:hxbmax2,&
         hxbmin3:hxbmax3)+csound(hxmin1:hxmax1,hxmin2:hxmax2,&
         hxmin3:hxmax3))/2d0
    enddo
    call twofl_get_temp_c_pert_from_etot(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, temp)
    do ii=1,ndim
      !TODO the following is copied
      !rho_c
      call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
         ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,rho_c_), ii, tmp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3))
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = c_hyp(rho_c_) * csound_dim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ii) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + c_shk(rho_c_) * (dxlevel(ii)**2) &
         *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ii)
      maxCoef = max(maxCoef,maxval(nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      !TH c  
      call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
         ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, temp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3), ii, tmp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3))
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = c_hyp(e_c_) * csound_dim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ii) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + c_shk(e_c_) * (dxlevel(ii)**2) &
         *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ii)
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = nu(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) * rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) * Rc/(twofl_gamma-1d0)
      maxCoef = max(maxCoef,maxval(nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      !visc c
      do jj=1,ndir
        call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
            vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,jj), ii,&
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
        nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = c_hyp(mom_c(jj)) * csound_dim(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           ii) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + c_shk(mom_c(jj)) * (dxlevel(ii)**2) &
           *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ii)
        nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) 
        maxCoef = max(maxCoef,maxval(nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)))
      enddo

      ! Ohmic
      do jj=1,ndir
        if(ii .ne. jj) then
          call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
              w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(jj)), ii,&
              tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
          nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = c_hyp(mag(jj)) * csound_dim(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             ii) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + c_shk(mag(jj)) * (dxlevel(ii)**2) &
             *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ii)
          maxCoef = max(maxCoef,maxval(nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)))
        endif
      enddo

    enddo 
  
      !TODO the following is copied, as charges, and as in add_source!
    ! neutrals
    call twofl_get_v_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,vel)
    call  twofl_get_csound_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,csound)
    csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) + sqrt(sum(vel(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)**2 ,dim=ndim+1))
    do ii=1,ndim
      call div_vel_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, vel, ii,&
          divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii))
      hxmin1=ixImin1+1;hxmin2=ixImin2+1;hxmin3=ixImin3+1;
      hxmax1=ixImax1-1;hxmax2=ixImax2-1;hxmax3=ixImax3-1;
      hxbmin1=hxmin1-kr(ii,1);hxbmin2=hxmin2-kr(ii,2);hxbmin3=hxmin3-kr(ii,3)
      hxbmax1=hxmax1-kr(ii,1);hxbmax2=hxmax2-kr(ii,2);hxbmax3=hxmax3-kr(ii,3);
      csound_dim(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,&
         ii) = (csound(hxbmin1:hxbmax1,hxbmin2:hxbmax2,&
         hxbmin3:hxbmax3)+csound(hxmin1:hxmax1,hxmin2:hxmax2,&
         hxmin3:hxmax3))/2d0
    enddo
    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    call twofl_get_temp_n_pert_from_etot(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, temp)
    do ii=1,ndim
      !rho_n
      call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
         ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,rho_n_), ii, tmp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3))
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = c_hyp(rho_n_) * csound_dim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ii) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + c_shk(rho_n_) * (dxlevel(ii)**2) &
         *divv(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,ii)
      maxCoef = max(maxCoef,maxval(nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      !TH n  
      call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
         ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, temp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3), ii, tmp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3))
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = c_hyp(e_n_) * csound_dim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ii) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + c_shk(e_n_) * (dxlevel(ii)**2) &
         *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ii)
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = nu(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) * rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) * Rn/(twofl_gamma-1d0)
      maxCoef = max(maxCoef,maxval(nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)))

      !visc n
      do jj=1,ndir
        call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
            vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,jj), ii,&
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
        nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = c_hyp(mom_n(jj)) * csound_dim(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           ii) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + c_shk(mom_n(jj)) * (dxlevel(ii)**2) &
           *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,ii)
        nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) 
        maxCoef = max(maxCoef,maxval(nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)))
      enddo
    enddo 

    dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**2/maxCoef,dtnew)
  end subroutine hyperdiffusivity_get_dt

  subroutine  add_source_hyperdiffusive(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,wCT,x)
    use mod_global_parameters
    use mod_hyperdiffusivity

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    double precision :: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),csound_dim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    integer :: ii,ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
       hxbmin1,hxbmin2,hxbmin3,hxbmax1,hxbmax2,hxbmax3,hxmin1,hxmin2,hxmin3,&
       hxmax1,hxmax2,hxmax3
    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call twofl_get_v_c(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,vel)
    call get_rhoc_tot(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,rho)
    call twofl_get_csound2_c_from_conserved(wCT,x,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       csound)
    csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = sqrt(csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)) + sqrt(twofl_mag_en_all(wCT,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3) /rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
    csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) + sqrt(sum(vel(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)**2 ,dim=ndim+1))
    do ii=1,ndim
      call div_vel_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, vel, ii,&
          divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii))
      hxmin1=ixImin1+1;hxmin2=ixImin2+1;hxmin3=ixImin3+1;
      hxmax1=ixImax1-1;hxmax2=ixImax2-1;hxmax3=ixImax3-1;
      hxbmin1=hxmin1-kr(ii,1);hxbmin2=hxmin2-kr(ii,2);hxbmin3=hxmin3-kr(ii,3)
      hxbmax1=hxmax1-kr(ii,1);hxbmax2=hxmax2-kr(ii,2);hxbmax3=hxmax3-kr(ii,3);
      csound_dim(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,&
         ii) = (csound(hxbmin1:hxbmax1,hxbmin2:hxbmax2,&
         hxbmin3:hxbmax3)+csound(hxmin1:hxmax1,hxmin2:hxmax2,&
         hxmin3:hxmax3))/2d0
    enddo
    call add_density_hyper_Source(rho_c_)
    call add_viscosity_hyper_Source(rho,mom_c(1), e_c_)
    call add_th_cond_c_hyper_Source(rho)
    call add_ohmic_hyper_Source()

    call twofl_get_v_n(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,vel)
    call  twofl_get_csound_n(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,csound)
    csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) + sqrt(sum(vel(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)**2 ,dim=ndim+1))
    do ii=1,ndim
      call div_vel_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, vel, ii,&
          divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii))
      hxmin1=ixImin1+1;hxmin2=ixImin2+1;hxmin3=ixImin3+1;
      hxmax1=ixImax1-1;hxmax2=ixImax2-1;hxmax3=ixImax3-1;
      hxbmin1=hxmin1-kr(ii,1);hxbmin2=hxmin2-kr(ii,2);hxbmin3=hxmin3-kr(ii,3)
      hxbmax1=hxmax1-kr(ii,1);hxbmax2=hxmax2-kr(ii,2);hxbmax3=hxmax3-kr(ii,3);
      csound_dim(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,&
         ii) = (csound(hxbmin1:hxbmax1,hxbmin2:hxbmax2,&
         hxbmin3:hxbmax3)+csound(hxmin1:hxmax1,hxmin2:hxmax2,&
         hxmin3:hxmax3))/2d0
    enddo
    call add_density_hyper_Source(rho_n_)
    call get_rhon_tot(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,rho)
    call add_viscosity_hyper_Source(rho,mom_n(1), e_n_)
    call add_th_cond_n_hyper_Source(rho)

    contains

    subroutine add_density_hyper_Source(index_rho)
      integer, intent(in) :: index_rho

      double precision :: nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
          tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

      do ii=1,ndim
        call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
            wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,index_rho), ii,&
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
        nu(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
           ixOOmin3:ixOOmax3) = c_hyp(index_rho) * &
           csound_dim(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,&
           ii) * dxlevel(ii) *  tmp(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
           ixOOmin3:ixOOmax3) + c_shk(index_rho) * (dxlevel(ii)**2) &
           *divv(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,ii)
        !print*, "IXOO HYP ", ixOO^L, " IDIMM ", ii
        call second_same_deriv(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
            nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
            wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,index_rho), ii,&
            tmp)

        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           index_rho) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           index_rho) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        !print*, "RHO ", index_rho, maxval(abs(tmp(ixO^S)))
      enddo
    end subroutine add_density_hyper_Source   

    subroutine add_th_cond_c_hyper_Source(var2)
      double precision, intent(in) :: var2(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
      double precision :: nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
          tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
          var(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
      call twofl_get_temp_c_pert_from_etot(wCT, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, var)
      do ii=1,ndim
        call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
            var(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3), ii,&
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
        nu(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
           ixOOmin3:ixOOmax3) = c_hyp(e_c_) * csound_dim(ixOOmin1:ixOOmax1,&
           ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,&
           ii) * dxlevel(ii) *  tmp(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
           ixOOmin3:ixOOmax3) + c_shk(e_c_) * (dxlevel(ii)**2) &
           *divv(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,ii)
        call second_same_deriv2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
            nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
            var2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) ,&
           var(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3), ii, tmp)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_c_) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * Rc/(twofl_gamma-1d0)
      !print*, "TH C ", maxval(abs(tmp(ixO^S)))
      enddo
    end subroutine add_th_cond_c_hyper_Source   

    subroutine add_th_cond_n_hyper_Source(var2)
      double precision, intent(in) :: var2(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
      double precision :: nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
          tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
          var(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
      call twofl_get_temp_n_pert_from_etot(wCT, x, ixImin1,ixImin2,ixImin3,&
         ixImax1,ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, var)
      do ii=1,ndim
        call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
            var(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3), ii,&
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
        nu(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
           ixOOmin3:ixOOmax3) = c_hyp(e_n_) * csound_dim(ixOOmin1:ixOOmax1,&
           ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,&
           ii) * dxlevel(ii) *  tmp(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
           ixOOmin3:ixOOmax3) + c_shk(e_n_) * (dxlevel(ii)**2) &
           *divv(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,ii)
        call second_same_deriv2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
            nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
            var2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) ,&
           var(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3), ii, tmp)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_n_) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * Rn/(twofl_gamma-1d0)
      !print*, "TH N ", maxval(abs(tmp(ixO^S)))
      enddo
    end subroutine add_th_cond_n_hyper_Source   

    subroutine add_viscosity_hyper_Source(rho,index_mom1, index_e)
      double precision, intent(in) :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)
      integer, intent(in) :: index_mom1, index_e

      double precision :: nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir,1:ndim), tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
      integer :: jj

      do jj=1,ndir
        do ii=1,ndim
          call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
              vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,jj), ii,&
              tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
          nu(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,jj,&
             ii) = c_hyp(index_mom1-1+jj) * csound_dim(ixOOmin1:ixOOmax1,&
             ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,&
             ii) * dxlevel(ii) *  tmp(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
             ixOOmin3:ixOOmax3) + c_shk(index_mom1-1+jj) * (dxlevel(ii)**2) &
             *divv(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,ii)
        enddo
      enddo
      
      do jj=1,ndir
        do ii=1,ndim
          call second_same_deriv2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
              nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,jj,ii),&
              rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
              vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,jj), ii,&
              tmp)
          call second_same_deriv2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
              nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,jj,ii),&
              wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
             index_mom1-1+jj), vel(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,jj), ii, tmp2)
          if(ii .eq. jj) then
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               index_mom1-1+jj) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,index_mom1-1+jj) + qdt * tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               index_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               index_e) + qdt * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)

          else
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               index_mom1-1+jj) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,index_mom1-1+jj) + 0.5*qdt * &
               tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               index_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               index_e) + 0.5*qdt * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
            call second_cross_deriv2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
                nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii,jj),&
                rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
                vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii), jj,&
                ii, tmp)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               index_mom1-1+jj) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,index_mom1-1+jj) + 0.5*qdt * &
               tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
            call second_cross_deriv2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
                nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,jj,ii),&
                wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               index_mom1-1+jj), vel(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3,jj), ii, jj, tmp2)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               index_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               index_e) + 0.5*qdt * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
          endif

        enddo
      enddo

    end subroutine add_viscosity_hyper_Source   

    subroutine add_ohmic_hyper_Source()
      double precision :: nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1:ndir,1:ndim), tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
      integer :: jj

      do jj=1,ndir
        do ii=1,ndim
          if(ii .ne. jj) then
            call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
                ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
                wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(jj)),&
                ii, tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
            nu(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,jj,&
               ii) = c_hyp(mag(jj)) * csound_dim(ixOOmin1:ixOOmax1,&
               ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,&
               ii) * dxlevel(ii) *  tmp(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
               ixOOmin3:ixOOmax3) + c_shk(mag(jj)) * (dxlevel(ii)**2) &
               *divv(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3,ii)
          endif
        enddo
      enddo
      
      do jj=1,ndir
        do ii=1,ndim
          if(ii .ne. jj) then
            !mag field
            call second_same_deriv(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
                nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,jj,ii),&
                wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(jj)),&
                ii, tmp)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(jj)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(jj)) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
            call second_cross_deriv(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
                nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii,jj),&
                 wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(ii)),&
                jj, ii, tmp)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(jj)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               mag(jj)) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
            !in the total energy
            call second_same_deriv(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
                nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,jj,ii),&
                wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(jj)),&
                ii, tmp)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               e_c_) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
            call second_cross_deriv2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
               ixImax3, ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
                nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ii,jj),&
                wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(jj)),&
                wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(ii)),&
                jj, ii, tmp)
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               e_c_) + qdt * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
          endif

        enddo
      enddo

    end subroutine add_ohmic_hyper_Source   

  end subroutine  add_source_hyperdiffusive

  function dump_hyperdiffusivity_coef_x(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
      nwc) result(wnew)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nwc
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) 
    double precision   :: wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        1:nwc)

    if(nw .ne. nwc) call mpistop("nw != nwc")
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:nw) =  dump_hyperdiffusivity_coef_dim(ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
        1) 

  end function dump_hyperdiffusivity_coef_x

  function dump_hyperdiffusivity_coef_y(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
      nwc) result(wnew)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nwc
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) 
    double precision   :: wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        1:nwc)

    if(nw .ne. nwc) call mpistop("nw != nwc")
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:nw) =  dump_hyperdiffusivity_coef_dim(ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
        2) 

  end function dump_hyperdiffusivity_coef_y

  function dump_hyperdiffusivity_coef_z(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
      nwc) result(wnew)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nwc
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) 
    double precision   :: wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        1:nwc)

    if(nw .ne. nwc) call mpistop("nw != nwc")
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:nw) =  dump_hyperdiffusivity_coef_dim(ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
        3) 

  end function dump_hyperdiffusivity_coef_z

  function dump_hyperdiffusivity_coef_dim(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOPmin1,ixOPmin2,ixOPmin3,ixOPmax1,ixOPmax2,ixOPmax3, w,&
      x, ii) result(wnew)
    use mod_global_parameters
    use mod_hyperdiffusivity
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOPmin1,ixOPmin2,ixOPmin3,ixOPmax1,ixOPmax2,ixOPmax3, ii
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) 
    double precision   :: wnew(ixOPmin1:ixOPmax1,ixOPmin2:ixOPmax2,&
       ixOPmin3:ixOPmax3, 1:nw)

    double precision :: nu(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       temp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)
    double precision :: csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),csound_dim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision              :: dxarr(ndim)
    integer :: ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, hxbmin1,&
       hxbmin2,hxbmin3,hxbmax1,hxbmax2,hxbmax3, hxmin1,hxmin2,hxmin3,hxmax1,&
       hxmax2,hxmax3,  jj, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3

    ! this is done because of save_physical_boundary = true
    ixOmin1=max(ixOPmin1,ixImin1+3);ixOmin2=max(ixOPmin2,ixImin2+3)
    ixOmin3=max(ixOPmin3,ixImin3+3);
    ixOmax1=min(ixOPmax1,ixImax1-3);ixOmax2=min(ixOPmax2,ixImax2-3)
    ixOmax3=min(ixOPmax3,ixImax3-3);

    wnew(ixOPmin1:ixOPmax1,ixOPmin2:ixOPmax2,ixOPmin3:ixOPmax3,1:nw) = 0d0

    ! charges
    call twofl_get_temp_c_pert_from_etot(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, temp)
    call twofl_get_v_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,vel)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,rho)
    call twofl_get_csound2_c_from_conserved(w,x,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       csound) 
    csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = sqrt(csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)) + sqrt(twofl_mag_en_all(w,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3) /rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
    csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) + sqrt(sum(vel(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)**2 ,dim=ndim+1))
    !for dim
    call div_vel_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, vel, ii,&
        divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
    hxmin1=ixImin1+1;hxmin2=ixImin2+1;hxmin3=ixImin3+1;
    hxmax1=ixImax1-1;hxmax2=ixImax2-1;hxmax3=ixImax3-1;
    hxbmin1=hxmin1-kr(ii,1);hxbmin2=hxmin2-kr(ii,2);hxbmin3=hxmin3-kr(ii,3)
    hxbmax1=hxmax1-kr(ii,1);hxbmax2=hxmax2-kr(ii,2);hxbmax3=hxmax3-kr(ii,3);
    csound_dim(hxmin1:hxmax1,hxmin2:hxmax2,&
       hxmin3:hxmax3) = (csound(hxbmin1:hxbmax1,hxbmin2:hxbmax2,&
       hxbmin3:hxbmax3)+csound(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3))/2d0

    !TODO the following is copied
    !rho_c
    call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
       ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,rho_c_), ii, tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3))
    nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = c_hyp(rho_c_) * &
       csound_dim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + c_shk(rho_c_) * (dxlevel(ii)**2) &
       *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_c_) = nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    !TH c  
    call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
       ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, temp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3), ii, tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3))
    nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = c_hyp(e_c_) * &
       csound_dim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + c_shk(e_c_) * (dxlevel(ii)**2) *divv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = nu(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * Rc/(twofl_gamma-1d0)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       e_c_) = nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    !visc c
    do jj=1,ndir
      call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
         ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, vel(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,jj), ii, tmp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3))
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = c_hyp(mom_c(jj)) * csound_dim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) + c_shk(mom_c(jj)) * &
         (dxlevel(ii)**2) *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = nu(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) * rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(jj)) = nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    enddo

    ! Ohmic
    do jj=1,ndir
      if(ii .ne. jj) then
        call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3,&
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(jj)), ii,&
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
        nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = c_hyp(mag(jj)) * csound_dim(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3) * dxlevel(ii) *  &
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + c_shk(mag(jj)) * (dxlevel(ii)**2) &
           *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(jj)) = nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      endif
    enddo

    !end for dim
  
    ! neutrals
    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    call twofl_get_temp_n_pert_from_etot(w, x, ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, temp)
    call twofl_get_v_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,vel)
    call  twofl_get_csound_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,csound)
    csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = csound(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) + sqrt(sum(vel(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)**2 ,dim=ndim+1))
    !for dim
    call div_vel_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOOmin1,ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, vel, ii,&
        divv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
    hxbmin1=ixOOmin1-kr(ii,1);hxbmin2=ixOOmin2-kr(ii,2)
    hxbmin3=ixOOmin3-kr(ii,3);hxbmax1=ixOOmax1-kr(ii,1)
    hxbmax2=ixOOmax2-kr(ii,2);hxbmax3=ixOOmax3-kr(ii,3);
    csound_dim(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
       ixOOmin3:ixOOmax3) = (csound(hxbmin1:hxbmax1,hxbmin2:hxbmax2,&
       hxbmin3:hxbmax3)+csound(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,&
       ixOOmin3:ixOOmax3))/2d0
    !rho_n
    call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
       ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,rho_n_), ii, tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3))
    nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = c_hyp(rho_n_) * &
       csound_dim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + c_shk(rho_n_) * (dxlevel(ii)**2) &
       *divv(ixOOmin1:ixOOmax1,ixOOmin2:ixOOmax2,ixOOmin3:ixOOmax3)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_n_) = nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    !TH n  
    call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
       ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, temp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3), ii, tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3))
    nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = c_hyp(e_n_) * &
       csound_dim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) + c_shk(e_n_) * (dxlevel(ii)**2) *divv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = nu(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * Rn/(twofl_gamma-1d0)
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       e_n_) = nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    !visc n
    do jj=1,ndir
      call hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOOmin1,&
         ixOOmin2,ixOOmin3,ixOOmax1,ixOOmax2,ixOOmax3, vel(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,jj), ii, tmp(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3))
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = c_hyp(mom_n(jj)) * csound_dim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) * dxlevel(ii) *  tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) + c_shk(mom_n(jj)) * &
         (dxlevel(ii)**2) *divv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = nu(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) * rho(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
      wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(jj)) = nu(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    enddo
    !end for dim

  end function dump_hyperdiffusivity_coef_dim

  function dump_coll_terms(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, nwc) result(wnew)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nwc
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) 
    double precision   :: wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        1:nwc)
    double precision   :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
 
    call get_alpha_coll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
        tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)= tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
    call get_gamma_ion_rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
        tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)= tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
    wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)= tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 

  end function dump_coll_terms

  subroutine get_gamma_ion_rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, gamma_rec,&
      gamma_ion)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)       :: gamma_rec(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3),gamma_ion(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    ! calculations are done in S.I. units
    double precision, parameter :: A = 2.91e-14, & !m3/s
                                K = 0.39, XX = 0.232, Eion = 13.6 ! eV      
    double precision, parameter :: ECHARGE=1.6022d-19 !C
    double precision        :: rho(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call twofl_get_pthermal_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(Rc * rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3))

    !transform to SI units
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * unit_temperature * kB_SI/ECHARGE !* BK/ECHARGE means  K to eV 
    !number electrons rho_c = n_e * MH, in normalized units MH=1 and n = rho
    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3) * unit_numberdensity  
    if(.not. SI_unit) then
      !1/cm^3 = 1e6/m^3
      rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * 1d6
    endif
    gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) /sqrt(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) * 2.6e-19
    gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = ((rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * A) /(XX + Eion/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))) * ((Eion/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))**K) * exp(-Eion/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    ! see Voronov table: valid for temp min = 1eV(approx 11605 K), Temp max = 20KeV
    !to normalized
    gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * unit_time  
    gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * unit_time  

    if (associated(usr_mask_gamma_ion_rec)) then
      call usr_mask_gamma_ion_rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,gamma_ion,&
          gamma_rec)
    end if
  end subroutine get_gamma_ion_rec

  subroutine get_alpha_coll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, alpha)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)       :: alpha(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    if(twofl_alpha_coll_constant) then
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = twofl_alpha_coll
    else
      call get_alpha_coll_plasma(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
          alpha)
    endif
    if (associated(usr_mask_alpha)) then
      call usr_mask_alpha(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,alpha)
    end if
  end subroutine get_alpha_coll

  subroutine get_alpha_coll_plasma(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, alpha)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)       :: alpha(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision        :: pe(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),rho(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    double precision :: sigma_in = 1e-19 ! m2,m2,m2
    ! make calculation in SI physical units

    call twofl_get_pthermal_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pe)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = pe(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(Rc * rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    call twofl_get_pthermal_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,pe)
    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho)
    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = pe(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(Rn * rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3))
    alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = (2d0/(mp_SI**(3d0/2) * &
       sqrt(dpi))*sqrt(0.5*(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))*unit_temperature*kB_SI) * sigma_in)*unit_time * &
       unit_density
    if(.not. SI_unit) then
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * 1d3 !this comes from unit_density: g/cm3,g/cm3,g/cm3 = 1e-3 kg/m3,kg/m3,kg/m3
    endif

  end subroutine get_alpha_coll_plasma

  subroutine calc_mult_factor1(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, step_dt, JJ, res)
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: step_dt
    double precision, intent(in) :: JJ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out) :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = step_dt/(1d0 + &
       step_dt * JJ(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))

  end subroutine calc_mult_factor1

  subroutine calc_mult_factor2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, step_dt, JJ, res)
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: step_dt
    double precision, intent(in) :: JJ(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out) :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = (1d0 - exp(-step_dt &
       * JJ(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))/JJ(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine calc_mult_factor2

  subroutine advance_implicit_grid(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, wout, x,&
      dtfactor,qdt)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qdt
    double precision, intent(in) :: dtfactor
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    ::  x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)       :: wout(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:nw)

    integer :: idir
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp3(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp4(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp5(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: v_c(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ndir), v_n(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ndir)
    double precision :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        alpha(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, allocatable :: gamma_rec(:,:,:), gamma_ion(:,:,:)

    !TODO latest changes sets already wout to w in implicit update (see where psb=psa) 
    ! commment out setting mag and density when they are not modified here

    ! copy vars at the indices which are not updated here: mag. field
    wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mag(:)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(:))

    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
    !update density
    if(twofl_coll_inc_ionrec) then
       allocate(gamma_ion(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
           gamma_rec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)) 
       call get_gamma_ion_rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, gamma_rec,&
           gamma_ion)
       tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) =  gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) +  gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       call calc_mult_factor(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, dtfactor * qdt,&
           tmp2, tmp3) 
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = (-gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_n_) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_c_) - tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    else
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_n_)
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_c_)
    endif

    call get_alpha_coll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, alpha)
    
    !-J11 + J12    for momentum and kinetic energy
    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) =  alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * (rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) +  rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))     
    if(twofl_coll_inc_ionrec) then
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    endif
    call calc_mult_factor(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, dtfactor * qdt, tmp2,&
        tmp3) 

    ! momentum update
    do idir=1,ndir

      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)* (-rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idir)) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir)))
      if(twofl_coll_inc_ionrec) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_n(idir)) + gamma_rec(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_c(idir))
      endif

      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idir)) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir)) - tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    enddo

    ! energy update
    
    ! kinetic energy update  
    if(.not. phys_internal_e) then
      ! E_tot includes kinetic energy
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = twofl_kin_en_n(w,&
         ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
         ixOmin3,ixOmax1,ixOmax2,ixOmax3) 
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = twofl_kin_en_c(w,&
         ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
         ixOmin3,ixOmax1,ixOmax2,ixOmax3) 
      tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) - tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) !E_tot - E_kin
      tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      if(phys_total_energy) then
        tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - twofl_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
      endif

      !!implicit update
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*(-rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(twofl_coll_inc_ionrec) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif

      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) - tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)

     else 
      tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) 
      tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) 
      ! calculate velocities, using the already updated variables
      call twofl_get_v_n(wout,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v_n)
      call twofl_get_v_c(wout,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v_c)
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) 
      if(twofl_coll_inc_ionrec) then
        tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif 

      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 0.5d0 * sum((v_c(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir) - v_n(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir))**2, dim=ndim+1) * dtfactor * qdt 
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
     endif

    !update internal energy
    if(twofl_coll_inc_te) then
      if(has_equi_pe_n0) then
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)= block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_n0_,0)*inv_gamma_1  
      endif
      if(has_equi_pe_c0) then
        tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_c0_,0)*inv_gamma_1
      endif
      if (twofl_equi_thermal) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *(-1d0/Rn*(rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_c_)) + 1d0/Rc*(rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) +tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_n_))) 
       endif 
      if(has_equi_pe_n0) then
        tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) 
      endif
      if(has_equi_pe_c0) then
        tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) 
      endif
      if (.not. twofl_equi_thermal) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) *(-rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/Rn * tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/Rc * tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      endif
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) =  alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * (rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/Rc +  rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/Rn)     
      if(twofl_coll_inc_ionrec) then
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) =  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rc + gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rn 
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rn * tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rc * tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif
      call calc_mult_factor(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, dtfactor * qdt,&
          tmp2, tmp3) 
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    endif
    if(twofl_coll_inc_ionrec) then
       deallocate(gamma_ion, gamma_rec) 
    endif
  end subroutine advance_implicit_grid

  !> Implicit solve of psb=psa+dtfactor*dt*F_im(psb)
  subroutine twofl_implicit_coll_terms_update(dtfactor,qdt,qtC,psb,psa)
    use mod_global_parameters
    use mod_ghostcells_update

    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer :: iigrid, igrid
    !print*, "IMPL call ", it

    call getbc(global_time,0.d0,psa,1,nw)
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);
      block=>psa(igrid)
      call advance_implicit_grid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, psa(igrid)%w,&
          psb(igrid)%w, psa(igrid)%x, dtfactor,qdt)
    end do
    !$OMP END PARALLEL DO

  end subroutine twofl_implicit_coll_terms_update 

  !> inplace update of psa==>F_im(psa)
  subroutine twofl_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: iigrid, igrid, level

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);
      block=>psa(igrid)
       call coll_terms(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,&
          ixMlo3,ixMhi1,ixMhi2,ixMhi3,psa(igrid)%w,psa(igrid)%x)
    end do
    !$OMP END PARALLEL DO

  end subroutine twofl_evaluate_implicit

  subroutine coll_terms(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    integer :: idir
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp3(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp4(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp5(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    !double precision :: v_c(ixI^S,ndir), v_n(ixI^S,ndir)
    double precision, allocatable :: v_c(:,:,:,:), v_n(:,:,:,:)
    double precision, allocatable :: rho_c1(:,:,:), rho_n1(:,:,:)
    double precision :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        alpha(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, allocatable :: gamma_rec(:,:,:), gamma_ion(:,:,:)

    ! copy density before overwrite 
    if(twofl_equi_thermal) then
       allocate(rho_n1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
           rho_c1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)) 
       rho_n1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          rho_n_) 
       rho_c1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          rho_c_) 
    endif

    ! get total density before overwrite density
    call get_rhon_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    call get_rhoc_tot(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
    if(phys_internal_e) then
      ! get velocity before overwrite momentum
       allocate(v_n(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ndir),&
           v_c(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ndir)) 
      call twofl_get_v_n(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v_n)
      call twofl_get_v_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v_c)
    else
      ! get ke before overwrite density and momentum
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = twofl_kin_en_n(w,&
         ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
         ixOmin3,ixOmax1,ixOmax2,ixOmax3) 
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = twofl_kin_en_c(w,&
         ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,&
         ixOmin3,ixOmax1,ixOmax2,ixOmax3) 
    endif

    !update density
    if(twofl_coll_inc_ionrec) then
       allocate(gamma_ion(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
           gamma_rec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)) 
       call get_gamma_ion_rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, gamma_rec,&
           gamma_ion)
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = -gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          rho_n_) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          rho_c_) = -tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
    else
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_n_) = 0d0 
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_c_) = 0d0
  
    endif

    call get_alpha_coll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, alpha)

    ! momentum update
    do idir=1,ndir

      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)* (-rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idir)) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir)))
      if(twofl_coll_inc_ionrec) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_n(idir)) + gamma_rec(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_c(idir))
      endif

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idir)) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir)) = -tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
    enddo

    ! energy update
    
    ! kinetic energy update  
    if(.not. phys_internal_e) then
      ! E_tot includes kinetic energy
      tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) - tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) !E_tot - E_kin
      tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      if(phys_total_energy) then
        tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - twofl_mag_en(w,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
      endif
      ! tmp4 = eint_n, tmp5 = eint_c
      ! tmp1 = ke_n, tmp2 = ke_c
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*(-rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(twofl_coll_inc_ionrec) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = -tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 

     else 
      tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) 
      tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) 
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) 
      if(twofl_coll_inc_ionrec) then
        tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif 

      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 0.5d0 * sum((v_c(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir) - v_n(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir))**2, dim=ndim+1) 
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
     endif

    !update internal energy
    if(twofl_coll_inc_te) then

      if(has_equi_pe_n0) then
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)= block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_n0_,0)*inv_gamma_1  
      endif
      if(has_equi_pe_c0) then
        tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_c0_,0)*inv_gamma_1
      endif
      if (twofl_equi_thermal) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *(-1d0/Rn*(rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*rho_c1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)) + 1d0/Rc*(rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) +tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*rho_n1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))) 
       endif 
      if(has_equi_pe_n0) then
        tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) 
      endif
      if(has_equi_pe_c0) then
        tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) 
      endif
      if (.not. twofl_equi_thermal) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *(-rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rn * tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rc * tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
      endif

      if(twofl_coll_inc_ionrec) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rn * tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rc * tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    endif
    if(twofl_coll_inc_ionrec) then
       deallocate(gamma_ion, gamma_rec) 
    endif
    if(phys_internal_e) then
       deallocate(v_n, v_c) 
    endif
    if(twofl_equi_thermal) then
       deallocate(rho_n1, rho_c1) 
    endif
    !set contribution to mag field
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1:ndir)) = 0d0 

  end subroutine coll_terms

  subroutine twofl_explicit_coll_terms_update(qdt,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,&
     wCT,x)
    use mod_global_parameters
  
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    integer :: idir
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp3(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp4(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp5(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: v_c(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ndir), v_n(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ndir)
    double precision :: rhon(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        rhoc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        alpha(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, allocatable :: gamma_rec(:,:,:), gamma_ion(:,:,:)

    call get_rhon_tot(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhon)
    call get_rhoc_tot(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoc)
    !update density
    if(twofl_coll_inc_ionrec) then
       allocate(gamma_ion(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
           gamma_rec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)) 
       call get_gamma_ion_rec(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, wCT, x, gamma_rec,&
           gamma_ion)
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = qdt *(-gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_n_) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_c_) - tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
    endif

    call get_alpha_coll(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, wCT, x, alpha)

    ! momentum update
    do idir=1,ndir

      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)* (-rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_n(idir)) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mom_c(idir)))
      if(twofl_coll_inc_ionrec) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom_n(idir)) + gamma_rec(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3) * wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom_c(idir))
      endif
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) =tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * qdt

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_n(idir)) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom_c(idir)) - tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
    enddo

    ! energy update
    
    ! kinetic energy update  
    if(.not. phys_internal_e) then
      ! E_tot includes kinetic energy
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = twofl_kin_en_n(wCT,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) 
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = twofl_kin_en_c(wCT,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3) 
      tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,e_n_) - tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) !E_tot - E_kin
      tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,e_c_) - tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      if(phys_total_energy) then
        tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - twofl_mag_en(wCT,ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
      endif

      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*(-rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3))
      if(twofl_coll_inc_ionrec) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) =tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * qdt

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) - tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) 

    else 
      tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) 
      tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) 
      call twofl_get_v_n(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v_n)
      call twofl_get_v_c(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v_c)
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) 
      if(twofl_coll_inc_ionrec) then
        tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif 

      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 0.5d0 * sum((v_c(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir) - v_n(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1:ndir))**2, dim=ndim+1) * qdt 
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) + tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
     endif

    !update internal energy
    if(twofl_coll_inc_te) then
      if(has_equi_pe_n0) then
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)= block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_n0_,0)*inv_gamma_1  
      endif
      if(has_equi_pe_c0) then
        tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=block%equi_vars(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,equi_pe_c0_,0)*inv_gamma_1
      endif
      if (twofl_equi_thermal) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *(-1d0/Rn*(rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,rho_c_)) + 1d0/Rc*(rhon(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3) * tmp5(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3) +tmp3(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)*wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_n_))) 
       endif 
      if(has_equi_pe_n0) then
        tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) 
      endif
      if(has_equi_pe_c0) then
        tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) 
      endif
      if (.not. twofl_equi_thermal) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) *(-rhoc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rn * tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + rhon(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rc * tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
      endif

      if(twofl_coll_inc_ionrec) then
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) - gamma_ion(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rn * tmp4(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + gamma_rec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/Rc * tmp5(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif

      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) =tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * qdt

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_n_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         e_c_)-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    endif
    if(twofl_coll_inc_ionrec) then
       deallocate(gamma_ion, gamma_rec) 
    endif
  end subroutine twofl_explicit_coll_terms_update

  subroutine Rfactor_c(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Rfactor)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: Rfactor(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    Rfactor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=Rc

  end subroutine Rfactor_c

end module mod_twofl_phys
