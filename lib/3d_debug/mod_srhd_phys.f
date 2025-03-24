!> Special Relativistic Hydrodynamics (with EOS) physics module
module mod_srhd_phys
  use mod_physics
  use mod_constants
  use mod_comm_lib, only: mpistop
  implicit none
  private

  !> Whether particles module is added
  logical, public, protected              :: srhd_particles = .false.

  !> Number of tracer species
  integer, public,protected              :: srhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public,protected              :: rho_
  integer, public,protected              :: d_ 

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> Index of the energy density 
  integer, public,protected              :: e_
  !> Index of the gas pressure should equal e_
  integer, public,protected              :: p_

  !> Index of the Lorentz factor
  integer, public,protected     :: lfac_

  !> Index of the inertia
  integer, public,protected     :: xi_

  !> Whether synge eos is used
  logical, public                         :: srhd_eos = .false.

  !> The adiabatic index and derived values
  double precision, public                :: srhd_gamma = 5.d0/3.0d0
  double precision, public                :: gamma_1,inv_gamma_1,&
     gamma_to_gamma_1

  !> The smallest allowed energy
  double precision, public             :: small_e
  !> The smallest allowed inertia
  double precision, public             :: small_xi

  !> Allows overruling default corner filling (for debug mode, otherwise corner primitives fail)
  logical, public, protected              :: srhd_force_diagonal = .false.

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.0d0

  !> parameters for NR in con2prim
  integer, public                  :: maxitnr   = 100
  double precision, public         :: absaccnr  = 1.0d-8
  double precision, public         :: tolernr   = 1.0d-9
  double precision, public         :: dmaxvel   = 1.0d-7
  double precision, public         :: lfacmax
  double precision, public :: minp, minrho, smalltau, smallxi

  ! Public methods
  public :: srhd_phys_init
  public :: srhd_get_pthermal
  public :: srhd_get_auxiliary
  public :: srhd_get_auxiliary_prim
  public :: srhd_get_csound2
  public :: srhd_to_conserved
  public :: srhd_to_primitive
  public :: srhd_check_params
  public :: srhd_check_w 
  public :: srhd_get_Geff_eos
  public :: srhd_get_enthalpy_eos
contains

  !> Read this module's parameters from a file
  subroutine srhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srhd_list/ srhd_n_tracer, srhd_eos, srhd_gamma, srhd_particles,&
        srhd_force_diagonal, SI_unit, He_abundance

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, srhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine srhd_read_params

  !> Write this module's parameters to a snapshot
  subroutine srhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = srhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)

  end subroutine srhd_write_info

  !> Initialize the module
  subroutine srhd_phys_init()
    use mod_global_parameters
    use mod_particles, only: particles_init
    integer :: itr,idir

    call srhd_read_params(par_files)

    physics_type = "srhd"
    phys_energy  = .true.
    phys_total_energy  = .true.
    phys_gamma = srhd_gamma

    ! unused physics options
    phys_internal_e = .false.
    phys_partial_ionization=.false.
    phys_trac=.false.

    use_particles = srhd_particles

    ! note: number_species is 1 for srhd
    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=1

    ! Determine flux variables
    rho_ = var_set_rho()
    d_=rho_

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    e_ = var_set_energy()
    p_ = e_

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ! derive units from basic units
    call srhd_physical_units()

    if (srhd_force_diagonal) then
       ! ensure corners are filled, otherwise divide by zero when getting primitives
       !  --> only for debug purposes
       phys_req_diagonal = .true.
    endif

    allocate(tracer(srhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, srhd_n_tracer
       tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! Set index for auxiliary variables
    ! MUST be after the possible tracers (which have fluxes)
    xi_  = var_set_auxvar('xi','xi')
    lfac_= var_set_auxvar('lfac','lfac')

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1

    ! dummy for now, no extra source terms precoded
    phys_add_source          => srhd_add_source
    phys_get_dt              => srhd_get_dt
    ! copied in from HD/MHD, for certain limiters
    phys_get_a2max           => srhd_get_a2max

    ! actual srhd routines
    phys_check_params        => srhd_check_params
    phys_check_w             => srhd_check_w
    phys_get_cmax            => srhd_get_cmax
    phys_get_cbounds         => srhd_get_cbounds
    phys_get_flux            => srhd_get_flux
    phys_add_source_geom     => srhd_add_source_geom
    phys_to_conserved        => srhd_to_conserved
    phys_to_primitive        => srhd_to_primitive
    phys_get_pthermal        => srhd_get_pthermal
    phys_get_auxiliary       => srhd_get_auxiliary
    phys_get_auxiliary_prim  => srhd_get_auxiliary_prim
    phys_get_pthermal        => srhd_get_pthermal
    phys_get_v               => srhd_get_v
    phys_write_info          => srhd_write_info
    phys_handle_small_values => srhd_handle_small_values

    ! Initialize particles module
    if (srhd_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

  end subroutine srhd_phys_init

  subroutine srhd_check_params
    use mod_global_parameters

    if (srhd_gamma <= 0.0d0 .or. srhd_gamma == 1.0d0) call mpistop &
       ("Error: srhd_gamma <= 0 or srhd_gamma == 1")
    ! additional useful values
    gamma_1=srhd_gamma-1.0d0
    inv_gamma_1=1.0d0/gamma_1
    gamma_to_gamma_1=srhd_gamma/gamma_1
   
    ! the following sets small_e and small_xi from small_density/small_pressure
    ! according to the srhd_eos used
    call srhd_get_smallvalues_eos

    if(mype==0)then
       write(*,*&
          )'------------------------------------------------------------'
       write(*,*)'Using EOS set via srhd_eos=',srhd_eos
       write(*,*)'Maximal lorentz factor (via dmaxvel) is=',lfacmax
       write(*,*)'Use fixes set through check/fix small values:',&
           check_small_values,fix_small_values
       write(*,*)'Controlled with small pressure/density:', small_pressure,&
          small_density
       write(*,*)'Derived small values: xi and e ',small_xi,small_e
       write(*,*&
          )'------------------------------------------------------------'
    endif

  end subroutine srhd_check_params

  subroutine srhd_physical_units
    use mod_global_parameters
    double precision :: mp,kB

    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      unit_velocity=c_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      unit_velocity=const_c
    end if
    if(unit_numberdensity*unit_length<=0.0d0)then
       call mpistop&
          ("Abort: must set positive values for unit length and numberdensity")

    endif
    ! we assume user sets: unit_numberdensity, unit_length, He_abundance
    ! then together with light speed c, all units fixed
    unit_density=(1.0d0+4.0d0*He_abundance)*mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.0d0+&
       3.0d0*He_abundance)*unit_numberdensity*kB)
    unit_time=unit_length/unit_velocity
    unit_mass = unit_density*unit_length**3

  end subroutine srhd_physical_units

  !> Returns logical argument flag T where values are not ok
  subroutine srhd_check_w(primitive, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, flag)
    use mod_global_parameters
    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    logical, intent(inout)       :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    flag=.false.

    ! NOTE: we should not check or use nwaux variables here
    if(primitive) then
       where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,rho_) = .true.
       where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          p_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,p_) = .true.
    else
       where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          d_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,d_) = .true.
       where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          e_) < small_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          e_) = .true.
    endif

  end subroutine srhd_check_w

  !> Returns logical argument flag T where auxiliary values are not ok
  subroutine srhd_check_w_aux(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, flag)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    logical, intent(inout)       :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    flag=.false.

    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       xi_) < small_xi) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       xi_) = .true.
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_) < one) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_) = .true.

    if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)))then
       write(*,*)'auxiliary xi too low: abort'
       call mpistop('auxiliary  check failed')
    endif
    if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)))then
       write(*,*)'auxiliary lfac too low: abort'
       call mpistop('auxiliary  check failed')
    endif

  end subroutine srhd_check_w_aux

  !> Set auxiliary variables lfac and xi from a primitive state
  !> only used when handle_small_values average on primitives
  subroutine srhd_get_auxiliary_prim(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) :: rho,rhoh,pth

    ! assume four-velocity in momentum vector (i.e. lfac*v)
    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)    = &
       sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(:))**2,&
        dim=ndim+1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_) = dsqrt(1.0d0+rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
    ! compute rho*h (enthalpy h) from density-pressure
    call srhd_get_enthalpy_eos(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       rho,pth,rhoh)

    ! fill auxiliary variable xi= lfac^2 rhoh
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2.0d0*rhoh(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine srhd_get_auxiliary_prim

  !> Compute auxiliary variables lfac and xi from a conservative state
  !> using srhd_con2prim to calculate enthalpy and lorentz factor
  subroutine srhd_get_auxiliary(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    implicit none

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)

    integer                        :: ix1,ix2,ix3,ierror,idir
    integer                        :: flag_error(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision               :: ssqr

    if(srhd_eos)then
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        ierror=0
        ssqr=0.0d0
        do idir=1,ndir
          ssqr= ssqr+w(ix1,ix2,ix3,mom(idir))**2
        enddo
        if(w(ix1,ix2,ix3,d_)<small_density)then
           print *,'entering con2prim with', w(ix1,ix2,ix3,lfac_),w(ix1,ix2,ix3,xi_), &
                    w(ix1,ix2,ix3,d_),ssqr,w(ix1,ix2,ix3,e_)
           print *,'in position',ix1,ix2,ix3,x(ix1,ix2,ix3,1:ndim)
           print *,small_density,small_e,small_pressure,small_xi
           if(check_small_values) call mpistop(&
              'small density on entry con2prim')
           if(fix_small_values) w(ix1,ix2,ix3,d_)=small_density
        endif
        if(w(ix1,ix2,ix3,e_)<small_e)then
           print *,'entering con2prim with', w(ix1,ix2,ix3,lfac_),w(ix1,ix2,ix3,xi_), &
                    w(ix1,ix2,ix3,d_),ssqr,w(ix1,ix2,ix3,e_)
           print *,'in position',ix1,ix2,ix3,x(ix1,ix2,ix3,1:ndim)
           print *,small_density,small_e,small_pressure,small_xi
           if(check_small_values) call mpistop(&
              'small energy on entry con2prim')
           if(fix_small_values) w(ix1,ix2,ix3,e_)=small_e
        endif
        call con2prim_eos(w(ix1,ix2,ix3,lfac_),w(ix1,ix2,ix3,xi_), w(ix1,ix2,&
           ix3,d_),ssqr,w(ix1,ix2,ix3,e_),ierror)
        flag_error(ix1,ix2,ix3) = ierror
      enddo
      enddo
      enddo
    else
      do ix3=ixOmin3,ixOmax3
      do ix2=ixOmin2,ixOmax2
      do ix1=ixOmin1,ixOmax1
        ierror=0
        ssqr=0.0d0
        do idir=1,ndir
          ssqr= ssqr+w(ix1,ix2,ix3,mom(idir))**2
        enddo
        if(w(ix1,ix2,ix3,d_)<small_density)then
           print *,'entering con2prim with', w(ix1,ix2,ix3,lfac_),w(ix1,ix2,ix3,xi_), &
                    w(ix1,ix2,ix3,d_),ssqr,w(ix1,ix2,ix3,e_)
           print *,'in position',ix1,ix2,ix3,x(ix1,ix2,ix3,1:ndim)
           print *,small_density,small_e,small_pressure,small_xi
           if(check_small_values) call mpistop(&
              'small density on entry con2prim')
           if(fix_small_values) w(ix1,ix2,ix3,d_)=small_density
        endif
        if(w(ix1,ix2,ix3,e_)<small_e)then
           print *,'entering con2prim with', w(ix1,ix2,ix3,lfac_),w(ix1,ix2,ix3,xi_), &
                    w(ix1,ix2,ix3,d_),ssqr,w(ix1,ix2,ix3,e_)
           print *,'in position',ix1,ix2,ix3,x(ix1,ix2,ix3,1:ndim)
           print *,small_density,small_e,small_pressure,small_xi
           if(check_small_values) call mpistop(&
              'small energy on entry con2prim')
           if(fix_small_values) w(ix1,ix2,ix3,e_)=small_e
        endif
        call con2prim(w(ix1,ix2,ix3,lfac_),w(ix1,ix2,ix3,xi_), w(ix1,ix2,ix3,&
           d_),ssqr,w(ix1,ix2,ix3,e_),ierror)
        flag_error(ix1,ix2,ix3) = ierror
      enddo
      enddo
      enddo
    endif

    if(check_small_values)then
     if(any(flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        ixOmin3:ixOmax3)/=0))then
         print *,flag_error(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
         call mpistop('Problem when getting auxiliaries')
    !    call srhd_handle_small_values(.false.,w,x,ixI^L,ixO^L,'srhd_get_auxiliary')
     end if 
    end if 

  end subroutine srhd_get_auxiliary

  !> Transform primitive variables into conservative ones
  subroutine srhd_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    integer                         :: idir,itr
    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) :: rhoh,rho,pth

    ! assume four-velocity in momentum vector (i.e. lfac*v)
    ! use rhoh slot for temporary array
    rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)    = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3, mom(:))**2, dim=ndim+1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_) = dsqrt(1.0d0+rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))

    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
    pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
    ! compute rho*h (enthalpy h) from density-pressure
    call srhd_get_enthalpy_eos(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
       rho,pth,rhoh)

    ! compute rhoh*lfac (recycle rhoh)
    rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)= rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       lfac_)
    ! fill auxiliary variable xi= lfac^2 rhoh
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)*rhoh(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! set conservative density: d = lfac * rho
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)*rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! Convert four-velocity (lfac*v) to momentum (xi*v=[rho*h*lfac^2]*v)
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir)) = rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))
    end do 

    ! set tau = xi-p-d energy variable
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)-pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,d_)

    do itr=1,srhd_n_tracer
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          tracer(itr)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          d_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,tracer(itr))
    end do

  end subroutine srhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine srhd_to_primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    integer                         :: idir,itr
    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) :: rho,rhoh,E
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: pth
    character(len=30)                  :: subname_loc

    ! get auxiliary variables lfac and xi from conserved set
    call srhd_get_auxiliary(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

    ! from d to rho (d=rho*lfac)
    rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,lfac_)

    ! compute pressure
    ! deduce rho*h from xi/lfac^2
    rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,lfac_)**2.0d0
    call srhd_get_pressure_eos(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho,rhoh,pth,E)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)=rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    ! from xi*v to U=lfac*v (four-velocity) 
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         lfac_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)
    end do
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)=pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    do itr=1,srhd_n_tracer
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          tracer(itr)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          tracer(itr)) /(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          lfac_))
    end do

  end subroutine srhd_to_primitive

  !> Calculate v vector from conservatives
  subroutine srhd_get_v(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    integer :: idir

    ! get v from xi*v
    do idir=1,ndir
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(idir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)
    end do 

  end subroutine srhd_get_v

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> here computed from conservative set WITH ADDED rho*h
  !> local version: does not do con2prim
  subroutine srhd_get_csound2_rhoh(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rhoh,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw),rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    double precision                :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    rho=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       d_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
    call srhd_get_csound2_eos(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho,rhoh,csound2)

  end subroutine srhd_get_csound2_rhoh

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> here computed from conservative set and uses con2prim!!!
  !> public version!
  subroutine srhd_get_csound2(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out)   :: csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    double precision                :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! get auxiliary variables lfac and xi from conserved set
    call srhd_get_auxiliary(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    ! quantify rho, rho*h
    rho  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       d_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
    rhoh = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       xi_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2.0d0
    call srhd_get_csound2_eos(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho,rhoh,csound2)

  end subroutine srhd_get_csound2

  !> Calculate thermal pressure p within ixO^L
  !> must follow after update conservative with auxiliaries
  subroutine srhd_get_pthermal(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, pth)
    use mod_global_parameters
    use mod_small_values, only: trace_small_values

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out)   :: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    integer                      :: iw, ix1,ix2,ix3
    double precision             :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! quantify rho, rho*h, tau and get pthermal
    rho  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       d_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
    rhoh = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       xi_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2.0d0
    call srhd_get_pressure_eos(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho,rhoh,pth,E)

  end subroutine srhd_get_pthermal

  !> dummy addsource subroutine
  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine srhd_add_source(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
     wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt,dtfactor
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw),wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

  end subroutine srhd_add_source

  !> dummy get_dt subroutine
  subroutine srhd_get_dt(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, dtnew, dx1,dx2,dx3, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:3)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

  end subroutine srhd_get_dt

  !> Calculate cmax_idim within ixO^L
  !> used especially for setdt CFL limit
  subroutine srhd_get_cmax(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, cmax)
    use mod_global_parameters

    integer, intent(in)                       :: ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)        :: csound2,tmp1,tmp2,v2
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)        :: vidim, cmin

    logical       :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:nw)

    !!call srhd_check_w_aux(ixI^L, ixO^L, w, flag)

    ! auxiliaries are filled here
    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,lfac_)**2.0d0
    v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=1.0d0-&
       1.0d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)**2
    call srhd_get_csound2_rhoh(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp1,csound2)
    vidim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3, mom(idim))/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3, xi_)
    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=vidim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2.0d0
    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.0d0-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) -tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(1.0d0-csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=dsqrt(csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(one-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=vidim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(one-csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))/(one-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))/(one-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    ! Limit by speed of light
    cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = max(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), - 1.0d0)
    cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = min(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),   1.0d0)
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = max(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), - 1.0d0)
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = min(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),   1.0d0)
    ! now take extremal value only for dt limit
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = max(dabs(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)),dabs(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)))

  end subroutine srhd_get_cmax

  subroutine srhd_get_a2max(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
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
       jxOmax3,kxOmin1,kxOmin2,kxOmin3,kxOmax1,kxOmax2,kxOmax3,i

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
         1:nwflux)=dabs(-w(kxOmin1:kxOmax1,kxOmin2:kxOmax2,kxOmin3:kxOmax3,&
         1:nwflux)+16.d0*w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
         1:nwflux)-30.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)+16.d0*w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
         1:nwflux)-w(gxOmin1:gxOmax1,gxOmin2:gxOmax2,gxOmin3:gxOmax3,&
         1:nwflux))
      a2max(i)=maxval(a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i,&
         1:nwflux))/12.d0/dxlevel(i)**2
    end do

  end subroutine srhd_get_a2max

  !> local version for recycling code when computing cmax-cmin
  subroutine srhd_get_cmax_loc(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vidim,csound2,v2,cmax,&
     cmin)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)             :: vidim(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(in), dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) :: csound2
    double precision, intent(in)             :: v2(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(out)          :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(out)          :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3):: tmp1,tmp2

    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=vidim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**2.0d0
    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=1.0d0-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) -tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(1.0d0-csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=dsqrt(csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(one-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=vidim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*(one-csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))/(one-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    ! Limit by speed of light
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = max(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), - 1.0d0)
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = min(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),   1.0d0)
    cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))/(one-v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))
    ! Limit by speed of light
    cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = max(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), - 1.0d0)
    cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = min(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),   1.0d0)

  end subroutine srhd_get_cmax_loc

  !> Estimating bounds for the minimum and maximum signal velocities
  !> here we will not use Hspeed at all (one species only)
  subroutine srhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     idim,Hspeed,cmax,cmin)
    use mod_global_parameters
    use mod_variables

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    ! conservative left and right status
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    ! primitive left and right status
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    double precision :: wmean(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       nw)
    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) :: csound2,tmp1,tmp2,tmp3
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: vidim,cmaxL,cmaxR,cminL,cminR,v2

    logical       :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:nw)

    select case(boundspeed)
    case(1) ! we do left-right first and take maximals
      !!call srhd_check_w(.true.,ixI^L, ixO^L, wLp, flag)
      !!call srhd_check_w_aux(ixI^L, ixO^L, wLp, flag)
      tmp1=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
      tmp2=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         xi_)/wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         lfac_)**2.0d0
      tmp3=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
      call srhd_get_csound2_prim_eos(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3,tmp1,tmp2,tmp3,csound2)
      vidim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mom(idim))/wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, lfac_)
      v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1.0d0-1.0d0/wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,lfac_)**2
      call srhd_get_cmax_loc(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vidim,csound2,v2,&
         cmaxL,cminL)

      !!call srhd_check_w(.true.,ixI^L, ixO^L, wRp, flag)
      !!call srhd_check_w_aux(ixI^L, ixO^L, wRp, flag)
      tmp1=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
      tmp2=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         xi_)/wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         lfac_)**2.0d0
      tmp3=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
      call srhd_get_csound2_prim_eos(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3,tmp1,tmp2,tmp3,csound2)
      vidim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mom(idim))/wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, lfac_)
      v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1.0d0-1.0d0/wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,lfac_)**2
      call srhd_get_cmax_loc(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vidim,csound2,v2,&
         cmaxR,cminR)

      if(present(cmin))then
        ! for HLL
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(cmaxL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           cmaxR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=min(cminL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           cminR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      else
        ! for TVDLF
        cmaxL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=max(cmaxL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),dabs(cminL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)))
        cmaxR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=max(cmaxR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3),dabs(cminR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)))
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(cmaxL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           cmaxR(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      endif
    case(2) ! this is cmaxmean from conservatives
      ! here we do arithmetic mean of conservative vars
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux))
      ! get auxiliary variables
      call srhd_get_auxiliary(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wmean,x)
      ! here tmp1 is rhoh
      tmp1=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         xi_)/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         lfac_)**2.0d0
      call srhd_get_csound2_rhoh(wmean,x,ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp1,&
         csound2)
      vidim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mom(idim))/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, xi_)
      v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.0d0-1.0d0/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,lfac_)**2
      call srhd_get_cmax_loc(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vidim,csound2,v2,&
         cmaxL,cminL)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=cmaxL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=cminL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(cmaxL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           dabs(cminL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))
      endif
    case(3) ! this is cmaxmean from primitives
      ! here we do arithmetic mean of primitive vars
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)=0.5d0*(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1:nwflux))
      ! get auxiliary variables for wmean (primitive array)
      call srhd_get_auxiliary_prim(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wmean)
      ! here tmp1 is rhoh
      tmp1=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
      tmp2=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         xi_)/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         lfac_)**2.0d0
      tmp3=wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
      call srhd_get_csound2_prim_eos(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
         ixOmax3,tmp1,tmp2,tmp3,csound2)
      vidim(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, mom(idim))/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, lfac_)
      v2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1.0d0-1.0d0/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,lfac_)**2
      call srhd_get_cmax_loc(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vidim,csound2,v2,&
         cmaxL,cminL)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=cmaxL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=cminL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)=max(cmaxL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
           dabs(cminL(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))
      endif
    end select

  end subroutine srhd_get_cbounds

  !> Calculate fluxes within ixO^L.
  subroutine srhd_get_flux(wC,wP,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,f)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    ! primitive w
    double precision, intent(in) :: wP(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nwflux)

    double precision             :: pth(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision             :: v(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    integer                      :: iw,idir

    pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=wP(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,p_)
    do idir=1,ndir
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) = wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(idir))/wP(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
    end do 

    ! Get flux of density d, namely D*v 
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)=v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)*wC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

    ! Get flux of tracer
    do iw=1,srhd_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         tracer(iw))=v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k [+pth if i==k]
    do idir=1,ndir
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))= v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idim)*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))
    end do 
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim))=pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(idim))

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*pth
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)=v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,idim)*(wC(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_) + pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3))

  end subroutine srhd_get_flux

  !> Add geometrical source terms to w
  subroutine srhd_add_source_geom(qdt, dtfactor, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      wCT, w, x)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_surface
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, dtfactor, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        source(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)
    integer                         :: idir, h1xmin1,h1xmin2,h1xmin3,h1xmax1,&
       h1xmax2,h1xmax3, h2xmin1,h2xmin2,h2xmin3,h2xmax1,h2xmax2,h2xmax3
    integer :: mr_,mphi_,vr_,vphi_,vtheta_ ! Polar var. names
    double precision :: exp_factor(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), del_exp_factor(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), exp_factor_primitive(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    select case (coordinate)

    case(Cartesian_expansion)
      !the user provides the functions of exp_factor and del_exp_factor
      if(associated(usr_set_surface)) call usr_set_surface(ixImin1,ixImin2,&
         ixImin3,ixImax1,ixImax2,ixImax3,x,block%dx,exp_factor,del_exp_factor,&
         exp_factor_primitive)
      ! get auxiliary variables lfac and xi from conserved set
      call srhd_get_auxiliary(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,x)
      call srhd_get_pthermal(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, source)
      source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*del_exp_factor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/exp_factor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) + qdt*source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    case (cylindrical)
          mr_   = mom(r_)
          ! get auxiliary variables lfac and xi from conserved set
          call srhd_get_auxiliary(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
             ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,x)
          call srhd_get_pthermal(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,&
             ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
              source)
          if (phi_ > 0) then
             mphi_ = mom(phi_)
             vphi_ = mom(phi_)-1
             vr_   = mom(r_)-1
             call srhd_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)
             source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mphi_)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,vphi_)
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mr_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mr_) + qdt * source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, r_)
             source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) = -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, mphi_) * v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,vr_)
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mphi_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mphi_) + qdt * source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, r_)
          else
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mr_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mr_) + qdt * source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3, r_)
          end if

    case (spherical)
       mr_   = mom(r_)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmin3=ixOmin3-kr(1,3);h1xmax1=ixOmax1-kr(1,1)
       h1xmax2=ixOmax2-kr(1,2);h1xmax3=ixOmax3-kr(1,3)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmin3=ixOmin3-kr(2,3);h2xmax1=ixOmax1-kr(2,1)
       h2xmax2=ixOmax2-kr(2,2);h2xmax3=ixOmax3-kr(2,3);
       ! s[mr]=((stheta*vtheta+sphi*vphi)+2*p)/r
       ! get auxiliary variables lfac and xi from conserved set
       call srhd_get_auxiliary(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,x)
       call srhd_get_pthermal(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
          ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, pth)
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1) - block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,h1xmin3:h1xmax3,&
           1)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
       if (ndir > 1) then
         call srhd_get_v(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,v)
         do idir = 2, ndir
           source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) = source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3, mom(idir))*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,idir)
         end do
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mr_) + qdt * source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)

       
       vr_   = mom(r_)-1
       ! s[mtheta]=-(stheta*vr)/r+cot(theta)*(sphi*vphi+p)/r
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1) * (block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, 2) - block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          h2xmin3:h2xmax3, 2)) / block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
       if (ndir == 3) then
          source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, mom(3))*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,ndir)) / dtan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, 2))
       end if
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) - (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, mom(2)) * v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, vr_)) 
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(2)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(2)) + qdt * source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)

       if (ndir == 3) then
         vtheta_   = mom(2)-1
         ! s[mphi]=-(sphi*vr)/r-cot(theta)*(sphi*vtheta)/r
         source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) = -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, mom(3)) * v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, vr_)) - (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, mom(3)) * v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, vtheta_)) / dtan(x(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3, 2))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(3)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mom(3)) + qdt * source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, 1)
       end if
      
    end select

  end subroutine srhd_add_source_geom

  !> handles bootstrapping
  subroutine srhd_handle_small_values(primitive, w, x, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      subname)
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

    integer :: n,idir
    logical :: flag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw),&
       flagall(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call srhd_check_w(primitive, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, flag)

    if (any(flag)) then
      select case (small_values_method)
      case ("replace")
        ! any faulty cell is replaced by physical lower limit
        flagall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,rho_).or.flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,e_)) 

        where(flagall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
           ! D or rho: no difference primitive-conservative
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              rho_) = small_density
           !w(ixO^S,lfac_)= 1.0d0
           !w(ixO^S,xi_)  = small_xi
        endwhere
        !do idir = 1, ndir
        !   where(flagall(ixO^S)) w(ixO^S, mom(idir)) = 0.0d0
        !end do
        if(primitive) then
            where(flagall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, p_) = small_pressure
        else
            where(flagall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, e_) = small_e 
        endif

      case ("average")
        ! note: in small_values_average we use 
        ! small_values_fix_iw(1:nw) and small_values_daverage       
        ! when fails, may use small_pressure/small_density
        if(primitive)then
           ! averaging for all primitive fields (p, lfac*v, tau))
           call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
              ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
               flag)
           ! update the auxiliaries from primitives
           call srhd_get_auxiliary_prim(ixImin1,ixImin2,ixImin3,ixImax1,&
              ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
              w)
        else
           ! do averaging of density d
           call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
              ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
               flag, d_)
           ! do averaging of energy tau
           call small_values_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
              ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x,&
               flag, e_)
           ! and now hope for the best....
        endif
      case default
        if(.not.primitive) then
          ! note that we throw error here, which assumes w is primitive
          write(*,*)&
              "handle_small_values default: note reporting conservatives!"
        end if
        call small_values_error(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, flag,&
            subname)
      end select
    end if

  end subroutine srhd_handle_small_values

  !> calculate effective gamma
  subroutine srhd_get_Geff_eos(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,varconserve,Geff)
    use mod_global_parameters, only: nw
  !================== IMPORTANT ==================!
  !This subroutine is used with conserved variables in w when varconserve=T
  !This subroutine is used with primitive variables in w when varconserve=F
  !   both cases assume updated auxiliary variables xi_ en lfac_
  !===============================================!
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    logical, intent(in)                :: varconserve
    double precision, intent(out)      :: Geff(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) :: pth,rho,E_th,E

    if (srhd_eos) then
      if (varconserve) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)-w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)-w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)
        rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,d_)/w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)
        E_th = pth*inv_gamma_1
        E    = E_th+dsqrt(E_th**2+rho**2)
        Geff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = srhd_gamma-half*gamma_1 *          &
           (one-(rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))**2)
      else
        ! primitives available
        E_th = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           p_)*inv_gamma_1
        E    = E_th+dsqrt(E_th**2+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,rho_)**2)
        Geff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = srhd_gamma-half*gamma_1 *          &
           (one-(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           rho_)/E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))**2)
      end if
    else
      Geff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = srhd_gamma
    endif

  end subroutine srhd_get_Geff_eos

  !> Compute the small value limits
  subroutine srhd_get_smallvalues_eos
    use mod_global_parameters, only: small_pressure, small_density
    implicit none
    ! local small values
    double precision :: LsmallE,Lsmallp,Lsmallrho

    ! the maximal allowed Lorentz factor
    lfacmax=one/dsqrt(one-(one-dmaxvel)**2)
    minrho=small_density
    minp=small_pressure
    if(small_density*small_pressure<=0.0d0)then
       call mpistop&
("must set finite values small-density/pressure for small value treatments")
    endif
    if(srhd_eos)then
       Lsmallp=(one+10.d0*small_pressure)*small_pressure
       Lsmallrho=(one+10.d0*small_density)*small_density
       !!Lsmallp=small_pressure
       !!Lsmallrho=small_density
       LsmallE=Lsmallp*inv_gamma_1+dsqrt((Lsmallp*inv_gamma_1)**2+&
          Lsmallrho**2)
       small_xi=half*((srhd_gamma+one)*LsmallE-&
          gamma_1*Lsmallrho*(Lsmallrho/LsmallE))
       small_e=small_xi-Lsmallp-Lsmallrho
    else
       small_xi=small_density+gamma_to_gamma_1*small_pressure
       small_e =small_pressure*inv_gamma_1
    endif
    smallxi=small_xi
    smalltau=small_e

  end subroutine srhd_get_smallvalues_eos

  !> Compute the enthalpy rho*h from rho and pressure p
  subroutine srhd_get_enthalpy_eos(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3,rho,p,rhoh)
    use mod_global_parameters
    integer, intent(in)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, intent(in)       :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision, intent(out)      :: rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    double precision, dimension(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) :: E_th,E
    integer :: ix1,ix2,ix3

    if(srhd_eos) then
     E_th = p*inv_gamma_1
     E    = E_th+dsqrt(E_th**2+rho**2)
     ! writing rho/E on purpose, for numerics 
     rhoh = half*((srhd_gamma+one)*E - gamma_1*rho*(rho/E))
    else
     rhoh = rho+gamma_to_gamma_1*p
    end if

    if (check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(rhoh(ix1,ix2,ix3)<small_xi) then
           write(*,*) "local pressure and density",p(ix1,ix2,ix3),rho(ix1,ix2,&
              ix3)
           write(*,*) "Error: small value of enthalpy rho*h=",rhoh(ix1,ix2,&
              ix3)," encountered when call srhd_get_enthalpy_eos"
           call mpistop&
              ('enthalpy below small_xi: stop (may need to turn on fixes)')
         end if
      enddo
      enddo
      enddo
    end if

    if (fix_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(rhoh(ix1,ix2,ix3)<small_xi) then
            rhoh(ix1,ix2,ix3)=small_xi
         endif
      enddo
      enddo
      enddo
    endif

  end subroutine srhd_get_enthalpy_eos

  !> Calculate thermal pressure p from density rho and enthalpy rho*h 
  !> will provide p (and E if srhd_eos)
  subroutine srhd_get_pressure_eos(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho,rhoh,p,E)
    use mod_global_parameters
    integer, intent(in)            :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)   :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision, intent(out)  :: p(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out)  :: E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer :: ix1,ix2,ix3

    if(srhd_eos) then
     E = (rhoh+dsqrt(rhoh**2+(srhd_gamma**2-one)*rho**2)) /(srhd_gamma+one)
     p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = half*gamma_1* &
        (E-rho*(rho/E))
    else 
     p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
        (rhoh-rho)/gamma_to_gamma_1
    end if 

    if (check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(p(ix1,ix2,ix3)<small_pressure) then
           write(*,*) "local enthalpy rho*h and density rho",rhoh(ix1,ix2,ix3),&
              rho(ix1,ix2,ix3)
           if(srhd_eos) write(*,*) 'E, rho^2/E, difference', E(ix1,ix2,ix3),&
              rho(ix1,ix2,ix3)**2/E(ix1,ix2,ix3),E(ix1,ix2,ix3)-rho(ix1,ix2,&
              ix3)**2/E(ix1,ix2,ix3)
           write(*,*) "Error: small value of gas pressure",p(ix1,ix2,ix3),&
              " encountered when call srhd_get_pressure_eos"
           call mpistop&
('pressure below small_pressure: stop (may need to turn on fixes)')
         end if
      enddo
      enddo
      enddo
    end if

    if (fix_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(p(ix1,ix2,ix3)<small_pressure) then
            p(ix1,ix2,ix3)=small_pressure
            if(srhd_eos)E(ix1,ix2,ix3)=max(small_e,E(ix1,ix2,ix3))
         endif
      enddo
      enddo
      enddo
    endif

  end subroutine srhd_get_pressure_eos

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> available rho - rho*h 
  subroutine srhd_get_csound2_eos(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho,rhoh,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision, intent(out)   :: csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    double precision                :: p(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision                :: E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer :: ix1,ix2,ix3

    call srhd_get_pressure_eos(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,rho,rhoh,p,E)
    if(srhd_eos) then
       csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=(p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*((srhd_gamma+one)+gamma_1*(rho(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3)/E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))**2))/(2.0d0*rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3))
    else
       csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=srhd_gamma*p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)/rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)
    end if

    if (check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(csound2(ix1,ix2,ix3)>=1.0d0.or.csound2(ix1,ix2,ix3)<=0.0d0) then
           write(*,*) "sound speed error with p - rho - rhoh",p(ix1,ix2,ix3),&
              rhoh(ix1,ix2,ix3),rho(ix1,ix2,ix3)
           if(srhd_eos) write(*,*) 'and E', E(ix1,ix2,ix3)
           write(*,*) "Error: value of csound2",csound2(ix1,ix2,ix3),&
              " encountered when call srhd_get_csound2_eos"
           call mpistop('sound speed stop (may need to turn on fixes)')
         end if
      enddo
      enddo
      enddo
    end if

    if (fix_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(csound2(ix1,ix2,ix3)>=1.0d0) then
            csound2(ix1,ix2,ix3)=1.0d0-1.0d0/lfacmax**2
         endif
         if(csound2(ix1,ix2,ix3)<=0.0d0) then
            csound2(ix1,ix2,ix3)=srhd_gamma*small_pressure/small_xi
         endif
      enddo
      enddo
      enddo
    endif

  end subroutine srhd_get_csound2_eos

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> available rho - rho*h - p
  subroutine srhd_get_csound2_prim_eos(ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
     ixOmax3,rho,rhoh,p,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3
    double precision, intent(in)    :: rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3),rhoh(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
       p(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision, intent(out)   :: csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    double precision                :: E(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer :: ix1,ix2,ix3

    if(srhd_eos) then
       E = (rhoh+dsqrt(rhoh**2+(srhd_gamma**2-one)*rho**2))/(srhd_gamma+one)
       csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=(p*((srhd_gamma+one)+&
          gamma_1*(rho/E)**2))/(2.0d0*rhoh)
    else
       csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)=srhd_gamma*p/rhoh
    end if

    if (check_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(csound2(ix1,ix2,ix3)>=1.0d0.or.csound2(ix1,ix2,ix3)<=0.0d0) then
           write(*,*) "sound speed error with p - rho - rhoh",p(ix1,ix2,ix3),&
              rhoh(ix1,ix2,ix3),rho(ix1,ix2,ix3)
           if(srhd_eos) write(*,*) 'and E', E(ix1,ix2,ix3)
           write(*,*) "Error: value of csound2",csound2(ix1,ix2,ix3),&
              " encountered when call srhd_get_csound2_prim_eos"
           call mpistop('sound speed stop (may need to turn on fixes)')
         end if
      enddo
      enddo
      enddo
    end if

    if (fix_small_values) then
      do ix3= ixOmin3,ixOmax3
      do ix2= ixOmin2,ixOmax2
      do ix1= ixOmin1,ixOmax1
         if(csound2(ix1,ix2,ix3)>=1.0d0) then
            csound2(ix1,ix2,ix3)=1.0d0-1.0d0/lfacmax**2
         endif
         if(csound2(ix1,ix2,ix3)<=0.0d0) then
            csound2(ix1,ix2,ix3)=srhd_gamma*small_pressure/small_xi
         endif
      enddo
      enddo
      enddo
    endif

  end subroutine srhd_get_csound2_prim_eos

  !> con2prim: (D,S**2,tau) --> compute auxiliaries lfac and xi
  subroutine con2prim_eos(lfac,xi,myd,myssqr,mytau,ierror)
    use mod_con2prim_vars

    double precision, intent(in)    :: myd, myssqr, mytau
    double precision, intent(inout) :: lfac, xi
    integer, intent(inout)          :: ierror

    ! .. local ..
    double precision:: f,df,lfacl
    !------------------------------------------------------------------

    ! Save the input-state in mod_con2prim_vars 
    d = myd; ssqr = myssqr; tau = mytau; 

    ierror=0

    ! Check if guess is close enough: gives f,df,lfacl
    if(xi>smallxi)then
    call funcd_eos(xi,f,df,lfacl,d,ssqr,tau,ierror)
    if (ierror == 0 .and. dabs(f/df)<absaccnr) then
       xi   = xi - f/df
       lfac = lfacl
       return
    else
       ierror = 0
    end if
    else
      write(*,*)'entering con2prim_eos with xi=',xi
    end if

    ! ierror=1 : must have D>=minrho, tau>=smalltau
    !if(d<minrho .or. tau<smalltau)then
    !   ierror=1
    !   return
    !endif

    call con2primHydro_eos(lfac,xi,d,ssqr,tau,ierror)

  end subroutine con2prim_eos

  subroutine funcd_eos(xi,f,df,mylfac,d,ssqr,tau,ierror)
    double precision, intent(in)  :: xi,d,ssqr,tau
    double precision, intent(out) :: f,df,mylfac
    integer, intent(inout)        :: ierror

    ! .. local ..
    double precision  :: dlfac
    double precision  :: vsqr,p,dpdxi
    !-----------------------------------------------------------------

    vsqr = ssqr/xi**2

    if (vsqr<one) then
       mylfac = one/dsqrt(one-vsqr)
       dlfac = -mylfac**3*ssqr/(xi**3)
       !===== Pressure, calculate using EOS =====!
       call FuncPressure_eos(xi,mylfac,d,dlfac,p,dpdxi)
       !=========================================!
       f  = xi-tau-d-p
       df = one-dpdxi
    else
       ! print *,'Erroneous input to funcd since vsqr=',vsqr,' >=1'
       ! print *,'input values d, ssqr, tau:',d,ssqr,tau
       ierror =6
       return
    end if

  end subroutine funcd_eos

  !> SRHD iteration solves for p via NR, and then gives xi as output
  subroutine con2primHydro_eos(lfac,xi,d,sqrs,tau,ierror)
    double precision, intent(out) :: xi,lfac
    double precision, intent(in)  :: d,sqrs,tau
    integer,intent(inout)         :: ierror

    ! .. local ..
    integer          :: ni,niiter,nit,n2it,ni3
    double precision :: pcurrent,pnew
    double precision :: er,er1,ff,df,dp,v2
    double precision :: pmin,lfac2inv,pLabs,pRabs,pprev
    double precision :: s2overcubeG2rh
    double precision :: xicurrent,h,dhdp
    double precision :: oldff1,oldff2,Nff
    double precision :: pleft,pright
    !---------------------------------------------------------------------

    ierror=0
    ! ierror=0 : ok
    !  we already checked D>=minrho, tau>=smalltau (ierror=1)
    !
    ! ierror<>0
    !
    ! ierror=2 : maxitnr reached without convergence
    ! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
    ! ierror=4 : final v^2=1 hence problem as lfac=1/0
    ! ierror=5 : nonmonotonic function f (as df=0)

    ! left and right brackets for p-range
    pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
    pLabs=max(minp,pmin)
    pRabs=1.0d99
    ! start value from input
    pcurrent=pLabs

    er1=one
    pprev=pcurrent

    ! Fudge Parameters
    oldff1=1.0d7  ! High number
    oldff2=1.0d9  ! High number bigger then oldff1
    n2it = 0
    nit  = 0

    LoopNR:  do ni=1,maxitnr
       nit = nit + 1
       !=== Relax NR iteration accuracy=======!
       if(nit>maxitnr/4)then
          ! mix pressure value for convergence
          pcurrent=half*(pcurrent+pprev)
          ! relax accuracy requirement
          er1=10.0d0*er1
          nit = nit - maxitnr/10
       endif
       !=======================================!

       niiter=ni
       xicurrent=tau+d+pcurrent

       if(xicurrent<smallxi) then
          !        print *,'stop: too small xi iterate:',xicurrent
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          ierror=3
          return
       endif

       v2=sqrs/xicurrent**2
       lfac2inv=one - v2
       if(lfac2inv>zero) then
          lfac=one/dsqrt(lfac2inv)
       else
          !        print *,'stop: negative or zero factor 1-v2:',lfac2inv
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          !        print *,'values for v2,xi:',v2,xicurrent
          ierror=4
          return
       endif

       s2overcubeG2rh=sqrs/(xicurrent**3)
       !== calculation done using the EOS ==!
       call FuncEnthalpy_eos(pcurrent,lfac2inv,d,sqrs,xicurrent,s2overcubeG2rh,&
          h,dhdp)
       !=======================================!
       ff=-xicurrent*lfac2inv + h
       df=- two*sqrs/xicurrent**2  + dhdp - lfac2inv

       if (ff*df==zero) then
          if (ff==zero) then
             exit ! zero found
          else
             !     print *,'stop: df becomes zero, non-monotonic f(p)!'
             ierror=5
             return
          endif
       else
          pnew=pcurrent-ff/df
          if (ff*df>zero) then
             ! pressure iterate has decreased
             ! restrict to left
             pnew=max(pnew,pLabs)
          else  ! ff*df<0
             ! pressure iterate has increased
             ! restrict to right
             pnew=min(pnew,pRabs)
          endif
       endif

       !===============================================!
       dp=pcurrent-pnew
       er=two*dabs(dp)/(pnew+pcurrent)
       if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
       !===============================================!

       ! For very small values of pressure, NR algorithm is not efficient to
       ! find root, use Euler algorithm to find precise value of pressure
       if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.ff &
          * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

          n2it=n2it+1
          if(n2it<=3) pcurrent=half*(pnew+pcurrent)
          if(n2it>3)then
             pright =pcurrent
             pleft=pprev
             pcurrent=half*(pleft+pright)
             Dicho:  do ni3=1,maxitnr
                !===================!
                xicurrent=tau+d+pcurrent
                v2=sqrs/xicurrent**2
                lfac2inv=one - v2

                if(lfac2inv>zero)then
                   lfac=one/dsqrt(lfac2inv)
                else
                   ierror=4
                   return
                endif
                !===================!

                !== calculation done using the EOS ==!
                call Bisection_Enthalpy_eos(pnew,lfac2inv,d,xicurrent,h)
                Nff=-xicurrent*lfac2inv + h
                !=======================================!
                !==== Iterate ====!
                if(ff * Nff < zero)then
                   pleft=pcurrent
                else
                   pright=pcurrent
                endif

                pcurrent=half*(pleft+pright)
                !==================!

                !=== The iteration converged ===!
                if(2.0d0*dabs(pleft-pright)/(pleft+pright)< absaccnr .or. &
                   dabs(ff)<absaccnr)then
                   pnew=pcurrent
                   exit LoopNR
                endif
                !==============================!

                !==============================!

                !=== conserve the last value of Nff ===!
                ff=Nff
                !======================================!
             enddo    Dicho
          endif

       else
          !====== There is no problems, continue the NR iteration ======!
          pprev=pcurrent
          pcurrent=pnew
          !=============================================================!
       endif


       !=== keep the values of the 2 last ff ===!
       oldff2=oldff1
       oldff1=ff
       !========================================!
    enddo LoopNR

    if(niiter==maxitnr)then
       ierror=2
       return
    endif

    if(pcurrent<minp) then
       ierror=3
       return
    endif

    !------------------------------!
    xi=tau+d+pcurrent
    v2=sqrs/xicurrent**2
    lfac2inv=one - v2
    if(lfac2inv>zero) then
       lfac=one/dsqrt(lfac2inv)
    else
       ierror=4
       return
    endif

  end subroutine con2primHydro_eos

  !> pointwise evaluations used in con2prim
  !> compute pointwise value for pressure p and dpdxi
  subroutine FuncPressure_eos(xicurrent,lfac,d,dlfacdxi,p,dpdxi)

    double precision, intent(in)         :: xicurrent,lfac,d,dlfacdxi
    double precision, intent(out)        :: p,dpdxi
    ! .. local ..
    double precision                     :: rho,h,E,dhdxi,rhotoE
    double precision                     :: dpdchi,dEdxi

    ! rhoh here called h
    h=xicurrent/(lfac**2)
    rho=d/lfac
    E = (h+dsqrt(h**2+(srhd_gamma**2-one)*rho**2)) /(srhd_gamma+one)
    ! output pressure
    rhotoE = rho/E
    p = half*gamma_1*(E-rho*rhotoE)

    dhdxi = one/(lfac**2)-2.0d0*xicurrent/(lfac**2)*dlfacdxi/lfac

    dEdxi=(dhdxi+(h*dhdxi-(srhd_gamma**2-one)*rho**2*dlfacdxi/lfac)/dsqrt(h**2+&
       (srhd_gamma**2-one)*rho**2))/(srhd_gamma+one)

    ! output pressure derivative to xi
    dpdxi=half*gamma_1*(2.0d0*rho*rhotoE*dlfacdxi/lfac+(one+rhotoE**2)*dEdxi)

  end subroutine FuncPressure_eos

  !> pointwise evaluations used in con2prim
  !> returns enthalpy rho*h (h) and derivative d(rho*h)/dp (dhdp)
  subroutine FuncEnthalpy_eos(pcurrent,lfac2inv,d,sqrs,xicurrent,dv2d2p,h,&
     dhdp)

    double precision, intent(in) :: pcurrent,lfac2inv,d,sqrs,xicurrent,dv2d2p
    double precision, intent(out):: h,dhdp

    ! local
    double precision:: rho,E_th,E,dE_thdp,dEdp

    rho=d*dsqrt(lfac2inv)
    E_th = pcurrent*inv_gamma_1
    E = (E_th + dsqrt(E_th**2+rho**2))
    !== Enthalpy ==!
    h = half*((srhd_gamma+one)*E-gamma_1*rho*(rho/E))
    !=== Derivative of thermal energy ===!
    dE_thdp = one*inv_gamma_1
    !=== Derivative of internal energy ===!
    dEdp = dE_thdp * (one+E_th/dsqrt(E_th**2+rho**2))+  &
       d**2*dv2d2p/dsqrt(E_th**2+rho**2)
    !====== Derivative of Enthalpy ======!
    dhdp = half*((srhd_gamma+one)*dEdp + &
       gamma_1*(rho*(rho/E))*(-2.0d0*dv2d2p/lfac2inv+dEdp/E))
  end subroutine FuncEnthalpy_eos

  !> pointwise evaluations used in con2prim
  !> returns enthalpy rho*h (h)
  subroutine Bisection_Enthalpy_eos(pcurrent,lfac2inv,d,xicurrent,h)

    double precision, intent(in) :: pcurrent,lfac2inv,d,xicurrent
    double precision, intent(out):: h

    ! local
    double precision:: rho,E_th,E

    rho=d*dsqrt(lfac2inv)
    E_th = pcurrent*inv_gamma_1
    E = (E_th + dsqrt(E_th**2+rho**2))
    !== Enthalpy ==!
    h = half*((srhd_gamma+one)*E-gamma_1*rho*(rho/E))

    return
  end subroutine Bisection_Enthalpy_eos

  !> con2prim: (D,S**2,tau) --> compute auxiliaries lfac and xi
  subroutine con2prim(lfac,xi,myd,myssqr,mytau,ierror)
    use mod_con2prim_vars

    double precision, intent(in)    :: myd, myssqr, mytau
    double precision, intent(inout) :: lfac, xi
    integer, intent(inout)          :: ierror

    ! .. local ..
    double precision:: f,df,lfacl
    !------------------------------------------------------------------

    ! Save the input-state in mod_con2prim_vars 
    d = myd; ssqr = myssqr; tau = mytau; 

    ierror=0

    ! Check if guess is close enough: gives f,df,lfacl
    if(xi>smallxi)then
    call funcd(xi,f,df,lfacl,d,ssqr,tau,ierror)
    if (ierror == 0 .and. dabs(f/df)<absaccnr) then
       xi   = xi - f/df
       lfac = lfacl
       return
    else
       ierror = 0
    end if
    else
       write(*,*) 'entering con2prim with xi=',xi
    end if

    ! ierror=1 : must have D>=minrho, tau>=smalltau
    !if(d<minrho .or. tau<smalltau)then
    !   ierror=1
    !   return
    !endif

    call con2primHydro(lfac,xi,d,ssqr,tau,ierror)

  end subroutine con2prim

  subroutine funcd(xi,f,df,mylfac,d,ssqr,tau,ierror)
    double precision, intent(in)  :: xi,d,ssqr,tau
    double precision, intent(out) :: f,df,mylfac
    integer, intent(inout)        :: ierror

    ! .. local ..
    double precision  :: dlfac
    double precision  :: vsqr,p,dpdxi
    !-----------------------------------------------------------------

    vsqr = ssqr/xi**2

    if (vsqr<one) then
       mylfac = one/dsqrt(one-vsqr)
       dlfac = -mylfac**3*ssqr/(xi**3)
       !===== Pressure, calculate using EOS =====!
       call FuncPressure(xi,mylfac,d,dlfac,p,dpdxi)
       !=========================================!
       f  = xi-tau-d-p
       df = one-dpdxi
    else
       ! print *,'Erroneous input to funcd since vsqr=',vsqr,' >=1'
       ! print *,'input values d, ssqr, tau:',d,ssqr,tau
       ierror =6
       return
    end if

  end subroutine funcd

  !> SRHD iteration solves for p via NR, and then gives xi as output
  subroutine con2primHydro(lfac,xi,d,sqrs,tau,ierror)
    double precision, intent(out) :: xi,lfac
    double precision, intent(in)  :: d,sqrs,tau
    integer,intent(inout)         :: ierror

    ! .. local ..
    integer          :: ni,niiter,nit,n2it,ni3
    double precision :: pcurrent,pnew
    double precision :: er,er1,ff,df,dp,v2
    double precision :: pmin,lfac2inv,pLabs,pRabs,pprev
    double precision :: s2overcubeG2rh
    double precision :: xicurrent,h,dhdp
    double precision :: oldff1,oldff2,Nff
    double precision :: pleft,pright
    !---------------------------------------------------------------------

    ierror=0
    ! ierror=0 : ok
    !  we already checked D>=minrho, tau>=smalltau (ierror=1)
    !
    ! ierror<>0
    !
    ! ierror=2 : maxitnr reached without convergence
    ! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
    ! ierror=4 : final v^2=1 hence problem as lfac=1/0
    ! ierror=5 : nonmonotonic function f (as df=0)

    ! left and right brackets for p-range
    pmin=dsqrt(sqrs)/(one-dmaxvel)-tau-d
    pLabs=max(minp,pmin)
    pRabs=1.0d99
    ! start value from input
    pcurrent=pLabs

    er1=one
    pprev=pcurrent

    ! Fudge Parameters
    oldff1=1.0d7  ! High number
    oldff2=1.0d9  ! High number bigger then oldff1
    n2it = 0
    nit  = 0

    LoopNR:  do ni=1,maxitnr
       nit = nit + 1
       !=== Relax NR iteration accuracy=======!
       if(nit>maxitnr/4)then
          ! mix pressure value for convergence
          pcurrent=half*(pcurrent+pprev)
          ! relax accuracy requirement
          er1=10.0d0*er1
          nit = nit - maxitnr/10
       endif
       !=======================================!

       niiter=ni
       xicurrent=tau+d+pcurrent

       if(xicurrent<smallxi) then
          !        print *,'stop: too small xi iterate:',xicurrent
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          ierror=3
          return
       endif

       v2=sqrs/xicurrent**2
       lfac2inv=one - v2
       if(lfac2inv>zero) then
          lfac=one/dsqrt(lfac2inv)
       else
          !        print *,'stop: negative or zero factor 1-v2:',lfac2inv
          !        print *,'for pressure iterate p',pcurrent
          !        print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
          !        print *,'iteration number:',ni
          !        print *,'values for d,s,tau,s2:',d,sqrs,tau,sqrs
          !        print *,'values for v2,xi:',v2,xicurrent
          ierror=4
          return
       endif

       s2overcubeG2rh=sqrs/(xicurrent**3)
       !== calculation done using the EOS ==!
       call FuncEnthalpy(pcurrent,lfac2inv,d,sqrs,xicurrent,s2overcubeG2rh,h,&
          dhdp)
       !=======================================!
       ff=-xicurrent*lfac2inv + h
       df=- two*sqrs/xicurrent**2  + dhdp - lfac2inv

       if (ff*df==zero) then
          if (ff==zero) then
             exit ! zero found
          else
             !     print *,'stop: df becomes zero, non-monotonic f(p)!'
             ierror=5
             return
          endif
       else
          pnew=pcurrent-ff/df
          if (ff*df>zero) then
             ! pressure iterate has decreased
             ! restrict to left
             pnew=max(pnew,pLabs)
          else  ! ff*df<0
             ! pressure iterate has increased
             ! restrict to right
             pnew=min(pnew,pRabs)
          endif
       endif

       !===============================================!
       dp=pcurrent-pnew
       er=two*dabs(dp)/(pnew+pcurrent)
       if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
       !===============================================!

       ! For very small values of pressure, NR algorithm is not efficient to
       ! find root, use Euler algorithm to find precise value of pressure
       if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.ff &
          * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

          n2it=n2it+1
          if(n2it<=3) pcurrent=half*(pnew+pcurrent)
          if(n2it>3)then
             pright =pcurrent
             pleft=pprev
             pcurrent=half*(pleft+pright)
             Dicho:  do ni3=1,maxitnr
                !===================!
                xicurrent=tau+d+pcurrent
                v2=sqrs/xicurrent**2
                lfac2inv=one - v2

                if(lfac2inv>zero)then
                   lfac=one/dsqrt(lfac2inv)
                else
                   ierror=4
                   return
                endif
                !===================!

                !== calculation done using the EOS ==!
                call Bisection_Enthalpy(pnew,lfac2inv,d,xicurrent,h)
                Nff=-xicurrent*lfac2inv + h
                !=======================================!
                !==== Iterate ====!
                if(ff * Nff < zero)then
                   pleft=pcurrent
                else
                   pright=pcurrent
                endif

                pcurrent=half*(pleft+pright)
                !==================!

                !=== The iteration converged ===!
                if(2.0d0*dabs(pleft-pright)/(pleft+pright)< absaccnr .or. &
                   dabs(ff)<absaccnr)then
                   pnew=pcurrent
                   exit LoopNR
                endif
                !==============================!

                !==============================!

                !=== conserve the last value of Nff ===!
                ff=Nff
                !======================================!
             enddo    Dicho
          endif

       else
          !====== There is no problems, continue the NR iteration ======!
          pprev=pcurrent
          pcurrent=pnew
          !=============================================================!
       endif


       !=== keep the values of the 2 last ff ===!
       oldff2=oldff1
       oldff1=ff
       !========================================!
    enddo LoopNR

    if(niiter==maxitnr)then
       ierror=2
       return
    endif

    if(pcurrent<minp) then
       ierror=3
       return
    endif

    !------------------------------!
    xi=tau+d+pcurrent
    v2=sqrs/xicurrent**2
    lfac2inv=one - v2
    if(lfac2inv>zero) then
       lfac=one/dsqrt(lfac2inv)
    else
       ierror=4
       return
    endif

  end subroutine con2primHydro

  !> pointwise evaluations used in con2prim
  !> compute pointwise value for pressure p and dpdxi
  subroutine FuncPressure(xicurrent,lfac,d,dlfacdxi,p,dpdxi)

    double precision, intent(in)         :: xicurrent,lfac,d,dlfacdxi
    double precision, intent(out)        :: p,dpdxi
    ! .. local ..
    double precision                     :: rho,h,E,dhdxi,rhotoE
    double precision                     :: dpdchi,dEdxi

    ! rhoh here called h
    h=xicurrent/(lfac**2)
    rho=d/lfac
    ! output pressure
    p = (h - rho)/gamma_to_gamma_1
    dpdchi = one/gamma_to_gamma_1
    dpdxi = dpdchi * one/lfac**2
    ! zero case dlfacdxi implies zero velocity (ssqr=0)
    if (dlfacdxi /= 0.0d0) dpdxi = dpdxi  + dpdchi * &
       ((d*lfac-2.0d0*xicurrent)/lfac**3) * dlfacdxi

  end subroutine FuncPressure

  !> pointwise evaluations used in con2prim
  !> returns enthalpy rho*h (h) and derivative d(rho*h)/dp (dhdp)
  subroutine FuncEnthalpy(pcurrent,lfac2inv,d,sqrs,xicurrent,dv2d2p,h,dhdp)

    double precision, intent(in) :: pcurrent,lfac2inv,d,sqrs,xicurrent,dv2d2p
    double precision, intent(out):: h,dhdp

    ! local
    double precision:: rho,E_th,E,dE_thdp,dEdp

    rho=d*dsqrt(lfac2inv)
    h = rho + gamma_to_gamma_1 * pcurrent
    dhdp = gamma_to_gamma_1 + d/dsqrt(lfac2inv)*sqrs/xicurrent**3
  end subroutine FuncEnthalpy

  !> pointwise evaluations used in con2prim
  !> returns enthalpy rho*h (h)
  subroutine Bisection_Enthalpy(pcurrent,lfac2inv,d,xicurrent,h)

    double precision, intent(in) :: pcurrent,lfac2inv,d,xicurrent
    double precision, intent(out):: h

    ! local
    double precision:: rho,E_th,E

    rho=d*dsqrt(lfac2inv)
    h = rho + gamma_to_gamma_1 * pcurrent

    return
  end subroutine Bisection_Enthalpy

end module mod_srhd_phys
