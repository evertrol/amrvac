!> Module for including dust species, which interact with the gas through a drag
!> force
module mod_dust
  use mod_global_parameters, only: std_len
  use mod_physics
  use mod_comm_lib, only: mpistop

  implicit none
  private

  !> The number of dust species
  integer, public, protected      :: dust_n_species = 0

  integer, protected              :: gas_rho_ = -1
  integer, allocatable, protected :: gas_mom(:)
  integer, protected              :: gas_e_   = -1

  !> Indices of the dust densities
  integer, allocatable, public, protected :: dust_rho(:)

  !> Indices of the dust momentum densities
  integer, allocatable, public, protected :: dust_mom(:, :)

  !> Size of each dust species, dimensionless expression
  double precision, allocatable, public :: dust_size(:)

  !> Internal density of each dust species, dimensionless expression
  double precision, allocatable, public :: dust_density(:)

  !> Reduction of stopping time timestep limit
  double precision :: dust_dtpar = 0.5d0

  !> Factor used in squared thermal velocity
  double precision :: gas_vtherm_factor = 3.0d0

  !> Dust temperature in K (if dust_temperature_type is constant)
  double precision :: dust_temperature = -1.0d0

  !> Dust drag coefficient for linear drag (for testing dust_method=linear)
  double precision :: dust_K_lineardrag = -1.0d0

  !> If dust_temperature_type is stellar, it will be calculated according to Tielens (2005),
  !> eqn. 5.44 using an input stellar luminosity in solar luminosities
  double precision :: dust_stellar_luminosity = -1.0d0

  !> Set small dust densities to zero to avoid numerical problems
  logical, public, protected :: dust_small_to_zero = .false.

  !> Minimum dust density as used when dust_small_to_zero=T
  double precision, public, protected :: dust_min_rho = -1.0d0

  !> Adding dust in sourcesplit manner or not
  logical :: dust_source_split = .false.

  !> This can be turned off for testing purposes, if F then gas uncouples from dust
  logical :: dust_backreaction = .true.

  !> What type of dust drag force to use. Can be 'Kwok', 'sticking', 'linear', 'usr' or 'none'.
  character(len=std_len), public, protected :: dust_method = 'Kwok'

  !> Can be 'graphite' or 'silicate', affects the dust temperature
  character(len=std_len) :: dust_species = 'graphite'

  !> Determines the dust temperature, can be 'constant', 'ism', or 'stellar'
  character(len=std_len) :: dust_temperature_type = 'constant'

  !> whether second order terms (relevant only when dust_n_species >=2) are included
  !> there are the terms  n2, ni2, d2 in Eqs 6,7,8 in amrvac 3.0 paper
  logical :: dust_implicit_second_order = .true.  

  !> whether fh is added for gas energy:  is only added in the impliict implementation, the explicit one was left as before
  logical :: dust_backreaction_fh = .false.  


  ! Public methods
  public :: dust_init
  public :: dust_get_dt
  public :: dust_get_flux
  public :: dust_get_cmax
  public :: dust_get_flux_prim
  public :: dust_get_cmax_prim
  public :: dust_add_source
  public :: dust_to_conserved
  public :: dust_to_primitive
  public :: dust_check_params
  public :: dust_check_w
  public :: set_dusttozero
  public :: dust_implicit_update
  public :: dust_evaluate_implicit


contains

  subroutine dust_init(g_rho, g_mom, g_energy)
    use mod_global_parameters

    integer, intent(in) :: g_rho
    integer, intent(in) :: g_mom(ndir)
    integer, intent(in) :: g_energy ! Negative value if not present
    integer             :: n, idir
    character(len=2)    :: dim

    call dust_read_params(par_files)

    if(dust_source_split) any_source_split=.true.
    allocate(gas_mom(ndir))
    gas_rho_ = g_rho
    gas_mom  = g_mom
    gas_e_   = g_energy

    allocate(dust_size(dust_n_species))
    allocate(dust_density(dust_n_species))
    dust_size(:) = -1.0d0
    dust_density(:) = -1.0d0

    allocate(dust_rho(dust_n_species))
    allocate(dust_mom(ndir, dust_n_species))

    ! Set index of dust densities
    do n = 1, dust_n_species
      dust_rho(n) = var_set_fluxvar("rhod", "rhod", n)
    end do

    ! Dust momentum
    do idir = 1, ndir
      write(dim, "(I0,A)") idir, "d"
      do n = 1, dust_n_species
        dust_mom(idir, n) = var_set_fluxvar("m"//dim, "v"//dim, n)
      end do
    end do

  end subroutine dust_init

  !> Read dust_list module parameters from a file
  subroutine dust_read_params(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /dust_list/ dust_n_species, dust_min_rho, dust_method,&
        dust_K_lineardrag, dust_small_to_zero, dust_source_split,&
        dust_temperature, dust_temperature_type, dust_backreaction, dust_dtpar,&
        gas_vtherm_factor, dust_stellar_luminosity,dust_implicit_second_order,&
        dust_backreaction_fh 

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, dust_list, end=111)
111   close(unitpar)
    end do

  end subroutine dust_read_params

  subroutine dust_check_params()
    use mod_usr_methods, only: usr_get_3d_dragforce,usr_dust_get_dt
    use mod_global_parameters, only : mype, SI_unit, use_imex_scheme

    if (dust_method == 'sticking') then
       if (SI_unit) call mpistop("Dust error: sticking assumes cgs units")
       if (dust_temperature_type == "constant") then
          if (dust_temperature < 0.0d0) then
             call mpistop("Dust error: dust_temperature (in K) < 0 or not set")

          end if
       else if (dust_temperature_type == "stellar") then
          if (dust_stellar_luminosity < 0.0d0) then
             call mpistop&
("Dust error: dust_stellar_luminosity (in solar) < 0 or not set")
          end if
       end if
    end if

    if (dust_method == 'linear') then
        if(dust_K_lineardrag<0.0d0) then
          call mpistop&
("With dust_method=='linear', you must set a positive dust_K_lineardrag")
       end if
    end if

    if (any(dust_size < 0.0d0)) call mpistop(&
       "Dust error: any(dust_size < 0) or not set")
    if (any(dust_density < 0.0d0)) call mpistop(&
       "Dust error: any(dust_density < 0) or not set")

    if (dust_method == 'usr') then
       if (.not. associated(usr_get_3d_dragforce) .or. .not. &
          associated(usr_dust_get_dt)) call &
          mpistop&
          ("Dust error:usr_get_3d_dragforce and usr_dust_get_dt not defined")
    end if

    if(.not. use_imex_scheme .and. ((dust_dtpar .ge. &
       1d0).or.(dust_dtpar.le.0))) then
      if(mype .eq. 0) print*,&
          "EXPLICIT source for dust requires 0<dt_dustpar < 1, set to 0.8"
      dust_dtpar = 0.8
    endif

  end subroutine dust_check_params

  subroutine dust_check_w(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,flag)
    use mod_global_parameters
    
    integer, intent(in)         :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in):: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(inout)      :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer                     :: n

    do n = 1, dust_n_species
       flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          dust_rho(n))=(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          dust_rho(n))<0.0d0)
    enddo

  end subroutine dust_check_w

  subroutine dust_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    integer                         :: n, idir

    if(fix_small_values .and. dust_small_to_zero) call set_dusttozero(ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3, w, x)

    do n = 1, dust_n_species
      ! Convert velocity to momentum
      do idir = 1, ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
            n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            dust_rho(n)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            dust_mom(idir, n))
      end do
    end do

  end subroutine dust_to_conserved

  subroutine dust_to_primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    integer                         :: n, idir

    do n = 1, dust_n_species
      ! Convert momentum to velocity
      do idir = 1, ndir
        where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            dust_rho(n)) > 0.0d0)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
              n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_mom(idir, n)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, dust_rho(n))
        elsewhere
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
              n)) = 0.0d0
        end where
      end do
    end do

    if(fix_small_values .and. dust_small_to_zero) call set_dusttozero(ixImin1,&
       ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,&
       ixOmax1,ixOmax2,ixOmax3, w, x)

  end subroutine dust_to_primitive

  subroutine dust_get_flux(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nwflux)
    integer                         :: n, idir

    do n = 1, dust_n_species
       where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           dust_rho(n)) > 0.0d0)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_mom(idim, n))
       elsewhere
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n)) = 0.0d0
       end where

       do idir = 1, ndir
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
            n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            dust_mom(idir, n)) * get_vdust(w, ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
            idim, n)
       end do
    end do

  end subroutine dust_get_flux

  subroutine dust_get_flux_prim(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nwflux)
    integer                         :: n, idir

    do n = 1, dust_n_species
       where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           dust_rho(n)) > 0.0d0)
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_mom(idim, n))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, dust_rho(n))
       elsewhere
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n)) = 0.0d0
       end where

       do idir = 1, ndir
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
            n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            dust_mom(idir, n)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, dust_rho(n)) * get_vdust_prim(w, ixImin1,ixImin2,&
           ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
           ixOmax2,ixOmax3, idim, n)
       end do
    end do

  end subroutine dust_get_flux_prim

  function get_vdust(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, n) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, n
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        dust_rho(n)) > 0.0d0)
      vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          dust_mom(idim, n)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, dust_rho(n))
    elsewhere
      vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    end where

  end function get_vdust

  function get_vdust_prim(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, n) result(vdust)
    use mod_global_parameters, only: nw
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, n
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision              :: vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        dust_rho(n)) > 0.0d0)
      vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          dust_mom(idim, n))
    elsewhere
      vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
    end where

  end function get_vdust_prim

  ! Force dust density to zero if dust_rho <= dust_min_rho
  subroutine set_dusttozero(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    logical                         :: flag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                         :: n, idir

    do n = 1, dust_n_species
      flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=(w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n)) <= dust_min_rho)
      where (flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            dust_rho(n)) = 0.0d0
      end where
      do idir = 1, ndir
        where (flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
              n)) = 0.0d0
        end where
      end do
    end do

  end subroutine set_dusttozero

  ! Calculate drag force based on Epstein's law
  ! From Kwok 1975, page 584 (between eqn 8 and 9)
  subroutine get_3d_dragforce(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, fdrag, ptherm,&
      vgas)
    use mod_global_parameters
    use mod_usr_methods, only: usr_get_3d_dragforce
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(out)   :: fdrag(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndir, 1:dust_n_species)
    double precision, intent(in)    :: ptherm(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), vgas(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        1:ndir)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: vt2, deltav, fd, vdust
    double precision                   :: alpha_T(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:dust_n_species)
    integer                            :: n, idir

    vt2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       gas_vtherm_factor*ptherm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        gas_rho_)

    select case( TRIM(dust_method) )
    case ('Kwok') ! assume sticking coefficient equals 0.25

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n)) > 0.0d0)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, dust_mom(idir, n)) / w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)

            ! 0.75 from sticking coefficient
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = 0.75d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, dust_rho(n))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, gas_rho_)*deltav(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3) / (dust_density(n) * &
               dust_size(n))

            ! 0.75 from spherical grainvolume
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = -fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)*0.75d0*dsqrt(vt2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3) + deltav(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2)
          elsewhere
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
          end where
          fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
              n) = fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
      end do

    case ('sticking') !Calculate sticking coefficient based on the gas and dust temperatures

      call get_sticking(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, alpha_T, ptherm)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n))>0.0d0)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,dust_mom(idir, n)) / w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = (one-alpha_T(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,n)) * w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n))*w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                gas_rho_) * deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) / (dust_density(n)*dust_size(n))
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = -fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)*0.75d0*dsqrt(vt2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3) + deltav(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2)
          else where
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
          end where
          fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
             n) = fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
      end do

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n))>0.0d0)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,dust_mom(idir, n))/w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)

            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = -dust_K_lineardrag*deltav(&
               ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          else where
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
          end where
          fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
             n) = fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
      end do

    case('usr')
      call usr_get_3d_dragforce(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x, fdrag,&
          ptherm, vgas, dust_n_species)
    case('none')
      fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, :, :) = 0.0d0
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select

  end subroutine get_3d_dragforce

  !> Get sticking coefficient alpha_T (always between 0 and 1)
  !>
  !> Uses Temperatures in K
  !> Equation from Decin et al. 2006
  subroutine get_sticking(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, alpha_T,&
      ptherm)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(out) :: alpha_T(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:dust_n_species)
    double precision, intent(in)  :: ptherm(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision              :: Tgas(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    integer                       :: n

    ! get the dust species T in K
    call get_tdust(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, alpha_T)

    ! convert dimensionless gas T to K
    Tgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = (ptherm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        gas_rho_))*unit_temperature

    do n = 1, dust_n_species
      alpha_T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         n) =  0.35d0 * dexp(-dsqrt((Tgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + alpha_T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,n))/5.0d2))+0.1d0
    end do

  end subroutine get_sticking

  !> Returns dust temperature (in K), either as constant or based on equ. 5.41,
  !> 5.42 and 5.44 from Tielens (2005)
  !>
  !> Note that this calculation assumes cgs!!!! 
  !>
  !> It takes as input the stellar luminosity in solar units in 'stellar' case
  !> or a fixed dust input temperature in Kelvin when 'constant' or does case 'ism'
  subroutine get_tdust(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, Td)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(out) :: Td(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:dust_n_species)
    double precision              :: G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                       :: n

    select case( trim(dust_temperature_type) )
    case( 'constant' )
      Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          :) = dust_temperature
    case( 'ism' )
      select case( trim(dust_species) )
      case( 'graphite' )
        do n = 1, dust_n_species
          Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              n) = 15.8d0*((0.0001d0/(dust_size(n)*unit_length))**0.06d0)
        end do
      case( 'silicate' )
        do n = 1, dust_n_species
          Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              n) = 13.6d0*((0.0001d0/(dust_size(n)*unit_length))**0.06d0)
        end do
      case default
        call mpistop( "=== Dust species undetermined===" )
      end select
    case( 'stellar' )
      select case(coordinate)
      case(spherical)
        G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = max(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, 1)*unit_length, smalldouble)
      !!!case(cylindrical) convert R,Z to spherical radial coordinate r here
      !!! but only ok for 2D (R,Z) or 2.5D (R,Z) case
      !!! G0(ixO^S) = max(dsqrt(sum(x(ixO^S,:)**2,dim=ndim+1))*unit_length, smalldouble)
      case default
        call mpistop('stellar case not available in this coordinate system')
      end select

      G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 2.1d4*(dust_stellar_luminosity/1.0d8)*((&
         3.0857d17/G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))**2)

      select case( trim(dust_species) )
      case( 'graphite' )
        do n = 1, dust_n_species
          Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              n) = 61.0d0*((0.0001d0/(dust_size(n)*unit_length))**0.06d0) &
             *(G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)**(one/5.8d0))
        end do
      case( 'silicate' )
        do n = 1, dust_n_species
          Td(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              n) = 50.0d0*((0.0001d0/(dust_size(n)*unit_length))**0.06d0) &
             *(G0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)**(one/6.0d0))
        end do
      case default
        call mpistop( "=== Dust species undetermined===" )
      end select
    case default
      call mpistop( "=== Dust temperature undetermined===" )
    end select

  end subroutine get_tdust

  !> w[iw]= w[iw]+qdt*S[wCT,  x] where S is the source based on wCT within ixO
  subroutine dust_add_source(qdt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, wCT, w, x,&
      qsourcesplit, active)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    double precision :: ptherm(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), vgas(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        1:ndir)
    double precision :: fdrag(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        1:ndir, 1:dust_n_species)
    integer          :: n, idir

    select case( TRIM(dust_method) )
    case( 'none' )
      !do nothing here
    case default !all regular dust methods here
      if (qsourcesplit .eqv. dust_source_split) then
        active = .true.

        call phys_get_pthermal(wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ptherm)
        do idir=1,ndir
          vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             idir)=wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             gas_mom(idir))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,gas_rho_)
        end do

        call get_3d_dragforce(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
            ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, wCT, x, fdrag,&
            ptherm, vgas)
        fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, 1:ndir,&
            1:dust_n_species) = fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, 1:ndir, 1:dust_n_species) * qdt

        do idir = 1, ndir

          do n = 1, dust_n_species
            if (dust_backreaction) then
               w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   gas_mom(idir))  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  ixOmin3:ixOmax3, gas_mom(idir)) + fdrag(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir, n)
               if (gas_e_ > 0) then
                 w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     gas_e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    ixOmin3:ixOmax3, gas_e_) + vgas(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                     idir)  * fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    ixOmin3:ixOmax3, idir, n)
               end if
            end if


            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
                n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                dust_mom(idir, n)) - fdrag(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, idir, n)
          end do
        end do

      endif
    end select

  end subroutine dust_add_source

  !> inplace update of psa==>F_im(psa)
  subroutine dust_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: iigrid, igrid, level

    !dust_method = 'none' not used

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);
      block=>psa(igrid)
       call dust_terms(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,&
          ixMlo3,ixMhi1,ixMhi2,ixMhi3,psa(igrid)%w,psa(igrid)%x)
    end do
    !$OMP END PARALLEL DO

  end subroutine dust_evaluate_implicit




  subroutine dust_terms(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        vgas(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, 1:ndir)
    double precision :: alpha(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        1:ndir, 1:dust_n_species)
    integer          :: n, idir

    do idir=1,ndir
      vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         gas_mom(idir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         gas_rho_)
    end do
    call get_alpha_dust(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, vgas,x, alpha)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, gas_e_)=0d0
    do idir = 1, ndir

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, gas_mom(idir))=0d0
      do n = 1, dust_n_species
        ! contribution for gas momentum
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, idir,n) * (  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,dust_rho(n)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, gas_mom(idir)) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,gas_rho_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, dust_mom(idir, n)))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
            n)) = -tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        if (dust_backreaction) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              gas_mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, gas_mom(idir)) + tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          if (gas_e_ > 0) then
            if(dust_backreaction_fh) then
              where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 dust_rho(n)) > 0d0)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                    gas_e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3, gas_e_) + alpha(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
                   n) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   gas_rho_) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3, dust_mom(idir,n))**2/w(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   dust_rho(n))) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,dust_rho(n)) * (w(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                    gas_mom(idir))**2/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,gas_rho_)))
              elsewhere
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                    gas_e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3, gas_e_) + alpha(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
                   n) * ( - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                   dust_rho(n)) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3, gas_mom(idir))**2/w(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3,gas_rho_)))
              endwhere
            else  
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  gas_e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  gas_e_) + vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3, idir)  * tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
            end if
          end if
        end if
      end do
    end do
  end subroutine dust_terms

    !> Implicit solve of psb=psa+dtfactor*dt*F_im(psb)
  subroutine dust_implicit_update(dtfactor,qdt,qtC,psb,psa)
    use mod_global_parameters
    !use mod_ghostcells_update

    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer :: iigrid, igrid

    !call getbc(global_time,0.d0,psa,1,nw)
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);
      block=>psa(igrid)
      call dust_advance_implicit_grid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
         ixGhi3, ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3, psa(igrid)%w,&
          psb(igrid)%w, psa(igrid)%x, dtfactor,qdt)
    end do
    !$OMP END PARALLEL DO

   end subroutine dust_implicit_update 

  subroutine dust_advance_implicit_grid(ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, wout,&
      x, dtfactor,qdt)
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

    integer                            :: n, m, idir
    double precision :: alpha(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        1:ndir, 1:dust_n_species)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: tmp3(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision    :: vgas(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndir)


    do idir = 1, ndir
      vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         gas_mom(idir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         gas_rho_)
    end do
    call get_alpha_dust(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, vgas, x, alpha)
    !TODO this is still neeed?
    wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1:nw) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:nw)

    do idir = 1, ndir
      ! d1 from Eq 6
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
      do n = 1, dust_n_species
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) +  alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, idir,n) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,gas_rho_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,dust_rho(n))) 

      enddo
      ! store D in tmp
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = 1d0 + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) * qdt 
      if(dust_implicit_second_order) then
        ! d2 from Eq 6
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
        do n = 1, dust_n_species
          do m = n+1, dust_n_species
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) +  alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, idir,n) * alpha(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,m) *(w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,gas_rho_) + w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,dust_rho(n))+w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,dust_rho(m)))
          enddo
        enddo
        ! multiplied at the end by rho_gas 
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,gas_rho_)*tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) * (qdt**2)
      endif



      do n = 1, dust_n_species
        ! ni1 from eq 7
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, idir,n) * (  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,dust_rho(n)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, gas_mom(idir)) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,gas_rho_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3, dust_mom(idir, n))) * qdt

        if(dust_implicit_second_order) then 
          ! ni2 from eq 7 
          tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
          do m = n+1, dust_n_species
              tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) +  alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3, idir,n) * alpha(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
                 m) * ( w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 dust_rho(n)) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3, dust_mom(idir, n)) +  w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  gas_mom(idir))) - (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,gas_rho_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,dust_rho(m))) * w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir, n)) )  
          enddo
          ! tmp3 multiplied at the end by rho_gas 
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,gas_rho_)* (qdt**2) 
        endif
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
           n)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            dust_mom(idir,n)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      enddo

      if (dust_backreaction) then
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
        !n1 from eq 8 
        do n = 1, dust_n_species
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, idir,n) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,gas_rho_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, dust_mom(idir,n)) - w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,dust_rho(n)) * w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, gas_mom(idir)))

        enddo
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = qdt *  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) 
        if(dust_implicit_second_order) then 
          !n2 from eq 8 
          tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
          do n = 1, dust_n_species
            do m = n+1, dust_n_species
               tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  ixOmin3:ixOmax3) = tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  ixOmin3:ixOmax3) + alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  ixOmin3:ixOmax3, idir,n) * alpha(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
                  m) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  gas_rho_) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  ixOmin3:ixOmax3, dust_mom(idir, n)) + w(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
                   m))) - (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                  dust_rho(n)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                  ixOmin3:ixOmax3,dust_rho(m)))* w(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2,ixOmin3:ixOmax3, gas_mom(idir)))  
            enddo
          enddo
          ! tmp3 multiplied at the end by rho_gas 
          tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + (qdt**2)*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)* w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,gas_rho_)
        endif
        ! store in tmp2 contribution to momentum
        ! so that it is used when dust_backreaction_fh = .false.
        tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) / tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
        wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            gas_mom(idir)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            gas_mom(idir)) + tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)

        ! kinetic energy update
         if (gas_e_ > 0) then
          if(dust_backreaction_fh) then 
            ! add work done by coll terms + FrictionalHeating
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
            do n = 1, dust_n_species
              ! 2*dust kinetic energy: dust rho can be 0
              where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 dust_rho(n)) > 0d0)
                tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)= w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3, dust_mom(idir,n))**2/w(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2,ixOmin3:ixOmax3,dust_rho(n))
              elsewhere
                tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
              endwhere
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) + alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3, idir,n) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,gas_rho_) * tmp3(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3) - w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 dust_rho(n)) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3, gas_mom(idir))**2/w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,gas_rho_)))
  
            enddo
            tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = qdt *  tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) 
            if(dust_implicit_second_order) then
              tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
              do n = 1, dust_n_species
                do m = n+1, dust_n_species
                    tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3) = tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3) + alpha(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
                       n) * alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3, idir,m) * (w(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                       gas_rho_) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3, dust_mom(idir,&
                        n))**2/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3,dust_rho(n)) + w(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_mom(idir,&
                       m))**2/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3,dust_rho(m))) - (w(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                       dust_rho(n)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3,dust_rho(m)))* w(ixOmin1:ixOmax1,&
                       ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                        gas_mom(idir))**2/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       ixOmin3:ixOmax3,gas_rho_))  
                enddo
              enddo
              ! tmp3 multiplied at the end by rho_gas 
              tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) + (qdt**2)*tmp3(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3)* w(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,ixOmin3:ixOmax3,gas_rho_)
            endif
            wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                gas_e_) = wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                gas_e_) + 0.5d0 * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) / tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
          else
            ! dust_backreaction_fh = .false.
            ! add only work done by coll term by multiplyting the contribution in mom eq. by velocity
            wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                gas_e_) = wout(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                gas_e_) + vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                idir) * tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          endif
         end if
      end if
    end do !1..ndir
        

  end subroutine dust_advance_implicit_grid

  ! copied from  get_3d_dragforce subroutine
  subroutine get_alpha_dust(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, vgas,x, alpha)
    use mod_global_parameters
    use mod_usr_methods, only: usr_get_3d_dragforce
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision,intent(in)    :: vgas(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndir)
    double precision, intent(out)   :: alpha(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndir, 1:dust_n_species)

    double precision    :: ptherm(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: vt2, deltav, fd, vdust
    double precision                   :: alpha_T(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:dust_n_species)
    integer                            :: n, idir

    call phys_get_pthermal(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ptherm)

    vt2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       gas_vtherm_factor*ptherm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
        gas_rho_)

    select case( TRIM(dust_method) )
    case ('Kwok') ! assume sticking coefficient equals 0.25

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n)) > 0.0d0)

            ! 0.75 from sticking coefficient
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = 0.75d0 / (dust_density(n) * &
               dust_size(n))

            ! 0.75 from spherical grainvolume
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, dust_mom(idir, n)) / w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)*0.75d0*dsqrt(vt2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3) + deltav(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2)
          elsewhere
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
          end where
          alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
              n) = fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
      end do

    case ('sticking') !Calculate sticking coefficient based on the gas and dust temperatures

      call get_sticking(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, alpha_T, ptherm)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n))>0.0d0)
            ! sticking
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = (one-alpha_T(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               n)) / (dust_density(n)*dust_size(n))
            ! 0.75 from spherical grainvolume
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,dust_mom(idir, n)) / w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)*0.75d0*dsqrt(vt2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3) + deltav(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3)**2)
          else where
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
          end where
          alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
             n) = fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
      end do

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n))>0.0d0)
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)     = dust_K_lineardrag/(w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,gas_rho_)*w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n)))
          else where
            fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
          end where
          alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, idir,&
             n) = fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        end do
      end do

    case('none')
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, :, :) = 0.0d0
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select

  end subroutine get_alpha_dust


  !> Get dt related to dust and gas stopping time (Laibe 2011)
  subroutine dust_get_dt(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, dtnew, dx1,dx2,dx3, x)
    use mod_global_parameters
    use mod_usr_methods, only: usr_dust_get_dt

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(inout) :: dtnew

    double precision                :: ptherm(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), vgas(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        1:ndir)
    double precision, dimension(1:dust_n_species):: dtdust
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)           :: vt2, deltav, tstop, vdust
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:dust_n_species) :: alpha_T
    integer                                    :: n, idir

    if(dust_dtpar .le. 0) return

    call phys_get_pthermal(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, ptherm)
    do idir = 1, ndir
      vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         gas_mom(idir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         gas_rho_)
    end do

    select case( TRIM(dust_method) )

    case( 'Kwok' ) ! assume sticking coefficient equals 0.25
      dtdust(:) = bigdouble

      vt2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = gas_vtherm_factor*ptherm(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, gas_rho_)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n))>0.0d0)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,dust_mom(idir, n))/w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
            tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)  = 4.0d0*(dust_density(n)*dust_size(n))/ &
               (3.0d0*(0.75d0)*dsqrt(vt2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) + deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)**2)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, dust_rho(n)) + w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, gas_rho_)))
          else where
            tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = bigdouble
          end where

          dtdust(n) = min(minval(tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)), dtdust(n))
        end do
      end do

      dtnew = min(minval(dust_dtpar*dtdust(:)), dtnew)

    case( 'sticking' ) !Calculate sticking coefficient based on the gas temperature
      dtdust(:) = bigdouble

      vt2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = gas_vtherm_factor*ptherm(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3, gas_rho_)

      ! Sticking coefficient
      call get_sticking(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, alpha_T, ptherm)

      do idir = 1, ndir
        do n = 1, dust_n_species
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              dust_rho(n))>0.0d0)
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)  = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,dust_mom(idir, n))/w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, dust_rho(n))
            deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) = vgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, idir)-vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
            tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)  = 4.0d0*(dust_density(n)*dust_size(n))/ &
               (3.0d0*(one-alpha_T(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,n))*dsqrt(vt2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3) + deltav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)**2)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3, dust_rho(n)) + w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3, gas_rho_)))
          else where
            tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = bigdouble
          end where

          dtdust(n) = min(minval(tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)), dtdust(n))
        end do
      end do

      dtnew = min(minval(dust_dtpar*dtdust(:)), dtnew)

    case('linear') !linear with Deltav, for testing (see Laibe & Price 2011)
      dtdust(:) = bigdouble

      do n = 1, dust_n_species
        where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            dust_rho(n))>0.0d0)
          tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)  = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, dust_rho(n))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, gas_rho_))/ &
             (dust_K_lineardrag*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, dust_rho(n)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, gas_rho_)))
        else where
          tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = bigdouble
        end where

        dtdust(n) = min(minval(tstop(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)), dtdust(n))
      end do

      dtnew = min(minval(dust_dtpar*dtdust(:)), dtnew)
    case('usr')
      dtdust(:) = bigdouble
      call usr_dust_get_dt(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, dtdust, dx1,dx2,dx3,&
          x, dust_n_species)
      dtnew = min(minval(dust_dtpar*dtdust(:)), dtnew)
    case('none')
      ! no dust timestep
    case default
      call mpistop( "=== This dust method has not been implemented===" )
    end select

    if (dtnew < dtmin) then
      write(unitterm,*)"-------------------------------------"
      write(unitterm,*&
         )"Warning: found DUST related time step too small! dtnew=", dtnew
      write(unitterm,*)"on grid with index:", block%igrid," grid level=",&
          block%level
      write(unitterm,*)"grid corners are=",rnode(rpxmin1_, block%igrid),&
          rnode(rpxmax1_, block%igrid),rnode(rpxmin2_, block%igrid),&
          rnode(rpxmax2_, block%igrid),rnode(rpxmin3_, block%igrid),&
          rnode(rpxmax3_, block%igrid)
      write(unitterm,*)" dtdust =", dtdust(:)
      write(unitterm,*)"on processor:", mype
      write(unitterm,*)"-------------------------------------"
    endif

  end subroutine dust_get_dt

  ! Note that cmax and cmin are assumed to be initialized
  subroutine dust_get_cmax(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, cmax,&
      cmin)
    use mod_global_parameters
    use mod_variables

    integer, intent(in)                       :: ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:nw), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:number_species)
    double precision                          :: vdust(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer                                   :: n

    do n = 1, dust_n_species
      vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = get_vdust(w,&
          ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
         ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, n)

      if (present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1) = min(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1),&
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1) = max(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1),&
            vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3))
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1) = max(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1),&
            abs(vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)))
      end if
    end do
  end subroutine dust_get_cmax

  ! Note that cmax and cmin are assumed to be initialized
  subroutine dust_get_cmax_prim(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, cmax,&
      cmin)
    use mod_global_parameters

    integer, intent(in)                       :: ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:nw), x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    double precision                          :: vdust(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer                                   :: n

    do n = 1, dust_n_species
      vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = get_vdust_prim(w, ixImin1,ixImin2,ixImin3,ixImax1,&
         ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
          idim, n)

      if (present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = min(cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3), vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = max(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3), vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = max(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3), abs(vdust(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)))
      end if
    end do
  end subroutine dust_get_cmax_prim

end module mod_dust
