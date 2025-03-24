!> Reaction-diffusion module (physics routines)
!>
!> Multiple reaction-diffusion systems are included: the Gray-Scott model, the
!> Schnakenberg model, the Brusselator model, the diffusive logistic equation,
!> an analytical testcase from "Numerical solution of time-dependent advection-
!> diffusion-reaction equations" by Hundsdorfer & Verwer, the Oregonator model,
!> the extended Brusselator model, and the diffusive Lorenz system. See the
!> documentation of the reaction-diffusion module for more information.
!>
!> IMEX methods are also supported. The implicit system is solved by a
!> multigrid solver coupled into MPI-AMRVAC.
!>
module mod_rd_phys
  use mod_multigrid_coupling
  use mod_comm_lib, only: mpistop


  implicit none
  private

  integer, protected, public :: u_ = 1
  integer, protected, public :: v_ = 2 !< For 2 or more equations
  integer, protected, public :: w_ = 3 !< For 3 or more equations

  !> Whether particles module is added
  logical, public, protected :: rd_particles = .false.

  !> Parameter with which to multiply the reaction timestep restriction
  double precision, public, protected :: dtreacpar = 0.5d0

  !> Name of the system to be solved
  character(len=20), public, protected :: equation_name = "gray-scott"
  integer            :: number_of_species  = 2
  integer            :: equation_type      = 1
  integer, parameter :: eq_gray_scott      = 1
  integer, parameter :: eq_schnakenberg    = 2
  integer, parameter :: eq_brusselator     = 3
  integer, parameter :: eq_logistic        = 4
  integer, parameter :: eq_analyt_hunds    = 5
  integer, parameter :: eq_belousov_fn     = 6 !Field-Noyes model, or Oregonator
  integer, parameter :: eq_ext_brusselator = 7
  integer, parameter :: eq_lorenz          = 8

  !> Diffusion coefficient for first species (u)
  double precision, public, protected :: D1 = 0.05d0
  !> Diffusion coefficient for second species (v) (if applicable)
  double precision, public, protected :: D2 = 1.0d0
  !> Diffusion coefficient for third species (w) (if applicable)
  double precision, public, protected :: D3 = 1.0d0

  !> Parameter for Schnakenberg model
  double precision, public, protected :: sb_alpha = 0.1305d0
  !> Parameter for Schnakenberg model
  double precision, public, protected :: sb_beta  = 0.7695d0
  !> Parameter for Schnakenberg model
  double precision, public, protected :: sb_kappa = 100.0d0

  !> Feed rate for Gray-Scott model
  double precision, public, protected :: gs_F = 0.046d0
  !> Kill rate for Gray-Scott model
  double precision, public, protected :: gs_k = 0.063d0

  !> Parameter for Brusselator model
  double precision, public, protected :: br_A = 4.5d0
  !> Parameter for Brusselator model
  double precision, public, protected :: br_B = 8.0d0
  !> Parameter for extended Brusselator model
  double precision, public, protected :: br_C = 1.0d0
  !> Parameter for extended Brusselator model
  double precision, public, protected :: br_D = 1.0d0

  !> Parameter for logistic model (Fisher / KPP equation)
  double precision, public, protected :: lg_lambda = 1.0d0

  !> Parameter for the Field-Noyes model of the Belousov-Zhabotinsky reaction
  double precision, public, protected :: bzfn_epsilon = 1.0d0
  !> Parameter for the Field-Noyes model of the Belousov-Zhabotinsky reaction
  double precision, public, protected :: bzfn_delta   = 1.0d0
  !> Parameter for the Field-Noyes model of the Belousov-Zhabotinsky reaction
  double precision, public, protected :: bzfn_lambda  = 1.0d0
  !> Parameter for the Field-Noyes model of the Belousov-Zhabotinsky reaction
  double precision, public, protected :: bzfn_mu      = 1.0d0

  !> Parameter for Lorenz system (Rayleigh number)
  double precision, public, protected :: lor_r     = 28.0d0
  !> Parameter for Lorenz system (Prandtl number)
  double precision, public, protected :: lor_sigma = 10.0d0
  !> Parameter for Lorenz system (aspect ratio of the convection rolls)
  double precision, public, protected :: lor_b     = 8.0d0 / 3.0d0

  !> Whether to handle the explicitly handled source term in split fashion
  logical :: rd_source_split = .false.

  !> Boundary condition information for the multigrid method
  type(mg_bc_t), public :: rd_mg_bc(3, mg_num_neighbors)

  ! Public methods
  public :: rd_phys_init

contains

  subroutine rd_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rd_list/ D1, D2, D3, sb_alpha, sb_beta, sb_kappa, gs_F, gs_k,&
        br_A, br_B, br_C, br_D, lg_lambda, bzfn_epsilon, bzfn_delta,&
        bzfn_lambda, bzfn_mu, lor_r, lor_sigma, lor_b, equation_name,&
        rd_particles, rd_source_split, dtreacpar

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, rd_list, end=111)
111    close(unitpar)
    end do

    select case (equation_name)
    case ("gray-scott")
       equation_type = eq_gray_scott
       number_of_species = 2
    case ("schnakenberg")
       equation_type = eq_schnakenberg
       number_of_species = 2
    case ("brusselator")
       equation_type = eq_brusselator
       number_of_species = 2
    case ("logistic")
       equation_type = eq_logistic
       number_of_species = 1
    case ("analyt_hunds")
       equation_type = eq_analyt_hunds
       number_of_species = 1
    case ("belousov_fieldnoyes")
       equation_type = eq_belousov_fn
       number_of_species = 3
    case ("ext_brusselator")
       equation_type = eq_ext_brusselator
       number_of_species = 3
    case ("lorenz")
       equation_type = eq_lorenz
       number_of_species = 3
    case default
       call mpistop&
          ("Unknown equation_name (valid: gray-scott, schnakenberg, ...)")
    end select

  end subroutine rd_params_read

  !> Write this module's parameters to a snapshot
  subroutine rd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 0
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er
    integer                             :: idim

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
  end subroutine rd_write_info

  subroutine rd_phys_init()
    use mod_global_parameters
    use mod_physics
    use mod_multigrid_coupling
    use mod_particles, only: particles_init

    call rd_params_read(par_files)

    physics_type = "rd"
    phys_energy  = .false.
    phys_req_diagonal = .false.
    use_particles = rd_particles

    allocate(start_indices(number_species),stop_indices(number_species))
    ! set the index of the first flux variable for species 1
    start_indices(1)=1
    ! Use the first variable as a density
    u_ = var_set_fluxvar("u", "u")
    if (number_of_species >= 2) then
        v_ = var_set_fluxvar("v", "v")
    end if
    if (number_of_species >= 3) then
        w_ = var_set_fluxvar("w", "w")
    end if

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    ! Disable flux conservation near AMR boundaries, since we have no fluxes
    fix_conserve_global = .false.

    phys_get_cmax     => rd_get_cmax
    phys_get_cbounds  => rd_get_cbounds
    phys_get_flux     => rd_get_flux
    phys_to_conserved => rd_to_conserved
    phys_to_primitive => rd_to_primitive
    phys_add_source   => rd_add_source
    phys_get_dt       => rd_get_dt
    phys_write_info   => rd_write_info
    phys_check_params => rd_check_params
    phys_implicit_update   => rd_implicit_update
    phys_evaluate_implicit => rd_evaluate_implicit

    ! Initialize particles module
    if (rd_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

  end subroutine rd_phys_init

  subroutine rd_check_params
    use mod_global_parameters
    integer :: n, i, iw, species_list(number_of_species)

    if (any(flux_method /= fs_source)) then
       ! there are no fluxes, only source terms in reaction-diffusion
       call mpistop("mod_rd requires flux_scheme = source")
    end if

    if (use_imex_scheme) then
       use_multigrid = .true.
       select case(number_of_species)
       case(1); species_list = [u_]
       case(2); species_list = [u_, v_]
       case(3); species_list = [u_, v_, w_]
       end select

       do i = 1, number_of_species
          iw = species_list(i)

          ! Set boundary conditions for the multigrid solver
          do n = 1, 2*ndim
             select case (typeboundary(iw, n))
             case (bc_symm)
                ! d/dx u = 0
                rd_mg_bc(i, n)%bc_type = mg_bc_neumann
                rd_mg_bc(i, n)%bc_value = 0.0_dp
             case (bc_asymm)
                ! u = 0
                rd_mg_bc(i, n)%bc_type = mg_bc_dirichlet
                rd_mg_bc(i, n)%bc_value = 0.0_dp
             case (bc_cont)
                ! d/dx u = 0
                rd_mg_bc(i, n)%bc_type = mg_bc_neumann
                rd_mg_bc(i, n)%bc_value = 0.0_dp
             case (bc_periodic)
                ! Nothing to do here
             case (bc_special)
                if (.not. associated(rd_mg_bc(i, n)%boundary_cond)) then
                   write(*, "(A,I0,A,I0,A)") "typeboundary(", iw, ",", n, &
                        ") is 'special', but the corresponding method " // &
                        "rd_mg_bc(i, n)%boundary_cond is not set"
                   call mpistop("rd_mg_bc(i, n)%boundary_cond not set")
                end if
             case default
                write(*,*) "rd_check_params warning: unknown boundary type"
                rd_mg_bc(i, n)%bc_type = mg_bc_dirichlet
                rd_mg_bc(i, n)%bc_value = 0.0_dp
             end select
          end do
       end do
    end if

  end subroutine rd_check_params

  subroutine rd_to_conserved(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:3)

    ! Do nothing (primitive and conservative are equal for rd module)
  end subroutine rd_to_conserved

  subroutine rd_to_primitive(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:3)

    ! Do nothing (primitive and conservative are equal for rd module)
  end subroutine rd_to_primitive

  subroutine rd_get_cmax(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,&
       ixOmax3, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:3)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
  end subroutine rd_get_cmax

  subroutine rd_get_cbounds(wLC, wRC, wLp, wRp, x, ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
      idim,Hspeed, cmax, cmin)
    use mod_global_parameters
    use mod_variables
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:3)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:number_species)
    double precision, intent(in)    :: Hspeed(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:number_species)

    if (present(cmin)) then
       cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) = 0.0d0
       cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) = 0.0d0
    else
       cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) = 0.0d0
    end if

  end subroutine rd_get_cbounds

  subroutine rd_get_dt(w, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, dtnew, dx1,dx2,dx3, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, 1:3)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(inout) :: dtnew
    double precision                :: maxrate
    double precision                :: maxD

    ! dt < dx^2 / (2 * ndim * diffusion_coeff)
    ! use dtdiffpar < 1 for explicit and > 1 for imex/split
    maxD = D1
    if (number_of_species >= 2) then
        maxD = max(maxD, D2)
    end if
    if (number_of_species >= 3) then
        maxD = max(maxD, D3)
    end if
    dtnew = dtdiffpar * minval([ dx1,dx2,dx3 ])**2 / (2 * ndim * maxD)

    ! Estimate time step for reactions
    select case (equation_type)
    case (eq_gray_scott)
       maxrate = max(maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           v_))**2 + gs_F, maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, v_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, u_)) - gs_F - gs_k)
    case (eq_schnakenberg)
       maxrate = max(maxval(abs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, v_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, u_) - 1)), maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, u_))**2)
    case (eq_brusselator)
       maxrate = max( maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           u_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           v_) - (br_B+1)), maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, u_)**2) )
    case (eq_ext_brusselator)
       maxrate = max( maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           u_)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           v_) - (br_B+1)) + br_C, maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, u_)**2) )
       maxrate = max(maxrate, br_D)
    case (eq_logistic)
       maxrate = lg_lambda*maxval(abs(1 - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, u_))) !abs for safety, normally u < 1
    case (eq_analyt_hunds)
       maxrate = maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           u_)*abs(1 - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           u_))) / D1
    case (eq_belousov_fn)
       maxrate = max(maxval(abs(1.0d0 - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, w_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, u_))) / bzfn_epsilon,&
           maxval(bzfn_lambda + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3, u_)) / bzfn_delta )
    case (eq_lorenz)
       ! det(J) = sigma(b(r-1) + x*(x*+y*))
       maxrate = max(lor_sigma, 1.0d0, lor_b)
    case default
       maxrate = one
       call mpistop("Unknown equation type")
    end select

    dtnew = min(dtnew, dtreacpar / maxrate)

  end subroutine rd_get_dt

  ! There is nothing to add to the transport flux in the transport equation
  subroutine rd_get_flux(wC, w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:3)
    double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nwflux)

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3, :) = 0.0d0
  end subroutine rd_get_flux

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine rd_add_source(qdt,dtfactor,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
     wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, dtfactor
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw),wCTprim(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision                :: lpl_u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), lpl_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), lpl_w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    ! here we add the reaction terms (always) and the diffusion if no imex is used
    if (qsourcesplit .eqv. rd_source_split) then
       if (.not.use_imex_scheme) then
          call rd_laplacian(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
              ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
              wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, u_), lpl_u)
          if (number_of_species >= 2) then
             call rd_laplacian(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
                 ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
                 wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, v_),&
                 lpl_v)
          end if
          if (number_of_species >= 3) then
             call rd_laplacian(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
                 ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
                 wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, w_),&
                 lpl_w)
          end if
       else
          ! for all IMEX scheme variants: only add the reactions
          lpl_u = 0.0d0
          lpl_v = 0.0d0
          lpl_w = 0.0d0
       end if

       select case (equation_type)
       case (eq_gray_scott)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + qdt * (D1 * lpl_u - wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, u_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, v_)**2 + gs_F * (1 - wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_)))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) + qdt * (D2 * lpl_v + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, u_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, v_)**2 - (gs_F + gs_k) * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, v_))
       case (eq_schnakenberg)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + qdt * (D1 * lpl_u + sb_kappa * (sb_alpha - &
             wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_)**2 * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_)))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) + qdt * (D2 * lpl_v + sb_kappa * (sb_beta - &
             wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_)**2 * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_)))
       case (eq_brusselator)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + qdt * (D1 * lpl_u + br_A - (br_B + 1) * &
             wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_)**2 * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) + qdt * (D2 * lpl_v + br_B * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_) - wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_)**2 * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, v_))
       case (eq_ext_brusselator)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + qdt * (D1 * lpl_u + br_A - (br_B + 1) * &
             wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_)**2 * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) - br_C * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + br_D * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) + qdt * (D2 * lpl_v + br_B * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_) - wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_)**2 * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, v_))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_) + qdt * (D3 * lpl_w + br_C * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_) - br_D * w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, w_))
       case (eq_logistic)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + qdt * (D1 * lpl_u + lg_lambda * w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_) * (1 - w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_)))
       case (eq_analyt_hunds)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + qdt * (D1 * lpl_u + 1.0d0/D1 * w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_)**2 * (1 - w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_)))
       case (eq_belousov_fn)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + qdt * (D1 * lpl_u + 1.0d0/bzfn_epsilon * (bzfn_lambda * &
             wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) - wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) - wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_)**2))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) + qdt * (D2 * lpl_v + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, u_) - wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, v_))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_) + qdt * (D3 * lpl_w + 1.0d0/bzfn_delta * (-bzfn_lambda * &
             wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_) - wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_) + bzfn_mu * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, v_)))
       case (eq_lorenz)
          ! xdot = sigma.(y-x)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              u_) + qdt * (D1 * lpl_u + lor_sigma * (wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, v_) - wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_)))
          ! ydot = r.x - y - x.z
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              v_) + qdt * (D2 * lpl_v + lor_r * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_) - wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, v_) - wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, u_)*wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, w_))
          ! zdot = x.y - b.z
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              w_) + qdt * (D3 * lpl_w + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, u_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, v_) - lor_b * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, w_))
       case default
          call mpistop("Unknown equation type")
       end select

       ! enforce getbc call after source addition
       active = .true.
    end if

  end subroutine rd_add_source

  !> Compute the Laplacian using a standard second order scheme. For now this
  !> method only works in slab geometries. Requires one ghost cell only.
  subroutine rd_laplacian(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,var,lpl)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)  :: var(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out) :: lpl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer                       :: idir, jxOmin1,jxOmin2,jxOmin3,jxOmax1,&
       jxOmax2,jxOmax3, hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3
    double precision              :: h_inv2

    if (slab) then
       lpl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0.0d0
       do idir = 1, ndim
          hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
          hxOmin3=ixOmin3-kr(idir,3);hxOmax1=ixOmax1-kr(idir,1)
          hxOmax2=ixOmax2-kr(idir,2);hxOmax3=ixOmax3-kr(idir,3);
          jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
          jxOmin3=ixOmin3+kr(idir,3);jxOmax1=ixOmax1+kr(idir,1)
          jxOmax2=ixOmax2+kr(idir,2);jxOmax3=ixOmax3+kr(idir,3);
          h_inv2 = 1/dxlevel(idir)**2
          lpl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) = lpl(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + h_inv2 * (var(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
             jxOmin3:jxOmax3) - 2 * var(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) + var(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             hxOmin3:hxOmax3))
       end do
    else
       call mpistop("rd_laplacian not implemented in this geometry")
    end if

  end subroutine rd_laplacian

  subroutine put_laplacians_onegrid(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)

    double precision                :: lpl_u(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), lpl_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3), lpl_w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    call rd_laplacian(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3, u_), lpl_u)
    if (number_of_species >= 2) then
       call rd_laplacian(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w(ixImin1:ixImax1,&
          ixImin2:ixImax2,ixImin3:ixImax3, v_), lpl_v)
    end if
    if (number_of_species >= 3) then
       call rd_laplacian(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w(ixImin1:ixImax1,&
          ixImin2:ixImax2,ixImin3:ixImax3, w_), lpl_w)
    end if

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,u_) = D1*lpl_u
    if (number_of_species >= 2) then
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v_) = D2*lpl_v
    end if
    if (number_of_species >= 3) then
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,w_) = D3*lpl_w
    end if

  end subroutine put_laplacians_onegrid
  
  !> inplace update of psa==>F_im(psa)
  subroutine rd_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: iigrid, igrid, level
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3

    !ixO^L=ixG^LL^LSUB1;
    ixOmin1=ixMlo1;ixOmin2=ixMlo2;ixOmin3=ixMlo3;ixOmax1=ixMhi1
    ixOmax2=ixMhi2;ixOmax3=ixMhi3;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
       dxlevel(3)=rnode(rpdx3_,igrid);
       call put_laplacians_onegrid(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,psa(igrid)%w)
    end do
    !$OMP END PARALLEL DO

  end subroutine rd_evaluate_implicit

  !> Implicit solve of psa=psb+dtfactor*dt*F_im(psa)
  subroutine rd_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_forest

    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer                      :: n
    double precision             :: res, max_residual, lambda

    integer                        :: iw_to,iw_from
    integer                        :: iigrid, igrid, id
    integer                        :: nc, lvl
    type(tree_node), pointer       :: pnode
    real(dp)                       :: fac

    ! Avoid setting a very restrictive limit to the residual when the time step
    ! is small (as the operator is ~ 1/(D * qdt))
    if (qdt < dtmin) then
        if(mype==0)then
            print *,'skipping implicit solve: dt too small!'
            print *,'Currently at time=',global_time,' time step=',qdt,' dtmin=',dtmin
        endif
        return
    endif
    max_residual = 1d-7/qdt

    mg%operator_type = mg_helmholtz
    call mg_set_methods(mg)

    if (.not. mg%is_allocated) call mpistop(&
       "multigrid tree not allocated yet")

    ! First handle the u variable ***************************************
    lambda           = 1/(dtfactor * qdt * D1)
    call helmholtz_set_lambda(lambda)
    mg%bc(:, mg_iphi) = rd_mg_bc(1, :)

    call mg_copy_to_tree(u_, mg_irhs, factor=-lambda, state_from=psb)
    call mg_copy_to_tree(u_, mg_iphi, state_from=psb)

    call mg_fas_fmg(mg, .true., max_res=res)
    do n = 1, 10
       call mg_fas_vcycle(mg, max_res=res)
       if (res < max_residual) exit
    end do

    call mg_copy_from_tree_gc(mg_iphi, u_, state_to=psa)
    ! Done with the u variable ***************************************

    ! Next handle the v variable ***************************************
    if (number_of_species >= 2) then
       lambda = 1/(dtfactor * qdt * D2)
       call helmholtz_set_lambda(lambda)
       mg%bc(:, mg_iphi) = rd_mg_bc(2, :)

       call mg_copy_to_tree(v_, mg_irhs, factor=-lambda, state_from=psb)
       call mg_copy_to_tree(v_, mg_iphi, state_from=psb)

       call mg_fas_fmg(mg, .true., max_res=res)
       do n = 1, 10
          call mg_fas_vcycle(mg, max_res=res)
          if (res < max_residual) exit
       end do

       call mg_copy_from_tree_gc(mg_iphi, v_, state_to=psa)
    end if
    ! Done with the v variable ***************************************

    ! Next handle the w variable ***************************************
    if (number_of_species >= 3) then
       lambda = 1/(dtfactor * qdt * D3)
       call helmholtz_set_lambda(lambda)

       call mg_copy_to_tree(w_, mg_irhs, factor=-lambda, state_from=psb)
       call mg_copy_to_tree(w_, mg_iphi, state_from=psb)

       call mg_fas_fmg(mg, .true., max_res=res)
       do n = 1, 10
          call mg_fas_vcycle(mg, max_res=res)
          if (res < max_residual) exit
       end do

       call mg_copy_from_tree_gc(mg_iphi, w_, state_to=psa)
    end if
    ! Done with the w variable ***************************************

  end subroutine rd_implicit_update

end module mod_rd_phys
