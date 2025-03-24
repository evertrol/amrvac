!>   The module add viscous source terms and check time step
!>
!>   Viscous forces in the momentum equations:
!>   d m_i/dt +=  - div (vc_mu * PI)
!>   !! Viscous work in the energy equation:
!>   !! de/dt    += - div (v . vc_mu * PI)
!>   where the PI stress tensor is
!>   PI_i,j = - (dv_j/dx_i + dv_i/dx_j) + (2/3)*Sum_k dv_k/dx_k
!>   where vc_mu is the dynamic viscosity coefficient (g cm^-1 s^-1).
module mod_viscosity
  use mod_comm_lib, only: mpistop
  implicit none

  !> Viscosity coefficient
  double precision, public :: vc_mu = 1.d0

  !> fourth order
  logical :: vc_4th_order = .false.

  !> source split or not
  logical :: vc_split= .false.

  !> whether to compute the viscous terms as
  !> fluxes (ie in the div on the LHS), or not (by default)
  logical :: viscInDiv= .false.

  !> Index of the density (in the w array)
  integer, private, parameter              :: rho_ = 1

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, private, protected              :: e_

  ! Public methods
  public :: visc_get_flux_prim

contains
  !> Read this module"s parameters from a file
  subroutine vc_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /vc_list/ vc_mu, vc_4th_order, vc_split, viscInDiv

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, vc_list, end=111)
111    close(unitpar)
    end do

  end subroutine vc_params_read

  !> Initialize the module
  subroutine viscosity_init(phys_wider_stencil,phys_req_diagonal)
    use mod_global_parameters
    integer, intent(inout) :: phys_wider_stencil
    logical, intent(inout) :: phys_req_diagonal
    integer :: nwx,idir

    call vc_params_read(par_files)

    if(vc_split) any_source_split=.true.

    ! Determine flux variables
    nwx = 1                  ! rho (density)

    allocate(mom(ndir))
    do idir = 1, ndir
       nwx    = nwx + 1
       mom(idir) = nwx       ! momentum density
    end do

    nwx = nwx + 1
    e_     = nwx          ! energy density

    if (viscInDiv) then
      ! to compute the derivatives from left and right upwinded values
      phys_wider_stencil = 1
      phys_req_diagonal = .true.  ! viscInDiv
    end if

  end subroutine viscosity_init

  subroutine viscosity_add_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x,energy,&
     qsourcesplit,active)
  ! Add viscosity source in isotropic Newtonian fluids to w within ixO
  ! neglecting bulk viscosity
  ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    integer:: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,idir,jdir,iw
    double precision:: lambda(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ndir,ndir),tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),v(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndir),vlambda(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,ndir)

    if (viscInDiv) return

    if(qsourcesplit .eqv. vc_split) then
      active = .true.
      ! standard case, textbook viscosity
      ! Calculating viscosity sources
      if(.not.vc_4th_order) then
        ! involves second derivatives, two extra layers
        ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmin3=ixOmin3-2;ixmax1=ixOmax1+2
        ixmax2=ixOmax2+2;ixmax3=ixOmax3+2;
        if( ixImin1>ixmin1 .or. ixImax1<ixmax1.or. ixImin2>ixmin2 .or. &
           ixImax2<ixmax2.or. ixImin3>ixmin3 .or. ixImax3<ixmax3)call &
           mpistop("error for viscous source addition, 2 layers needed")
        ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
        ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;
      else
        ! involves second derivatives, four extra layers
        ixmin1=ixOmin1-4;ixmin2=ixOmin2-4;ixmin3=ixOmin3-4;ixmax1=ixOmax1+4
        ixmax2=ixOmax2+4;ixmax3=ixOmax3+4;
        if( ixImin1>ixmin1 .or. ixImax1<ixmax1.or. ixImin2>ixmin2 .or. &
           ixImax2<ixmax2.or. ixImin3>ixmin3 .or. ixImax3<ixmax3)&
          call mpistop("error for viscous source addition"//&
          "requested fourth order gradients: 4 layers needed")
        ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmin3=ixOmin3-2;ixmax1=ixOmax1+2
        ixmax2=ixOmax2+2;ixmax3=ixOmax3+2;
      end if

      ! get velocity
      do idir=1,ndir
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           idir)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           mom(idir))/wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           rho_)
      end do

      ! construct lambda tensor: lambda_ij = gradv_ij + gradv_ji
      ! initialize
      lambda=zero

      !next construct
      do idim=1,ndim; do idir=1,ndir
      ! Calculate velocity gradient tensor within ixL: gradv= grad v,
      ! thus gradv_ij=d_j v_i
        tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,idir)
        ! Correction for Christoffel terms in non-cartesian
        if (coordinate==cylindrical .and. idim==r_  .and. idir==phi_  ) &
           tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3) = tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)/x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)
        if (coordinate==spherical) then
          if     (idim==r_  .and. (idir==2 .or. idir==phi_)) then
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3) = tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3)/x(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3,1)

          elseif (idim==2  .and. idir==phi_) then
            tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3)/dsin(x(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3,2))

          endif
        endif
        call gradient(tmp,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idim,tmp2)
        ! Correction for Christoffel terms in non-cartesian
        if (coordinate==cylindrical .and. idim==r_  .and. idir==phi_  ) &
           tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp2(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)*x(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3,1)
        if (coordinate==cylindrical .and. idim==phi_ .and. idir==phi_ ) &
           tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp2(ixmin1:ixmax1,&
           ixmin2:ixmax2,ixmin3:ixmax3)+v(ixmin1:ixmax1,ixmin2:ixmax2,&
           ixmin3:ixmax3,r_)/x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)
        if (coordinate==spherical) then
          if (idim==r_  .and. (idir==2 .or. idir==phi_)) then
            tmp2(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3) = tmp2(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3)*x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)

          elseif (idim==2  .and. idir==phi_ ) then
            tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp2(ixmin1:ixmax1,&
               ixmin2:ixmax2,ixmin3:ixmax3)*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,2))
          elseif (idim==2   .and. idir==2   ) then
            tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp2(ixmin1:ixmax1,&
               ixmin2:ixmax2,ixmin3:ixmax3)+v(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,r_)/x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               1)
          elseif (idim==phi_.and. idir==phi_) then
            tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp2(ixmin1:ixmax1,&
               ixmin2:ixmax2,ixmin3:ixmax3)+v(ixmin1:ixmax1,ixmin2:ixmax2,&
               ixmin3:ixmax3,r_)/x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               1)+v(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               2)/(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
               1)*dtan(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)))

          endif
        endif
        lambda(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim,&
           idir)= lambda(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idim,&
           idir)+ tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)
        lambda(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idir,&
           idim)= lambda(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idir,&
           idim)+ tmp2(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)
      enddo; enddo;

      ! Multiply lambda with viscosity coefficient and dt
      lambda(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1:ndir,&
         1:ndir)=lambda(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1:ndir,&
         1:ndir)*vc_mu*qdt

      !calculate div v term through trace action separately
      ! rq : it is safe to use the trace rather than compute the divergence
      !      since we always retrieve the divergence (even with the
      !      Christoffel terms)
      tmp=0.d0
      do idir=1,ndir
         tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp(ixmin1:ixmax1,&
            ixmin2:ixmax2,ixmin3:ixmax3)+lambda(ixmin1:ixmax1,ixmin2:ixmax2,&
            ixmin3:ixmax3,idir,idir)
      end do
      tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=tmp(ixmin1:ixmax1,&
         ixmin2:ixmax2,ixmin3:ixmax3)/3.d0

      !substract trace from diagonal elements
      do idir=1,ndir
         lambda(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idir,&
            idir)=lambda(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idir,&
            idir)-tmp(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)
      enddo

      ! dm/dt= +div(mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v * kr)
      ! hence m_j=m_j+d_i tensor_ji
      do idir=1,ndir
        do idim=1,ndim
              tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)=lambda(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3,idir,idim)
              ! Correction for divergence of a tensor
              if (coordinate==cylindrical .and. idim==r_ .and. (idir==r_ .or. &
                 idir==z_)) tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3) = tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)*x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)
              if (coordinate==cylindrical .and. idim==r_ .and. idir==phi_      &
                         ) tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3) = tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                 ixmin3:ixmax3)*x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                 1)**two
              if (coordinate==spherical) then
                if (idim==r_ .and. idir==r_                 ) &
                   tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3) = tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3)*x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                   1)**two
                if (idim==r_ .and. (idir==2 .or. idir==phi_)) &
                   tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3) = tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3)*x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                   1)**3.d0

                if (idim==2  .and. (idir==r_ .or. idir==2))   &
                   tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3) = tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3)*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3,2))
                if (idim==2  .and. idir==phi_               ) &
                   tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3) = tmp(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3)*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,&
                   ixmin3:ixmax3,2))**two

              endif
              call gradient(tmp,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
                 ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,&
                 tmp2)
              ! Correction for divergence of a tensor
              if (coordinate==cylindrical .and. idim==r_ .and. (idir==r_ .or. &
                 idir==z_)) tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,1)
              if (coordinate==cylindrical .and. idim==r_ .and. idir==phi_      &
                         ) tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3,1)**two)
              if (coordinate==spherical) then
                if (idim==r_ .and. idir==r_                 ) &
                   tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,1)**two)
                if (idim==r_ .and. (idir==2 .or. idir==phi_)) &
                   tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,1)**3.d0)

                if (idim==2  .and. (idir==r_ .or. idir==2))   &
                   tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)/(dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,2)))
                if (idim==2  .and. idir==phi_               ) &
                   tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3) = tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)/(dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,2))**two)

              endif
              w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                 mom(idir))+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 ixOmin3:ixOmax3)
        enddo
        ! Correction for geometrical terms in the div of a tensor
        if (coordinate==cylindrical .and. idir==r_  ) w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))-lambda(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_,phi_)/x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
        if (coordinate==spherical   .and. idir==r_  ) w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))-(lambda(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,2,2)+lambda(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_,phi_))/x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)

        if (coordinate==spherical   .and. idir==2   ) w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))-lambda(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,phi_,phi_)/(x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)/dtan(x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))

      end do

      if(energy) then
        ! de/dt= +div(v.dot.[mu*[d_j v_i+d_i v_j]-(2*mu/3)* div v *kr])
        ! thus e=e+d_i v_j tensor_ji
        vlambda=0.d0
        do idim=1,ndim
          do idir=1,ndir
             vlambda(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                idim)=vlambda(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                idim)+v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                idir)*lambda(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
                idir,idim)
          end do
        end do
        call divvector(vlambda,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp2)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           e_)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end if
    end if

  end subroutine viscosity_add_source

  subroutine viscosity_get_dt(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,&
     x)
    ! Check diffusion time limit for dt < dtdiffpar * dx**2 / (mu/rho)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision:: dtdiff_visc, dxinv2(1:ndim)
    integer:: idim

    ! Calculate the kinematic viscosity tmp=mu/rho

    tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=vc_mu/w(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

    dxinv2(1)=one/dx1**2;dxinv2(2)=one/dx2**2;dxinv2(3)=one/dx3**2;
    do idim=1,ndim
       dtdiff_visc=dtdiffpar/maxval(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)*dxinv2(idim))
       ! limit the time step
       dtnew=min(dtnew,dtdiff_visc)
    enddo

  end subroutine viscosity_get_dt

  ! viscInDiv
  ! Get the viscous stress tensor terms in the idim direction
  ! Beware : a priori, won't work for ndir /= ndim
  ! Rq : we work with primitive w variables here
  ! Rq : ixO^L is already extended by 1 unit in the direction we work on

  subroutine visc_get_flux_prim(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim, f,&
      energy)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:3)
    double precision, intent(inout) :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, nwflux)
    logical, intent(in) :: energy
    integer                         :: idir, i
    double precision :: v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir)

    double precision                :: divergence(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    double precision:: lambda(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ndir) !, tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) !gradV(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,ndir,ndir)

    if (.not. viscInDiv) return

    do i=1,ndir
     v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,i)=w(ixImin1:ixImax1,&
        ixImin2:ixImax2,ixImin3:ixImax3,i+1)
    enddo
    call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divergence)

    call get_crossgrad(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w,idim,lambda)
    lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idim) = lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       idim) - (2.d0/3.d0) * divergence(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    ! Compute the idim-th row of the viscous stress tensor
    do idir = 1, ndir
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(idir)) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          mom(idir)) - vc_mu * lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,idir)
      if (energy) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          e_) = f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          e_) - vc_mu * lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir) * v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idir)
    enddo

  end subroutine visc_get_flux_prim

  ! Compute the cross term ( d_i v_j + d_j v_i in Cartesian BUT NOT IN
  ! CYLINDRICAL AND SPHERICAL )
  subroutine get_crossgrad(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w,idim,cross)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(out)   :: cross(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndir)
    integer :: idir
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    if (ndir/=ndim) call mpistop&
       ("This formula are probably wrong for ndim/=ndir")
    ! Beware also, we work w/ the angle as the 3rd component in cylindrical
    ! and the colatitude as the 2nd one in spherical
    cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:)=zero
    tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=zero
    select case(coordinate)
    case (Cartesian,Cartesian_stretched)
      call cart_cross_grad(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w,idim,cross)
    case (cylindrical)
      if (idim==1) then
        ! for rr and rz
        call cart_cross_grad(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w,idim,cross)
        ! then we overwrite rth w/ the correct expression
        
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th (rq : already contains 1/r)
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2))/x(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,1) !v_th / r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp) !d_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)*x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)
       
      elseif (idim==2) then
        ! thr (idem as above)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2))/x(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,1) !v_th / r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp) !d_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)*x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)
        ! thth
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2)) !v_th
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=two*tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)+two*v(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)/x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1) !+ 2 vr/r
        !thz
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(3)) !v_z
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th
       
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2)) !v_th
        
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp) !d_z
       
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        
      elseif (idim==3) then
        ! for zz and rz
        call cart_cross_grad(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w,idim,cross)
        ! then we overwrite zth w/ the correct expression
        !thz
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(3)) !v_z
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2)) !v_th
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp) !d_z
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
       
      endif
    case (spherical)
      if (idim==1) then
        ! rr (normal, simply 2 * dr vr)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp) !d_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)=two*tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        !rth
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th (rq : already contains 1/r)
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2))/x(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,1) !v_th / r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp) !d_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)*x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)
       
        
        !rph
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp) !d_phi (rq : contains 1/rsin(th))
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(3))/x(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,1) !v_phi / r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp) !d_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)*x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)
       
      elseif (idim==2) then
        ! thr
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th (rq : already contains 1/r)
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2))/x(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,1) !v_th / r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp) !d_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)*x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)
        ! thth
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2)) !v_th
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=two*tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)+two*v(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)/x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1) !+ 2 vr/r
       
        
        !thph
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2)) !v_th
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp) !d_phi (rq : contains 1/rsin(th))
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(3))/dsin(x(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,2)) !v_ph / sin(th)
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)*dsin(x(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,2))
       
        
      elseif (idim==3) then
        !phr
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp) !d_phi
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(3))/x(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,1) !v_phi / r
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp) !d_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)*x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)
        !phth
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2)) !v_th
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp) !d_phi
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(3))/dsin(x(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,2)) !v_ph / sin(th)
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           2)+tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)*dsin(x(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,2))
        !phph
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(3)) !v_ph
        call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
           ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp) !d_phi
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)=two*tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(1)) !v_r
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)+two*v(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)/x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1) !+ 2 vr/r
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,mom(2)) !v_th
        cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)=cross(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           3)+two*v(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)/(x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
           1)*dtan(x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,2))) !+ 2 vth/(rtan(th))
       
      endif
    case default
      call mpistop("Unknown geometry specified")
    end select

  end subroutine get_crossgrad

  !> yields d_i v_j + d_j v_i for a given i, OK in Cartesian and for some
  !> tensor terms in cylindrical (rr & rz) and in spherical (rr)
  subroutine cart_cross_grad(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w,idim,cross)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:3)
    double precision, intent(out)   :: cross(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,ndir)
    integer :: idir
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,mom(idim))
    do idir=1,ndir
      call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,tmp)
      cross(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    enddo
    do idir=1,ndir
      v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,mom(idir))
      call gradient(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idim,tmp)
      cross(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)=cross(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         idir)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    enddo

  end subroutine cart_cross_grad

  subroutine visc_add_source_geom(qdt, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, wCT, w, x)
    use mod_global_parameters
    use mod_geometry
    ! w and wCT conservative variables here
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    double precision :: v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir), vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        divergence(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer          :: i

    if (.not. viscInDiv) return

    select case (coordinate)
    case (cylindrical)
      ! get the velocity components
      do i=1,ndir
       v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          i)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          mom(i))/wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_)
      enddo
      ! thth tensor term - - -
        ! 1st the cross grad term

      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2) !v_th
      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp1) !d_th
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=two*(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         1)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)) !2 ( d_th v_th / r + vr/r )
        ! 2nd the divergence
      call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divergence)
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) - (2.d0/3.d0) * divergence(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      ! s[mr]=-thth/radius
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))-qdt*vc_mu*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
      ! rth tensor term - - -
      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1) !v_r
      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp1) !d_th
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2)/x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1) !v_th / r
      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp1) !d_r
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)+tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
      ! s[mphi]=+rth/radius
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2))+qdt*vc_mu*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)

    case (spherical)
      ! get the velocity components
      do i=1,ndir
       v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          i)=wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          mom(i))/wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,rho_)
      enddo
      ! thth tensor term - - -
      ! 1st the cross grad term
      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2) !v_th

      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp1) !d_th
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=two*(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)) !2 ( 1/r * d_th v_th + vr/r )
      ! 2nd the divergence
      call divvector(v,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divergence)
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) - (2.d0/3.d0) * divergence(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      ! s[mr]=-thth/radius
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))-qdt*vc_mu*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)

      ! phiphi tensor term - - -
      ! 1st the cross grad term
      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,3) !v_ph

      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp1) !d_phi
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=two*(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)+v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)+v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)*dtan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))) !2 ( 1/rsinth * d_ph v_ph + vr/r + vth/rtanth )

      ! 2nd the divergence
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) - (2.d0/3.d0) * divergence(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      ! s[mr]=-phiphi/radius
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))-qdt*vc_mu*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
      ! s[mth]=-cotanth*phiphi/radius

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2))-qdt*vc_mu*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)*dtan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))

      ! rth tensor term - - -
      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1) !v_r
      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp) !d_th (rq : already contains 1/r)
      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2)/x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1) !v_th / r
      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp1) !d_r
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)+tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
      ! s[mth]=+rth/radius
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2))+qdt*vc_mu*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
      ! rphi tensor term - - -
      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1) !v_r

      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp) !d_phi (rq : contains 1/rsin(th))

      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,3)/x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1) !v_phi / r
      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,tmp1) !d_r
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)+tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
      ! s[mphi]=+rphi/radius
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(3))+qdt*vc_mu*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
      ! phith tensor term - - -
      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2) !v_th

      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,3,tmp) !d_phi


      vv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=v(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,3)/dsin(x(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,2)) !v_ph / sin(th)
      call gradient(vv,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,2,tmp1) !d_th
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)+tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2))
      ! s[mphi]=+cotanth*phith/radius
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(3))+qdt*vc_mu*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)*dtan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))

    end select

  end subroutine visc_add_source_geom

end module mod_viscosity
