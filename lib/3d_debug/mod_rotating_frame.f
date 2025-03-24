!> Module for including rotating frame in (magneto)hydrodynamics simulations
!> The rotation vector is assumed to be along z direction
!> (both in cylindrical and spherical)

module mod_rotating_frame
  implicit none

  !> Rotation frequency of the frame
  double precision :: omega_frame

  !> Index of the density (in the w array)
  integer, private, parameter :: rho_ = 1
  
contains
  !> Read this module's parameters from a file
  subroutine rotating_frame_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    
    namelist /rotating_frame_list/ omega_frame
    
    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, rotating_frame_list, end=111)
111    close(unitpar)
    end do
    
  end subroutine rotating_frame_params_read
  
  !> Initialize the module
  subroutine rotating_frame_init()
    use mod_global_parameters
    integer :: nwx,idir
    
    call rotating_frame_params_read(par_files)
    
  end subroutine rotating_frame_init
  
  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine rotating_frame_add_source(qdt,dtfactor,ixImin1,ixImin2,ixImin3,&
     ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
     wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry
    use mod_physics, only: phys_energy, phys_internal_e
    use mod_comm_lib, only: mpistop
    
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, dtfactor,x(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer                         :: idir

    ! .. local ..
    double precision :: rotating_terms(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), frame_omega(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: work(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    select case (coordinate)
    case (cylindrical)

      ! S[mrad] = 2*mphi*Omegaframe + rho*r*Omegaframe**2
      rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = omega_frame**2 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,r_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,iw_rho)

      if (phi_ > 0) then
        rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) + 2.d0 * omega_frame *wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mom(phi_))
      end if

      if(local_timestep) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            iw_mom(r_)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            iw_mom(r_)) + block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*dtfactor * rotating_terms(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            iw_mom(r_)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
            iw_mom(r_)) + qdt * rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)
      endif
      ! S[mphi] = -2*mrad*Omegaframe
      if (phi_ > 0) then
        rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3) = - 2.0d0*omega_frame * wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mom(r_))
        if(local_timestep) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(phi_)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, iw_mom(phi_)) + block%dt(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)*dtfactor * &
             rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        else
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(phi_)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3, iw_mom(phi_)) + qdt * &
             rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        endif
      end if

      ! S[etot] = mrad*r*Omegaframe**2
      if (phys_energy .and. (.not.phys_internal_e)) then
        if(local_timestep) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_e) + block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3) *dtfactor * omega_frame**2 * x(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_) * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mom(r_))
        else
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_e) + qdt * omega_frame**2 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,r_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,iw_mom(r_))
        endif
      endif

    case (spherical)
       frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = omega_frame * dsin(x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))

       ! S[mrad] = 2*mphi*Omegaframe + rho*r*Omegaframe**2
       rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3) = frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)**2 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,r_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3,iw_rho)

       if (phi_ > 0) then
         rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) = rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) + 2.d0 * frame_omega(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3) * wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,iw_mom(phi_))
       end if
        if(local_timestep) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(r_)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(r_)) + block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*dtfactor * rotating_terms(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        else
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(r_)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(r_)) + qdt * rotating_terms(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        endif
       
       ! S[mtheta] = cot(theta) * S[mrad], reuse above rotating_terms
        if(local_timestep) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(2)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(2)) + block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*dtfactor * rotating_terms(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)/ tan(x(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, 2))
        else
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(2)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(2)) + qdt * rotating_terms(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)/ tan(x(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3, 2))
        endif
       ! S[mphi] = -2*Omegaframe * (mrad + cot(theta)*mtheta)
       if (phi_ > 0) then
         rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) = -2.d0*frame_omega(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3)* wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_mom(r_))- 2.d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, iw_mom(2)) * frame_omega(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3)/ tan(x(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,ixOmin3:ixOmax3, 2))
        if(local_timestep) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(3)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
              iw_mom(3)) + block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*dtfactor * rotating_terms(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,ixOmin3:ixOmax3)
        else
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_mom(3)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             iw_mom(3)) + qdt * rotating_terms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)
        endif
       end if
      

       ! S[etot] = r*Omegaframe**2 * (mrad + cot(theta)*mtheta)
       if (phys_energy .and. (.not.phys_internal_e)) then
         work(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) = frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)**2 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,r_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,iw_mom(r_))
         
         work(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) = work(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3) + frame_omega(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)**2 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,r_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, iw_mom(2))/ tan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3, 2))
        
          if(local_timestep) then
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw_e) + block%dt(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)*dtfactor* work(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
          else
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                iw_e) + qdt * work(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)
        endif
       endif

    case default
       call mpistop("Rotating frame not implemented in this geometry")
    end select
    
  end subroutine rotating_frame_add_source

end module mod_rotating_frame
