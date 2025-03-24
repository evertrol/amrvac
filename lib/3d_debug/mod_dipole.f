!> sets up a magnetic dipole in a 3D cartesian box
!  parameters for dipole moment: 
!      strength mu_dipole (factor 4pi included)
!      polar angle in degrees theta_dipole [0 to 90]
!      azimuthal angle in degrees phi_dipole [0 to 360]
!  center of dipole at x1_dip, x2_dip, x3_dip
module mod_dipole
  implicit none
  double precision :: mu_dipole,theta_dipole,phi_dipole,x1_dip,x2_dip,x3_dip
  double precision :: mx, my, mz

contains


  subroutine dipole(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,B)
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: B(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)

    double precision :: xsincos(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),xsinsin(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),xcos(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: rr0(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       mdotr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
 
    B=0.0d0
    mx=mu_dipole*dsin(theta_dipole*dpi/180.0d0)*dcos(phi_dipole*dpi/180.d0)
    my=mu_dipole*dsin(theta_dipole*dpi/180.0d0)*dsin(phi_dipole*dpi/180.0d0)
    mz=mu_dipole*dcos(theta_dipole*dpi/180.0d0)
    rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=dsqrt((x(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)-x1_dip)**2 +(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)-x2_dip)**2 +(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)-x3_dip)**2)
    if(any(rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)<smalldouble)) call &
       mpistop("dipole center at cell center")
    !where(rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)<smalldouble) 
    !   rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=bigdouble
    !endwhere
    xsincos(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)-x1_dip)/rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    xsinsin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)-x2_dip)/rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    xcos(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   =(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)-x3_dip)/rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    mdotr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=mx*xsincos(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+my*xsinsin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)+mz*xcos(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)=(3.0d0*xsincos(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*mdotr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-mx)/rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**3
    B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)=(3.0d0*xsinsin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*mdotr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-my)/rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**3
    B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)=(3.0d0*xcos(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)   *mdotr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-mz)/rr0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**3

  end subroutine dipole

end module mod_dipole
