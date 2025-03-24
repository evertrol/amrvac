!> Subroutines for Roe-type Riemann solver for HD
module mod_hyperdiffusivity


  implicit none

  !all public
  !private
  !public :: hyperdiffusivity_init


contains

  subroutine hyperdiffusivity_init()
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: phys_req_diagonal

    !print*, slab_uniform !!THIS IS FALSE, why??

    !if(coordinate .ne. Cartesian .or. .not. slab_uniform) then
    if(coordinate .ne. Cartesian) then
      call mpistop&
         ("Hyperdiffusivity only implemented for Cartesian uniform grid")
    endif

    nghostcells = max(nghostcells,3)
    phys_req_diagonal = .true.

  end subroutine hyperdiffusivity_init


  subroutine hyp_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, var, idimm, nu_hyp)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, idimm
    integer, intent(out)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, intent(in)       :: var(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out)       :: nu_hyp(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    double precision        :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: hxmin1,hxmin2,hxmin3,hxmax1,hxmax2,hxmax3, hx1fmin1,hx1fmin2,&
       hx1fmin3,hx1fmax1,hx1fmax2,hx1fmax3, hx1bmin1,hx1bmin2,hx1bmin3,&
       hx1bmax1,hx1bmax2,hx1bmax3, hx2bmin1,hx2bmin2,hx2bmin3,hx2bmax1,&
       hx2bmax2,hx2bmax3

    double precision, parameter :: eps=1e-8


   hxmin1=ixImin1+2*kr(idimm,1);hxmin2=ixImin2+2*kr(idimm,2)
   hxmin3=ixImin3+2*kr(idimm,3); 
   hxmax1=ixImax1-kr(idimm,1);hxmax2=ixImax2-kr(idimm,2)
   hxmax3=ixImax3-kr(idimm,3); 

   hx1fmin1=hxmin1+kr(idimm,1);hx1fmin2=hxmin2+kr(idimm,2)
   hx1fmin3=hxmin3+kr(idimm,3);hx1fmax1=hxmax1+kr(idimm,1)
   hx1fmax2=hxmax2+kr(idimm,2);hx1fmax3=hxmax3+kr(idimm,3);  
   hx1bmin1=hxmin1-kr(idimm,1);hx1bmin2=hxmin2-kr(idimm,2)
   hx1bmin3=hxmin3-kr(idimm,3);hx1bmax1=hxmax1-kr(idimm,1)
   hx1bmax2=hxmax2-kr(idimm,2);hx1bmax3=hxmax3-kr(idimm,3);  
   hx2bmin1=hxmin1-2*kr(idimm,1);hx2bmin2=hxmin2-2*kr(idimm,2)
   hx2bmin3=hxmin3-2*kr(idimm,3);hx2bmax1=hxmax1-2*kr(idimm,1)
   hx2bmax2=hxmax2-2*kr(idimm,2);hx2bmax3=hxmax3-2*kr(idimm,3);  

  !store d3 in tmp
  tmp(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3)= abs(3d0 * (var(hxmin1:hxmax1,&
     hxmin2:hxmax2,hxmin3:hxmax3) - var(hx1bmin1:hx1bmax1,hx1bmin2:hx1bmax2,&
     hx1bmin3:hx1bmax3)) - (var(hx1fmin1:hx1fmax1,hx1fmin2:hx1fmax2,&
     hx1fmin3:hx1fmax3)-var(hx2bmin1:hx2bmax1,hx2bmin2:hx2bmax2,&
     hx2bmin3:hx2bmax3)))
  !store d1 in tmp2
  tmp2(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3)= abs((var(hxmin1:hxmax1,&
     hxmin2:hxmax2,hxmin3:hxmax3) - var(hx1bmin1:hx1bmax1,hx1bmin2:hx1bmax2,&
     hx1bmin3:hx1bmax3)))

  ixOmin1=hxmin1+kr(idimm,1);ixOmin2=hxmin2+kr(idimm,2)
  ixOmin3=hxmin3+kr(idimm,3);
  ixOmax1=hxmax1-kr(idimm,1);ixOmax2=hxmax2-kr(idimm,2)
  ixOmax3=hxmax3-kr(idimm,3);
  hx1fmin1=ixOmin1+kr(idimm,1);hx1fmin2=ixOmin2+kr(idimm,2)
  hx1fmin3=ixOmin3+kr(idimm,3);hx1fmax1=ixOmax1+kr(idimm,1)
  hx1fmax2=ixOmax2+kr(idimm,2);hx1fmax3=ixOmax3+kr(idimm,3);  
  hx1bmin1=ixOmin1-kr(idimm,1);hx1bmin2=ixOmin2-kr(idimm,2)
  hx1bmin3=ixOmin3-kr(idimm,3);hx1bmax1=ixOmax1-kr(idimm,1)
  hx1bmax2=ixOmax2-kr(idimm,2);hx1bmax3=ixOmax3-kr(idimm,3);  

  nu_hyp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3) = max(tmp(hx1bmin1:hx1bmax1,hx1bmin2:hx1bmax2,&
     hx1bmin3:hx1bmax3),tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
     tmp(hx1fmin1:hx1fmax1,hx1fmin2:hx1fmax2,&
     hx1fmin3:hx1fmax3))/max(tmp2(hx1bmin1:hx1bmax1,hx1bmin2:hx1bmax2,&
     hx1bmin3:hx1bmax3),tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
     tmp2(hx1fmin1:hx1fmax1,hx1fmin2:hx1fmax2,hx1fmin3:hx1fmax3),eps)  

  !print*, "HYP IXO ", ixO^L 

  end subroutine hyp_coeff



  subroutine div_vel_coeff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, vel, idimm, nu_vel)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, idimm
    integer, intent(out)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, intent(in)       :: vel(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    double precision, intent(out)       :: nu_vel(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: hx1fmin1,hx1fmin2,hx1fmin3,hx1fmax1,hx1fmax2,hx1fmax3, hx1bmin1,&
       hx1bmin2,hx1bmin3,hx1bmax1,hx1bmax2,hx1bmax3

    integer :: ii

   ixOmin1=ixImin1+1;ixOmin2=ixImin2+1;ixOmin3=ixImin3+1; 
   ixOmax1=ixImax1-1;ixOmax2=ixImax2-1;ixOmax3=ixImax3-1; 

   ixOmin1=ixOmin1+kr(idimm,1);ixOmin2=ixOmin2+kr(idimm,2)
   ixOmin3=ixOmin3+kr(idimm,3); 

   nu_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0 

  ! ndim should always be smaller or equal to ndir !
   do ii = 1, ndim
    hx1bmin1=ixOmin1-kr(ii,1);hx1bmin2=ixOmin2-kr(ii,2)
    hx1bmin3=ixOmin3-kr(ii,3);hx1bmax1=ixOmax1-kr(ii,1)
    hx1bmax2=ixOmax2-kr(ii,2);hx1bmax3=ixOmax3-kr(ii,3);
    if(ii==idimm) then
      nu_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = nu_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + (vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,ii) - vel(hx1bmin1:hx1bmax1,hx1bmin2:hx1bmax2,&
         hx1bmin3:hx1bmax3,ii))/dxlevel(ii) 
    else
      hx1fmin1=ixOmin1+kr(ii,1);hx1fmin2=ixOmin2+kr(ii,2)
      hx1fmin3=ixOmin3+kr(ii,3);hx1fmax1=ixOmax1+kr(ii,1)
      hx1fmax2=ixOmax2+kr(ii,2);hx1fmax3=ixOmax3+kr(ii,3);  
      nu_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = nu_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) + (vel(hx1fmin1:hx1fmax1,hx1fmin2:hx1fmax2,&
         hx1fmin3:hx1fmax3,ii) - vel(hx1bmin1:hx1bmax1,hx1bmin2:hx1bmax2,&
         hx1bmin3:hx1bmax3,ii))/(2d0*dxlevel(ii)) 
    endif

    where(nu_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) < 0d0)
      nu_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = -nu_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)
    elsewhere
      nu_vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0d0
    endwhere

   enddo 

  !print*, "DIV_VEL IXO ", ixO^L 

  end subroutine div_vel_coeff



  !var has cell centered values
  !nu_hyper is defined at the interfaces
  subroutine second_same_deriv(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nu_hyper, var, idimm,&
      res)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, idimm
    integer, intent(out)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, intent(in)       :: nu_hyper(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3),var(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out)       :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    integer :: hxfmin1,hxfmin2,hxfmin3,hxfmax1,hxfmax2,hxfmax3, hxbmin1,&
       hxbmin2,hxbmin3,hxbmax1,hxbmax2,hxbmax3

    ixOmin1=ixImin1+3;ixOmin2=ixImin2+3;ixOmin3=ixImin3+3;
    ixOmax1=ixImax1-3;ixOmax2=ixImax2-3;ixOmax3=ixImax3-3;

    hxfmin1=ixOmin1+kr(idimm,1);hxfmin2=ixOmin2+kr(idimm,2)
    hxfmin3=ixOmin3+kr(idimm,3);hxfmax1=ixOmax1+kr(idimm,1)
    hxfmax2=ixOmax2+kr(idimm,2);hxfmax3=ixOmax3+kr(idimm,3);
    hxbmin1=ixOmin1-kr(idimm,1);hxbmin2=ixOmin2-kr(idimm,2)
    hxbmin3=ixOmin3-kr(idimm,3);hxbmax1=ixOmax1-kr(idimm,1)
    hxbmax2=ixOmax2-kr(idimm,2);hxbmax3=ixOmax3-kr(idimm,3);

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       1d0/(dxlevel(idimm)**2)*(nu_hyper(hxfmin1:hxfmax1,hxfmin2:hxfmax2,&
       hxfmin3:hxfmax3) * (var(hxfmin1:hxfmax1,hxfmin2:hxfmax2,&
       hxfmin3:hxfmax3)-var(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))-nu_hyper(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * (var(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-var(hxbmin1:hxbmax1,hxbmin2:hxbmax2,hxbmin3:hxbmax3)))
   !print*, "SECOND SAME DERIV IXO ", ixO^L 

  end subroutine second_same_deriv

  !var has cell centered values
  !var2 has cell centered values
  !nu_hyper is defined at the interfaces
  subroutine second_same_deriv2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nu_hyper, var2,&
      var, idimm, res)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, idimm
    integer, intent(out)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, intent(in)       :: nu_hyper(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3),var(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), var2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(out)       :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    integer :: hxfmin1,hxfmin2,hxfmin3,hxfmax1,hxfmax2,hxfmax3, hxbmin1,&
       hxbmin2,hxbmin3,hxbmax1,hxbmax2,hxbmax3

    ixOmin1=ixImin1+3;ixOmin2=ixImin2+3;ixOmin3=ixImin3+3;
    ixOmax1=ixImax1-3;ixOmax2=ixImax2-3;ixOmax3=ixImax3-3;

    hxfmin1=ixOmin1+kr(idimm,1);hxfmin2=ixOmin2+kr(idimm,2)
    hxfmin3=ixOmin3+kr(idimm,3);hxfmax1=ixOmax1+kr(idimm,1)
    hxfmax2=ixOmax2+kr(idimm,2);hxfmax3=ixOmax3+kr(idimm,3);
    hxbmin1=ixOmin1-kr(idimm,1);hxbmin2=ixOmin2-kr(idimm,2)
    hxbmin3=ixOmin3-kr(idimm,3);hxbmax1=ixOmax1-kr(idimm,1)
    hxbmax2=ixOmax2-kr(idimm,2);hxbmax3=ixOmax3-kr(idimm,3);

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       1d0/(2d0*dxlevel(idimm)**2)*(nu_hyper(hxfmin1:hxfmax1,hxfmin2:hxfmax2,&
       hxfmin3:hxfmax3) * (var2(hxfmin1:hxfmax1,hxfmin2:hxfmax2,&
       hxfmin3:hxfmax3)+var2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) *  (var(hxfmin1:hxfmax1,hxfmin2:hxfmax2,&
       hxfmin3:hxfmax3)-var(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))-nu_hyper(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) * (var2(hxbmin1:hxbmax1,hxbmin2:hxbmax2,&
       hxbmin3:hxbmax3)+var2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)) * (var(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)-var(hxbmin1:hxbmax1,hxbmin2:hxbmax2,hxbmin3:hxbmax3)))
    !print*, "SECOND SAME DERIV2 IXO ", ixO^L 

  end subroutine second_same_deriv2

  !idimm inner derivative, idimm2 outer
  ! deriv_idimm2(nu * deriv_idimm (u) )
  !var has cell centered values
  !nu_hyper is defined at the interfaces
  subroutine second_cross_deriv(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nu_hyper, var,&
      idimm, idimm2, res)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, idimm,idimm2
    integer, intent(out)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, intent(in)       :: nu_hyper(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3),var(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision, intent(out)       :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    integer :: hxfimin1,hxfimin2,hxfimin3,hxfimax1,hxfimax2,hxfimax3, hxbimin1,&
       hxbimin2,hxbimin3,hxbimax1,hxbimax2,hxbimax3, hxfifjmin1,hxfifjmin2,&
       hxfifjmin3,hxfifjmax1,hxfifjmax2,hxfifjmax3, hxbibjmin1,hxbibjmin2,&
       hxbibjmin3,hxbibjmax1,hxbibjmax2,hxbibjmax3, hxfibjmin1,hxfibjmin2,&
       hxfibjmin3,hxfibjmax1,hxfibjmax2,hxfibjmax3, hxbifjmin1,hxbifjmin2,&
       hxbifjmin3,hxbifjmax1,hxbifjmax2,hxbifjmax3

    ixOmin1=ixImin1+3;ixOmin2=ixImin2+3;ixOmin3=ixImin3+3;
    ixOmax1=ixImax1-3;ixOmax2=ixImax2-3;ixOmax3=ixImax3-3;

    hxfimin1=ixOmin1+kr(idimm2,1);hxfimin2=ixOmin2+kr(idimm2,2)
    hxfimin3=ixOmin3+kr(idimm2,3);hxfimax1=ixOmax1+kr(idimm2,1)
    hxfimax2=ixOmax2+kr(idimm2,2);hxfimax3=ixOmax3+kr(idimm2,3);
    hxbimin1=ixOmin1-kr(idimm2,1);hxbimin2=ixOmin2-kr(idimm2,2)
    hxbimin3=ixOmin3-kr(idimm2,3);hxbimax1=ixOmax1-kr(idimm2,1)
    hxbimax2=ixOmax2-kr(idimm2,2);hxbimax3=ixOmax3-kr(idimm2,3);

    hxfifjmin1=hxfimin1+kr(idimm,1);hxfifjmin2=hxfimin2+kr(idimm,2)
    hxfifjmin3=hxfimin3+kr(idimm,3);hxfifjmax1=hxfimax1+kr(idimm,1)
    hxfifjmax2=hxfimax2+kr(idimm,2);hxfifjmax3=hxfimax3+kr(idimm,3);
    hxfibjmin1=hxfimin1-kr(idimm,1);hxfibjmin2=hxfimin2-kr(idimm,2)
    hxfibjmin3=hxfimin3-kr(idimm,3);hxfibjmax1=hxfimax1-kr(idimm,1)
    hxfibjmax2=hxfimax2-kr(idimm,2);hxfibjmax3=hxfimax3-kr(idimm,3);

    hxbifjmin1=hxbimin1+kr(idimm,1);hxbifjmin2=hxbimin2+kr(idimm,2)
    hxbifjmin3=hxbimin3+kr(idimm,3);hxbifjmax1=hxbimax1+kr(idimm,1)
    hxbifjmax2=hxbimax2+kr(idimm,2);hxbifjmax3=hxbimax3+kr(idimm,3);
    hxbibjmin1=hxbimin1-kr(idimm,1);hxbibjmin2=hxbimin2-kr(idimm,2)
    hxbibjmin3=hxbimin3-kr(idimm,3);hxbibjmax1=hxbimax1-kr(idimm,1)
    hxbibjmax2=hxbimax2-kr(idimm,2);hxbibjmax3=hxbimax3-kr(idimm,3);

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       1d0/(8d0*dxlevel(idimm) * dxlevel(idimm2))*((nu_hyper(&
       hxfifjmin1:hxfifjmax1,hxfifjmin2:hxfifjmax2,&
       hxfifjmin3:hxfifjmax3) + nu_hyper(hxfimin1:hxfimax1,hxfimin2:hxfimax2,&
       hxfimin3:hxfimax3)) * (var(hxfifjmin1:hxfifjmax1,hxfifjmin2:hxfifjmax2,&
       hxfifjmin3:hxfifjmax3)-var(hxfibjmin1:hxfibjmax1,hxfibjmin2:hxfibjmax2,&
       hxfibjmin3:hxfibjmax3))-(nu_hyper(hxbifjmin1:hxbifjmax1,&
       hxbifjmin2:hxbifjmax2,hxbifjmin3:hxbifjmax3) + &
       nu_hyper(hxbimin1:hxbimax1,hxbimin2:hxbimax2,&
       hxbimin3:hxbimax3)) * (var(hxbifjmin1:hxbifjmax1,hxbifjmin2:hxbifjmax2,&
       hxbifjmin3:hxbifjmax3)-var(hxbibjmin1:hxbibjmax1,hxbibjmin2:hxbibjmax2,&
       hxbibjmin3:hxbibjmax3)))

    !print*, "SECOND CROSS DERIV IXO ", ixO^L 

  end subroutine second_cross_deriv


  !idimm inner derivative, idimm2 outer
  ! deriv_idimm2(nu * var2 * deriv_idimm (u) )
  !var has cell centered values
  !var2 has cell centered values
  !nu_hyper is defined at the interfaces (with numbering index left center): center_i-1 interface_i center_i
  subroutine second_cross_deriv2(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, nu_hyper, var2,&
      var, idimm, idimm2, res)
    use mod_global_parameters
    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, idimm,idimm2
    integer, intent(out)                :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3
    double precision, intent(in)       :: nu_hyper(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3),var(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),var2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, intent(out)       :: res(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    integer :: hxfimin1,hxfimin2,hxfimin3,hxfimax1,hxfimax2,hxfimax3, hxbimin1,&
       hxbimin2,hxbimin3,hxbimax1,hxbimax2,hxbimax3, hxfifjmin1,hxfifjmin2,&
       hxfifjmin3,hxfifjmax1,hxfifjmax2,hxfifjmax3, hxbibjmin1,hxbibjmin2,&
       hxbibjmin3,hxbibjmax1,hxbibjmax2,hxbibjmax3, hxfibjmin1,hxfibjmin2,&
       hxfibjmin3,hxfibjmax1,hxfibjmax2,hxfibjmax3, hxbifjmin1,hxbifjmin2,&
       hxbifjmin3,hxbifjmax1,hxbifjmax2,hxbifjmax3

    ixOmin1=ixImin1+3;ixOmin2=ixImin2+3;ixOmin3=ixImin3+3;
    ixOmax1=ixImax1-3;ixOmax2=ixImax2-3;ixOmax3=ixImax3-3;

    hxfimin1=ixOmin1+kr(idimm2,1);hxfimin2=ixOmin2+kr(idimm2,2)
    hxfimin3=ixOmin3+kr(idimm2,3);hxfimax1=ixOmax1+kr(idimm2,1)
    hxfimax2=ixOmax2+kr(idimm2,2);hxfimax3=ixOmax3+kr(idimm2,3);
    hxbimin1=ixOmin1-kr(idimm2,1);hxbimin2=ixOmin2-kr(idimm2,2)
    hxbimin3=ixOmin3-kr(idimm2,3);hxbimax1=ixOmax1-kr(idimm2,1)
    hxbimax2=ixOmax2-kr(idimm2,2);hxbimax3=ixOmax3-kr(idimm2,3);

    hxfifjmin1=hxfimin1+kr(idimm,1);hxfifjmin2=hxfimin2+kr(idimm,2)
    hxfifjmin3=hxfimin3+kr(idimm,3);hxfifjmax1=hxfimax1+kr(idimm,1)
    hxfifjmax2=hxfimax2+kr(idimm,2);hxfifjmax3=hxfimax3+kr(idimm,3);
    hxfibjmin1=hxfimin1-kr(idimm,1);hxfibjmin2=hxfimin2-kr(idimm,2)
    hxfibjmin3=hxfimin3-kr(idimm,3);hxfibjmax1=hxfimax1-kr(idimm,1)
    hxfibjmax2=hxfimax2-kr(idimm,2);hxfibjmax3=hxfimax3-kr(idimm,3);

    hxbifjmin1=hxbimin1+kr(idimm,1);hxbifjmin2=hxbimin2+kr(idimm,2)
    hxbifjmin3=hxbimin3+kr(idimm,3);hxbifjmax1=hxbimax1+kr(idimm,1)
    hxbifjmax2=hxbimax2+kr(idimm,2);hxbifjmax3=hxbimax3+kr(idimm,3);
    hxbibjmin1=hxbimin1-kr(idimm,1);hxbibjmin2=hxbimin2-kr(idimm,2)
    hxbibjmin3=hxbimin3-kr(idimm,3);hxbibjmax1=hxbimax1-kr(idimm,1)
    hxbibjmax2=hxbimax2-kr(idimm,2);hxbibjmax3=hxbimax3-kr(idimm,3);

    res(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
       1d0/(8d0*dxlevel(idimm) * dxlevel(idimm2))*((nu_hyper(&
       hxfifjmin1:hxfifjmax1,hxfifjmin2:hxfifjmax2,&
       hxfifjmin3:hxfifjmax3) + nu_hyper(hxfimin1:hxfimax1,hxfimin2:hxfimax2,&
       hxfimin3:hxfimax3)) * var2(hxfimin1:hxfimax1,hxfimin2:hxfimax2,&
       hxfimin3:hxfimax3) * (var(hxfifjmin1:hxfifjmax1,hxfifjmin2:hxfifjmax2,&
       hxfifjmin3:hxfifjmax3)-var(hxfibjmin1:hxfibjmax1,hxfibjmin2:hxfibjmax2,&
       hxfibjmin3:hxfibjmax3))-(nu_hyper(hxbifjmin1:hxbifjmax1,&
       hxbifjmin2:hxbifjmax2,hxbifjmin3:hxbifjmax3) + &
       nu_hyper(hxbimin1:hxbimax1,hxbimin2:hxbimax2,&
       hxbimin3:hxbimax3)) * var2(hxbimin1:hxbimax1,hxbimin2:hxbimax2,&
       hxbimin3:hxbimax3) * (var(hxbifjmin1:hxbifjmax1,hxbifjmin2:hxbifjmax2,&
       hxbifjmin3:hxbifjmax3)-var(hxbibjmin1:hxbibjmax1,hxbibjmin2:hxbibjmax2,&
       hxbibjmin3:hxbibjmax3)))

    !print*, "SECOND CROSS DERIV2 IXO ", ixO^L 

  end subroutine second_cross_deriv2

end module mod_hyperdiffusivity
