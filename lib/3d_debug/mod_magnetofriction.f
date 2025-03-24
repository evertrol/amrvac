!> module mod_magnetofriction.t
!> Purpose: use magnetofrictional method to relax 3D magnetic field to
!>          force-free field
!> 01.04.2016 developed by Chun Xia and Yang Guo
!> 04.10.2017 modulized by Chun Xia
!> Usage:
!> in amrvac.par:
!>   &methodlist
!>    time_stepper='onestep' ! time marching scheme, or 'twostep','threestep'
!>    flux_method=13*'cd4' ! or 'tvdlf', 'fd'
!>    limiter= 13*'koren' ! or 'vanleer','cada3','mp5' so on
!>   /
!>   &meshlist
!>    ditregrid=20 ! set iteration interval for adjusting AMR
!>   /
!>   &mhd_list
!>    mhd_magnetofriction=.true.
!>   /
!>   &mf_list
!>    mf_it_max=60000  ! set the maximum iteration number
!>    mf_ditsave=20000 ! set iteration interval for data output
!>    mf_cc=0.3     ! stability coefficient controls numerical stability
!>    mf_cy=0.2     ! frictional velocity coefficient
!>    mf_cdivb=0.01 ! divb cleaning coefficient controls diffusion speed of divb
!>   /
module mod_magnetofriction
  implicit none

  !> stability coefficient controls numerical stability
  double precision :: mf_cc
  !> frictional velocity coefficient
  double precision :: mf_cy, mf_cy_max
  !> divb cleaning coefficient controls diffusion speed of divb
  double precision :: mf_cdivb, mf_cdivb_max
  !> TVDLF dissipation coefficient controls the dissipation term
  double precision :: mf_tvdlfeps, mf_tvdlfeps_min
  !> time in magnetofriction process
  double precision :: tmf
  !> maximal speed for fd scheme
  double precision :: cmax_mype
  !> maximal speed for fd scheme
  double precision :: cmax_global
  !> maximal limit of magnetofrictional velocity in cm s^-1 (Pomoell 2019 A&A)
  double precision, public :: mf_vmax = 3.d6

  !> Index of the density (in the w array)
  integer, private, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, private, protected :: mom(:)

  !> Indices of the magnetic field
  integer, allocatable, private, protected :: mag(:)

  integer :: mf_ditsave
  integer :: mf_it_max
  integer :: mf_it
  logical :: mf_advance
  logical :: fix_conserve_at_step = .true.

contains
  !> Read this module"s parameters from a file
  subroutine mf_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mf_list/ mf_ditsave, mf_it_max, mf_it, mf_cc, mf_cy, mf_cy_max,&
        mf_cdivb, mf_cdivb_max, mf_tvdlfeps, mf_tvdlfeps_min

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mf_list, end=111)
111    close(unitpar)
    end do

  end subroutine mf_params_read

  !> Initialize the module
  subroutine magnetofriction_init()
    use mod_global_parameters
    use mod_usr_methods, only: usr_before_main_loop

    rho_=iw_rho
    allocate(mom(ndir))
    mom=iw_mom
    allocate(mag(ndir))
    mag=iw_mag

    mf_it=0          ! set the initial iteration number
    mf_it_max=60000  ! set the maximum iteration number
    mf_ditsave=20000 ! set iteration interval for data output
    mf_cc=0.3d0      ! stability coefficient controls numerical stability
    mf_cy=0.2d0 !frictional velocity coefficient. The default value is mf_cy itself.
    mf_cy_max=mf_cy  ! maximum of the frictional velocity coefficient
    mf_cdivb=0.01d0 !divb cleaning coefficient controls diffusion speed of divb
    mf_cdivb_max=mf_cdivb ! maximum of the divb cleaning coefficient
    mf_tvdlfeps=1.d0 ! coefficient to control the TVDLF dissipation
    mf_tvdlfeps_min = mf_tvdlfeps !minimum of the TVDLF dissipation coefficient
    ! get dimensionless maximal mf velocity limit
    mf_vmax=mf_vmax/unit_velocity

    call mf_params_read(par_files)

    if(mf_cy_max < mf_cy) mf_cy_max = mf_cy
    if(mf_cdivb_max < mf_cdivb) mf_cdivb_max = mf_cdivb
    if(mf_tvdlfeps_min > mf_tvdlfeps) mf_tvdlfeps_min = mf_tvdlfeps

    usr_before_main_loop => magnetofriction

  end subroutine magnetofriction_init

  subroutine magnetofriction
    use mod_global_parameters
    use mod_physics
    use mod_ghostcells_update
    use mod_input_output
    use mod_amr_grid, only: resettree
    use mod_comm_lib, only: mpistop

    double precision :: dvolume(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
       dsurface(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),dvone
    double precision :: dtfff,dtfff_pe,dtnew,dx1,dx2,dx3
    double precision :: cwsin_theta_new,cwsin_theta_old
    double precision :: sum_jbb,sum_jbb_ipe,sum_j,sum_j_ipe,sum_l_ipe,sum_l
    double precision :: f_i_ipe,f_i,volumepe,volume,tmpt,time_in
    double precision, external :: integral_grid
    integer :: i,iigrid, igrid, idims,ix1,ix2,ix3,hxMlo1,hxMlo2,hxMlo3,hxMhi1,&
       hxMhi2,hxMhi3,fhmf,tmpit,i1,i2,i3
    logical :: patchwi(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3),&
        stagger_flag

    ! not do fix conserve and getbc for staggered values if stagger is used
    stagger_flag=stagger_grid
    stagger_grid=.false.

    time_in=MPI_WTIME()
    if(mype==0) write(*,*)&
        'Evolving to force-free field using magnetofricitonal method...'
    if(prolongprimitive) call mpistop&
       ('use prolongprimitive=.false. in MF module')
    mf_advance=.false.
    dtfff=1.d-2
    tmpt=global_time
    tmpit=it
    tmf=global_time
    i=mf_it
    if(snapshotini==0 .and. i==0) then
      call saveamrfile(1)
      call saveamrfile(2)
    end if
    mf_advance=.true.
    ! point bc mpi datatype to partial type for magnetic field
    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1
    ! create bc mpi datatype for ghostcells update
    call create_bc_mpi_datatype(mag(1),ndir)
    ! point bc mpi datatype to partial type for velocity field
    type_send_srl=>type_send_srl_p2
    type_recv_srl=>type_recv_srl_p2
    type_send_r=>type_send_r_p2
    type_recv_r=>type_recv_r_p2
    type_send_p=>type_send_p_p2
    type_recv_p=>type_recv_p_p2
    ! create bc mpi datatype for ghostcells update
    call create_bc_mpi_datatype(mom(1),ndir)
    ! convert conservative variables to primitive ones which are used during MF
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call phys_to_primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
          ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x)
    end do
    ! calculate magnetofrictional velocity
    call mf_velocity_update(dtfff)
    ! update velocity in ghost cells
    call getbc(tmf,0.d0,ps,mom(1),ndir)
    ! calculate initial values of metrics
    if(i==0) then
      call metrics
      call printlog_mf
    end if
    ! magnetofrictional loops
    do
      ! calculate time step based on Cmax= Alfven speed + abs(frictional speed)
      dtfff_pe=bigdouble
      cmax_mype=zero
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)
        dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
        dxlevel(3)=rnode(rpdx3_,igrid);
        call getdtfff_courant(ps(igrid)%w,ps(igrid)%x,ixGlo1,ixGlo2,ixGlo3,&
           ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
           dtnew)
        dtfff_pe=min(dtfff_pe,dtnew)
      end do
      call MPI_ALLREDUCE(dtfff_pe,dtfff,1,MPI_DOUBLE_PRECISION,MPI_MIN, icomm,&
         ierrmpi)
      call MPI_ALLREDUCE(cmax_mype,cmax_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
         icomm,ierrmpi)

      ! =======
      ! evolve
      ! =======
      call advectmf(1,ndim,tmf,dtfff)

      if(i>=10000)  then
        mf_tvdlfeps = 0.9998d0 * mf_tvdlfeps
        mf_tvdlfeps = max(mf_tvdlfeps_min,mf_tvdlfeps)
      end if
      if(i<=60000) then
        mf_cy=1.0001d0*mf_cy
        mf_cy = min(mf_cy_max,mf_cy)
      end if
      if(i>=100000) then
        mf_cdivb=1.0000018d0*mf_cdivb
        mf_cdivb=min(mf_cdivb_max,mf_cdivb)
      end if

      i=i+1
      tmf=tmf+dtfff
      if(mod(i,10)==0) then
        ! calculate metrics
        call metrics
        call printlog_mf
      end if
      if(mod(i,mf_ditsave)==0) then
        it=i
        global_time=tmf
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call phys_to_conserved(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ps(igrid)%w,&
             ps(igrid)%x)
        end do
        mf_advance=.false.
        call saveamrfile(1)
        call saveamrfile(2)
        do iigrid=1,igridstail; igrid=igrids(iigrid);
           call phys_to_primitive(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
              ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ps(igrid)%w,&
              ps(igrid)%x)
        end do
        mf_advance=.true.
        if(mype==0) then
          write(*,*) "itmf=",i
          write(*,*) '<CW sin theta>:',cwsin_theta_new
          write(*,*) '<f_i>:',f_i
          write(*,*)&
              '----------------------------------------------------------'
        end if
      end if
      ! reconstruct AMR grid every 10 step
      if(mod(i,ditregrid)==0 .and. refine_max_level>1) call resettree
      if (i>=mf_it_max) then
        if(mod(i,10)/=0) then
          ! calculate metrics
          call metrics
          call printlog_mf
        end if
        if(mype==0) then
          write (*,*) 'Reach maximum iteration step!'
          write (*,*) 'The total iteration step is:', i
        end if
        exit
      end if
    enddo
    ! point bc mpi data type back to full type for MHD
    type_send_srl=>type_send_srl_f
    type_recv_srl=>type_recv_srl_f
    type_send_r=>type_send_r_f
    type_recv_r=>type_recv_r_f
    type_send_p=>type_send_p_f
    type_recv_p=>type_recv_p_f
    bcphys=.true.
    ! set velocity back to zero and convert primitive variables back to conservative ones
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,mom(:))=zero
       call phys_to_conserved(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
          ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x)
    end do
    global_time=tmpt
    it=tmpit
    if (mype==0) call MPI_FILE_CLOSE(fhmf,ierrmpi)
    mf_advance=.false.
    ! restore stagger_grid value
    stagger_grid=stagger_flag
    if(mype==0) write(*,*) 'Magnetofriction phase took : ',MPI_WTIME()-time_in,&
       ' sec'
    contains

      subroutine metrics

        sum_jbb_ipe = 0.d0
        sum_j_ipe = 0.d0
        sum_l_ipe = 0.d0
        f_i_ipe = 0.d0
        volumepe=0.d0
        dsurface=zero
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          block=>ps(igrid)
          dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             ixMlo3:ixMhi3)=block%dvolume(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             ixMlo3:ixMhi3)
          do idims=1,ndim
            hxMlo1=ixMlo1-kr(idims,1);hxMlo2=ixMlo2-kr(idims,2)
            hxMlo3=ixMlo3-kr(idims,3);hxMhi1=ixMhi1-kr(idims,1)
            hxMhi2=ixMhi2-kr(idims,2);hxMhi3=ixMhi3-kr(idims,3);
            dsurface(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
               ixMlo3:ixMhi3)=dsurface(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
               ixMlo3:ixMhi3)+block%surfaceC(hxMlo1:hxMhi1,hxMlo2:hxMhi2,&
               hxMlo3:hxMhi3,idims)
          end do
          dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
          dxlevel(3)=rnode(rpdx3_,igrid);
          call mask_inner(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixMlo1,&
             ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,ps(igrid)%x)
          sum_jbb_ipe = sum_jbb_ipe+integral_grid_mf(ixGlo1,ixGlo2,ixGlo3,&
             ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
             ps(igrid)%w,ps(igrid)%x,1,patchwi)
          sum_j_ipe   = sum_j_ipe+integral_grid_mf(ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
             ps(igrid)%w,ps(igrid)%x,2,patchwi)
          f_i_ipe=f_i_ipe+integral_grid_mf(ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
             ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,ps(igrid)%w,&
             ps(igrid)%x,3,patchwi)
          sum_l_ipe   = sum_l_ipe+integral_grid_mf(ixGlo1,ixGlo2,ixGlo3,ixGhi1,&
             ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
             ps(igrid)%w,ps(igrid)%x,4,patchwi)
        end do
        call MPI_ALLREDUCE(sum_jbb_ipe,sum_jbb,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           icomm,ierrmpi)
        call MPI_ALLREDUCE(sum_j_ipe,sum_j,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           icomm,ierrmpi)
        call MPI_ALLREDUCE(f_i_ipe,f_i,1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,&
           ierrmpi)
        call MPI_ALLREDUCE(sum_l_ipe,sum_l,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           icomm,ierrmpi)
        call MPI_ALLREDUCE(volumepe,volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
           icomm,ierrmpi)
        ! current- and volume-weighted average of the sine of the angle
        ! between the magnetic field B and the current density J
        cwsin_theta_new = sum_jbb/sum_j
        ! volume-weighted average of the absolute value of the fractional
        ! magnetic flux change
        f_i = f_i/volume
        sum_j=sum_j/volume
        sum_l=sum_l/volume
      end subroutine metrics

      subroutine mask_inner(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

        integer, intent(in)         :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
           ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
        double precision, intent(in):: w(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,nw),x(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndim)
        double precision            :: xOmin1,xOmin2,xOmin3,xOmax1,xOmax2,&
           xOmax3
        integer                     :: ix1,ix2,ix3

        xOmin1 = xprobmin1 + 0.05d0*(xprobmax1-xprobmin1)
        xOmin2 = xprobmin2 + 0.05d0*(xprobmax2-xprobmin2)
        xOmin3 = xprobmin3 + 0.05d0*(xprobmax3-xprobmin3)
        xOmax1 = xprobmax1 - 0.05d0*(xprobmax1-xprobmin1)
        xOmax2 = xprobmax2 - 0.05d0*(xprobmax2-xprobmin2)
        xOmax3 = xprobmax3 - 0.05d0*(xprobmax3-xprobmin3)
        if(slab) then
          xOmin3 = xprobmin3
        else
          xOmin1 = xprobmin1
        end if

        do ix3=ixOmin3,ixOmax3
        do ix2=ixOmin2,ixOmax2
        do ix1=ixOmin1,ixOmax1
            if( x(ix1,ix2,ix3,1) > xOmin1 .and. x(ix1,ix2,ix3,&
               1) < xOmax1  .and.  x(ix1,ix2,ix3,2) > xOmin2 .and. x(ix1,ix2,&
               ix3,2) < xOmax2  .and.  x(ix1,ix2,ix3,3) > xOmin3 .and. x(ix1,&
               ix2,ix3,3) < xOmax3 ) then
              patchwi(ix1,ix2,ix3)=.true.
              volumepe=volumepe+dvolume(ix1,ix2,ix3)
            else
              patchwi(ix1,ix2,ix3)=.false.
            endif
        end do
        end do
        end do

      end subroutine mask_inner

      subroutine printlog_mf
        integer :: amode, status(MPI_STATUS_SIZE)
        character(len=800) :: filename,filehead
        character(len=2048) :: line,datastr
        logical, save :: logmfopened=.false.

        if(mype==0) then
          if(.not.logmfopened) then
            ! generate filename
            write(filename,"(a,a)") TRIM(base_filename), "_mflog.csv"

            amode=ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)
            amode=ior(amode,MPI_MODE_APPEND)
            call MPI_FILE_OPEN(MPI_COMM_SELF,filename,amode,MPI_INFO_NULL,fhmf,&
               ierrmpi)
            logmfopened=.true.
filehead="  itmf,  dt,  <f_i>,  <CW sin theta>,  <Current>,  <Lorenz force>"
            call MPI_FILE_WRITE(fhmf,filehead,len_trim(filehead),&
                MPI_CHARACTER,status,ierrmpi)
            call MPI_FILE_WRITE(fhmf,achar(10),1,MPI_CHARACTER,status,ierrmpi)
          end if
          line=''
          write(datastr,'(i6,a)') i,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6,a)') dtfff,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6,a)') f_i,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6,a)') cwsin_theta_new,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6,a)') sum_j,','
          line=trim(line)//trim(datastr)
          write(datastr,'(es13.6)') sum_l
          line=trim(line)//trim(datastr)//new_line('A')
          call MPI_FILE_WRITE(fhmf,line,len_trim(line),MPI_CHARACTER,status,&
             ierrmpi)
        end if

      end subroutine printlog_mf

      function integral_grid_mf(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
         ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,iw,&
         patchwi)
        use mod_geometry

        integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
           ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iw
        double precision, intent(in)       :: x(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
        double precision, intent(in)       :: w(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3,nw+nwauxio)
        logical, intent(in) :: patchwi(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3)

        double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
           ixImin3:ixImax3,1:ndir) :: bvec,qvec,current
        double precision :: integral_grid_mf,tmp(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3),b_mag(ixImin1:ixImax1,&
           ixImin2:ixImax2,ixImin3:ixImax3)
        integer :: ix1,ix2,ix3,i,idirmin,idir,jdir,kdir

        integral_grid_mf=0.d0
        select case(iw)
         case(1)
          ! Sum(dvolume*|JxB|/|B|)
          if(B0field) then
            bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               :)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               mag(:))+block%b0(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3,mag(:),0)
          else
            bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               :)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(:))
          endif
          call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
          ! calculate Lorentz force
          qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)=zero
          do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
             if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,jdir)*bvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,kdir)
                if(lvc(idir,jdir,kdir)==1)then
                   qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                      idir)=qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3,idir)+tmp(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                else
                   qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                      idir)=qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3,idir)-tmp(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                endif
             endif
          enddo; enddo; enddo

          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
             if(patchwi(ix1,ix2,ix3)) integral_grid_mf=integral_grid_mf+&
                sqrt(sum(qvec(ix1,ix2,ix3,:)**2)/sum(bvec(ix1,ix2,ix3,&
                :)**2))*dvolume(ix1,ix2,ix3)
          end do
          end do
          end do
         case(2)
          ! Sum(dvolume*|J|)
          call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
             if(patchwi(ix1,ix2,ix3)) integral_grid_mf=integral_grid_mf+&
                sqrt(sum(current(ix1,ix2,ix3,:)**2))*dvolume(ix1,ix2,ix3)
          end do
          end do
          end do
         case(3)
          ! f_i solenoidal property of B: (dvolume |div B|)/(dsurface |B|)
          ! Sum(dvolume*f_i)
          if(B0field) then
            bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               :)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               mag(:))+block%b0(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3,mag(:),0)
          else
            bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               :)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(:))
          endif
          call divvector(bvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,tmp)
          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
             if(patchwi(ix1,ix2,ix3)) integral_grid_mf=integral_grid_mf+&
                abs(tmp(ix1,ix2,ix3))*dvolume(ix1,ix2,&
                ix3)**2/sqrt(sum(bvec(ix1,ix2,ix3,:)**2))/dsurface(ix1,ix2,&
                ix3)
          end do
          end do
          end do
         case(4)
          ! Sum(|JxB|)
          if(B0field) then
            bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               :)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               mag(:))+block%b0(ixImin1:ixImax1,ixImin2:ixImax2,&
               ixImin3:ixImax3,mag(:),0)
          else
            bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
               :)=w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mag(:))
          endif
          call curlvector(bvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,current,idirmin,1,&
             ndir)
          ! calculate Lorentz force
          qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)=zero
          do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
             if(lvc(idir,jdir,kdir)/=0)then
                tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,jdir)*bvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   ixOmin3:ixOmax3,kdir)
                if(lvc(idir,jdir,kdir)==1)then
                   qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                      idir)=qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3,idir)+tmp(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                else
                   qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                      idir)=qvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3,idir)-tmp(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3)
                endif
             endif
          enddo; enddo; enddo

          do ix3=ixOmin3,ixOmax3
          do ix2=ixOmin2,ixOmax2
          do ix1=ixOmin1,ixOmax1
             if(patchwi(ix1,ix2,ix3)) integral_grid_mf=integral_grid_mf+&
                sqrt(sum(qvec(ix1,ix2,ix3,:)**2))*dvolume(ix1,ix2,ix3)
          end do
          end do
          end do
        end select
        return
      end function integral_grid_mf

  end subroutine magnetofriction

  subroutine mf_velocity_update(dtfff)
    use mod_global_parameters

    double precision, intent(in) :: dtfff
    integer :: i,iigrid, igrid
    double precision :: vhatmax,vhatmax_pe,vhatmaxgrid

    vhatmax_pe=smalldouble
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
      dxlevel(3)=rnode(rpdx3_,igrid);
      call vhat(ps(igrid)%w,ps(igrid)%x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,&
         ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,vhatmaxgrid)
      vhatmax_pe=max(vhatmax_pe,vhatmaxgrid)
    end do
    call MPI_ALLREDUCE(vhatmax_pe,vhatmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
        icomm,ierrmpi)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
      dxlevel(3)=rnode(rpdx3_,igrid);
      ! calculate frictional velocity
      call frictional_velocity(ps(igrid)%w,ps(igrid)%x,ixGlo1,ixGlo2,ixGlo3,&
         ixGhi1,ixGhi2,ixGhi3,ixMlo1,ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,&
         vhatmax,dtfff)
    end do

  end subroutine mf_velocity_update

  subroutine vhat(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,vhatmaxgrid)
    ! Calculate v_hat
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out) :: vhatmaxgrid

    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,7-2*ndir:3),tmp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),dxhm
    double precision :: dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: idirmin,idir,jdir,kdir

    call get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mom(:))=0.d0
    ! calculate Lorentz force
    do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
       if(lvc(idir,jdir,kdir)/=0)then
          if(B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,jdir)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,mag(kdir))+block%b0(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,ixOmin3:ixOmax3,kdir,0))
          else
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,jdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               ixOmin3:ixOmax3,mag(kdir))
          endif
          if(lvc(idir,jdir,kdir)==1)then
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom(idir))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          else
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mom(idir))-tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          endif
       endif
    enddo; enddo; enddo

    ! 1/B**2
    if(B0field) then
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/(sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))+block%b0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:,0))**2,dim=ndim+1)+smalldouble)
    else
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=1.d0/(sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)+smalldouble)
    endif

    if(slab_uniform) then
      dxhm=dble(ndim)/(1.0d0/dxlevel(1)+1.0d0/dxlevel(2)+1.0d0/dxlevel(3))
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=dxhm*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do
    else
      dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=dble(ndim)/sum(1.d0/block%dx(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,:),dim=ndim+1)
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do
    end if
    vhatmaxgrid=maxval(sqrt(sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1)))

  end subroutine vhat

  subroutine frictional_velocity(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qvmax,qdt)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim),qdt,qvmax
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: dxhm,disbd(6),bfzone1,bfzone2,bfzone3
    double precision :: dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: ix1,ix2,ix3, idir
    logical :: buffer

    if(slab_uniform) then
      dxhm=dble(ndim)/(1.0d0/dxlevel(1)+1.0d0/dxlevel(2)+1.0d0/dxlevel(3))
      dxhm=mf_cc*mf_cy/qvmax*dxhm/qdt
      ! dxhm=mf_cc*mf_cy/qvmax
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(:))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(:))*dxhm
    else
      dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=dble(ndim)/sum(1.d0/block%dx(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,:),dim=ndim+1)
      dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=mf_cc*mf_cy/qvmax*dxhms(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)/qdt
      ! dxhm=mf_cc*mf_cy/qvmax
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mom(idir))*dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      end do
    end if


    bfzone1=0.05d0*(xprobmax1-xprobmin1)
    bfzone2=0.05d0*(xprobmax2-xprobmin2)
    bfzone3=0.05d0*(xprobmax3-xprobmin3)
    do ix3=ixOmin3,ixOmax3
do ix2=ixOmin2,ixOmax2
do ix1=ixOmin1,ixOmax1
       disbd(1)=x(ix1,ix2,ix3,1)-xprobmin1
       disbd(2)=xprobmax1-x(ix1,ix2,ix3,1)
       disbd(3)=x(ix1,ix2,ix3,2)-xprobmin2
       disbd(4)=xprobmax2-x(ix1,ix2,ix3,2)
       disbd(5)=x(ix1,ix2,ix3,3)-xprobmin1
       disbd(6)=xprobmax3-x(ix1,ix2,ix3,3)

       if(slab) then
         if(disbd(1)<bfzone1) then
           w(ix1,ix2,ix3,mom(:))=(1.d0-((bfzone1-disbd(1))/bfzone1)**2)*w(ix1,&
              ix2,ix3,mom(:))
         endif
       else
         if(disbd(5)<bfzone3) then
           w(ix1,ix2,ix3,mom(:))=(1.d0-((bfzone3-disbd(5))/bfzone3)**2)*w(ix1,&
              ix2,ix3,mom(:))
         endif
       end if
       if(disbd(2)<bfzone1) then
         w(ix1,ix2,ix3,mom(:))=(1.d0-((bfzone1-disbd(2))/bfzone1)**2)*w(ix1,&
            ix2,ix3,mom(:))
       endif
       if(disbd(3)<bfzone2) then
         w(ix1,ix2,ix3,mom(:))=(1.d0-((bfzone2-disbd(3))/bfzone2)**2)*w(ix1,&
            ix2,ix3,mom(:))
       endif
       if(disbd(4)<bfzone2) then
         w(ix1,ix2,ix3,mom(:))=(1.d0-((bfzone2-disbd(4))/bfzone2)**2)*w(ix1,&
            ix2,ix3,mom(:))
       endif
       if(disbd(6)<bfzone3) then
         w(ix1,ix2,ix3,mom(:))=(1.d0-((bfzone3-disbd(6))/bfzone3)**2)*w(ix1,&
            ix2,ix3,mom(:))
       endif
    end do
end do
end do


    ! saturate mf velocity at mf_vmax
    dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=sqrt(sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mom(:))**2,dim=ndim+1))/mf_vmax+1.d-12
    dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)=dtanh(dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3))/dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(idir))*dxhms(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    end do
  end subroutine frictional_velocity

  subroutine advectmf(idimmin,idimmax,qt,qdt)
    !  integrate all grids by one step of its delta(global_time)
    ! This subroutine is in VAC terminology equivalent to
    ! `advect' (with the difference that it will `advect' all grids)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: idimmin,idimmax
    double precision, intent(in) :: qt, qdt

    integer :: iigrid, igrid

    call init_comm_fix_conserve(idimmin,idimmax,ndir)
    fix_conserve_at_step = mf_advance .and. levmax>levmin

    do iigrid=1,igridstail; igrid=igrids(iigrid);
      ps1(igrid)%w=ps(igrid)%w
    end do

    istep=0

    select case (t_stepper)
     case (onestep)
       call advect1mf(flux_method,qdt,one,idimmin,idimmax,qt,ps1,qt,ps)
     case (twostep)
       ! predictor step
       fix_conserve_at_step = .false.
       call advect1mf(typepred1,qdt,half,idimmin,idimmax,qt,ps,qt,ps1)
       ! corrector step
       fix_conserve_at_step = mf_advance .and. levmax>levmin
       call advect1mf(flux_method,qdt,one,idimmin,idimmax,qt+half*qdt,ps1,qt,&
          ps)
     case (threestep)
       ! three step Runge-Kutta in accordance with Gottlieb & Shu 1998
       call advect1mf(flux_method,qdt,one,idimmin,idimmax,qt,ps,qt,ps1)

       do iigrid=1,igridstail; igrid=igrids(iigrid);
          ps2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
             1:nwflux)=0.75d0*ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             ixGlo3:ixGhi3,1:nwflux)+0.25d0*ps1(igrid)%w(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
          if (nw>nwflux) ps2(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             ixGlo3:ixGhi3,nwflux+1:nw) = ps(igrid)%w(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,ixGlo3:ixGhi3,nwflux+1:nw)
       end do

       call advect1mf(flux_method,qdt,0.25d0,idimmin,idimmax,qt+qdt,ps1,&
          qt+dt*0.25d0,ps2)

       do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
          ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
             1:nwflux)=1.0d0/3.0d0*ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
             ixGlo3:ixGhi3,1:nwflux)+2.0d0/3.0d0*ps2(igrid)%w(ixGlo1:ixGhi1,&
             ixGlo2:ixGhi2,ixGlo3:ixGhi3,1:nwflux)
       end do

       call advect1mf(flux_method,qdt,2.0d0/3.0d0,idimmin,idimmax,qt+qdt/2.0d0,&
          ps2,qt+qdt/3.0d0,ps)
     case default
       call mpistop("unkown time_stepper in advectmf")
    end select

  end subroutine advectmf

  subroutine advect1mf(method,dtin,dtfactor,idimmin,idimmax,qtC,psa,qt,psb)
    ! Integrate all grids by one partial step
    ! This subroutine is equivalent to VAC's `advect1', but does
    ! the advection for all grids
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve

    integer, intent(in) :: idimmin,idimmax
    type(state) :: psa(max_blocks)! Compute fluxes based on this state
    type(state) :: psb(max_blocks) ! update on this state
    double precision, intent(in) :: dtin,dtfactor, qtC, qt
    integer, intent(in) :: method(nlevelshi)

    double precision :: qdt
    integer :: iigrid, igrid, level, i1,i2,i3
    logical :: setigrid

    istep=istep+1

    ! loop over all grids to arrive at equivalent
    qdt=dtfactor*dtin
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       block=>ps(igrid)
       level=node(plevel_,igrid)

       call process1_gridmf(method(level),igrid,qdt,ixGlo1,ixGlo2,ixGlo3,&
          ixGhi1,ixGhi2,ixGhi3,idimmin,idimmax,qtC,psa(igrid)%w,qt,&
          psb(igrid)%w)
    end do

    ! opedit: Send flux for all grids, expects sends for all
    ! nsend_fc(^D), set in connectivity.t.

    if (fix_conserve_at_step) then
      call recvflux(idimmin,idimmax)
      call sendflux(idimmin,idimmax)
      call fix_conserve(psb,idimmin,idimmax,mag(1),ndir)
    end if
    ! point bc mpi datatype to partial type for magnetic field
    type_send_srl=>type_send_srl_p1
    type_recv_srl=>type_recv_srl_p1
    type_send_r=>type_send_r_p1
    type_recv_r=>type_recv_r_p1
    type_send_p=>type_send_p_p1
    type_recv_p=>type_recv_p_p1
    ! update B in ghost cells
    call getbc(qt+qdt,qdt,psb,mag(1),ndir)
    ! calculate magnetofrictional velocity
    call mf_velocity_update(qdt)
    ! point bc mpi datatype to partial type for velocity field
    type_send_srl=>type_send_srl_p2
    type_recv_srl=>type_recv_srl_p2
    type_send_r=>type_send_r_p2
    type_recv_r=>type_recv_r_p2
    type_send_p=>type_send_p_p2
    type_recv_p=>type_recv_p_p2
    ! update magnetofrictional velocity in ghost cells
    call getbc(qt+qdt,qdt,psb,mom(1),ndir)

  end subroutine advect1mf

  subroutine process1_gridmf(method,igrid,qdt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
     ixGmax2,ixGmax3,idimmin,idimmax,qtC,wCT,qt,w)
    ! This subroutine is equivalent to VAC's `advect1' for one grid
    use mod_global_parameters
    use mod_fix_conserve
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: method
    integer, intent(in) :: igrid, ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
       ixGmax3, idimmin,idimmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision :: wCT(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       1:nw), w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:nw)
    double precision :: dx1,dx2,dx3, fC(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndir,1:ndim)
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3

    dx1=rnode(rpdx1_,igrid);dx2=rnode(rpdx2_,igrid);dx3=rnode(rpdx3_,igrid);
    dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid)
    dxlevel(3)=rnode(rpdx3_,igrid);
    fC=0.d0

    ixOmin1=ixGmin1+nghostcells;ixOmin2=ixGmin2+nghostcells
    ixOmin3=ixGmin3+nghostcells;ixOmax1=ixGmax1-nghostcells
    ixOmax2=ixGmax2-nghostcells;ixOmax3=ixGmax3-nghostcells;
    select case (method)
     case (fs_cd4)
       !================================
       ! 4th order central difference
       !================================
       call centdiff4mf(qdt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,&
          wCT,qt,w,fC,dx1,dx2,dx3,ps(igrid)%x)
     case (fs_tvdlf)
       !================================
       ! TVDLF
       !================================
       call tvdlfmf(qdt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,&
          wCT,qt,w,fC,dx1,dx2,dx3,ps(igrid)%x)
     case (fs_hancock)
       ! hancock predict (first) step for twostep tvdlf and tvdmu scheme
       call hancockmf(qdt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,&
          wCT,qt,w,dx1,dx2,dx3,ps(igrid)%x)
     case (fs_fd)
       !================================
       ! finite difference
       !================================
       call fdmf(qdt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixOmin1,&
          ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,wCT,qt,w,&
          fC,dx1,dx2,dx3,ps(igrid)%x)
    case default
       call mpistop("unknown flux scheme in advect1_gridmf")
    end select

    if (fix_conserve_at_step) then
      call store_flux(igrid,fC,idimmin,idimmax,ndir)
    end if

  end subroutine process1_gridmf

  subroutine upwindLRmf(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,ixRmin1,ixRmin2,ixRmin3,&
     ixRmax1,ixRmax2,ixRmax3,idim,w,wCT,wLC,wRC,x)
    ! Determine the upwinded wLC(ixL) and wRC(ixR) from w.
    ! the wCT is only used when PPM is exploited.
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixLmin1,ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3, ixRmin1,ixRmin2,&
       ixRmin3,ixRmax1,ixRmax2,ixRmax3, idim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: w, wCT
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wLC, wRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim) :: x

    integer :: jxRmin1,jxRmin2,jxRmin3,jxRmax1,jxRmax2,jxRmax3, ixCmin1,&
       ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, jxCmin1,jxCmin2,jxCmin3,&
       jxCmax1,jxCmax2,jxCmax3, iw
    double precision   :: ldw(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        rdw(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        dwC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    if (type_limiter(block%level) == limiter_mp5) then
       call MP5limiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixLmin1,&
          ixLmin2,ixLmin3,ixLmax1,ixLmax2,ixLmax3,idim,w,wLC,wRC)
    else if (type_limiter(block%level) == limiter_ppm) then
       call PPMlimiter(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixMlo1,&
          ixMlo2,ixMlo3,ixMhi1,ixMhi2,ixMhi3,idim,w,wCT,wLC,wRC)
    else
       jxRmin1=ixRmin1+kr(idim,1);jxRmin2=ixRmin2+kr(idim,2)
       jxRmin3=ixRmin3+kr(idim,3);jxRmax1=ixRmax1+kr(idim,1)
       jxRmax2=ixRmax2+kr(idim,2);jxRmax3=ixRmax3+kr(idim,3);
       ixCmax1=jxRmax1;ixCmax2=jxRmax2;ixCmax3=jxRmax3
       ixCmin1=ixLmin1-kr(idim,1);ixCmin2=ixLmin2-kr(idim,2)
       ixCmin3=ixLmin3-kr(idim,3);
       jxCmin1=ixCmin1+kr(idim,1);jxCmin2=ixCmin2+kr(idim,2)
       jxCmin3=ixCmin3+kr(idim,3);jxCmax1=ixCmax1+kr(idim,1)
       jxCmax2=ixCmax2+kr(idim,2);jxCmax3=ixCmax3+kr(idim,3);

       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,ixCmin3:jxCmax3,&
                iw)=dlog10(w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,ixCmin3:jxCmax3,&
                iw))
             wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
                iw)=dlog10(wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
                iw))
             wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                iw)=dlog10(wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                iw))
          end if

          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             jxCmin3:jxCmax3,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw)

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idim,&
             type_limiter(block%level),ldw,rdw)
          wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
             iw)=wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
             iw)+half*ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3)
          wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
             iw)=wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
             iw)-half*rdw(jxRmin1:jxRmax1,jxRmin2:jxRmax2,jxRmin3:jxRmax3)

          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,ixCmin3:jxCmax3,&
                iw)=10.0d0**w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,ixCmin3:jxCmax3,&
                iw)
             wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,ixLmin3:ixLmax3,&
                iw)=10.0d0**wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                ixLmin3:ixLmax3,iw)
             wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,ixRmin3:ixRmax3,&
                iw)=10.0d0**wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                ixRmin3:ixRmax3,iw)
          end if
       end do

    endif

  end subroutine upwindLRmf

  subroutine getfluxmf(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,idim,f)
    ! Calculate lux f_idim[idir] within ixO^L.
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idir, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision,intent(out)    :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    if (idim==idir) then
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=zero
    else
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idim))*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idir))-w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(idim))*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(idir))
      if (B0field) then
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=f(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(idim))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idir,idim)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(idir))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,idim,idim)
      end if
    end if

  end subroutine getfluxmf

  subroutine tvdlfmf(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,wCT,&
     qt,wnew,fC,dx1,dx2,dx3,x)
    ! method=='tvdlf'  --> 2nd order TVD-Lax-Friedrich scheme.
    ! method=='tvdlf1' --> 1st order TVD-Lax-Friedrich scheme.
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    double precision, intent(in)                         :: qdt, qtC, qt, dx1,&
       dx2,dx3
    integer, intent(in)                                  :: ixImin1,ixImin2,&
       ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,&
       ixOmax2,ixOmax3, idimmin,idimmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), intent(in) ::  x
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)               :: wCT, wnew
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir,1:ndim)        :: fC

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wLC, wRC, wmean
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)      :: fLC, fRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)      :: cmaxC
    double precision :: dxinv(1:ndim), inv_volume(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: idims, idir, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ixCmin1,ixCmin2,ixCmin3,&
       ixCmax1,ixCmax2,ixCmax3, ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,&
       ixCRmax3, jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3, kxCmin1,&
       kxCmin2,kxCmin3,kxCmax1,kxCmax2,kxCmax3, kxRmin1,kxRmin2,kxRmin3,&
       kxRmax1,kxRmax2,kxRmax3

    ! The flux calculation contracts by one in the idim direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmin3=ixOmin3;ixmax1=ixOmax1
    ixmax2=ixOmax2;ixmax3=ixOmax3;
    do idims= idimmin,idimmax
       ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
       ixmin3=ixmin3-2*kr(idims,3);ixmax1=ixmax1+2*kr(idims,1)
       ixmax2=ixmax2+2*kr(idims,2);ixmax3=ixmax3+2*kr(idims,3);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImin3>ixmin3.or.ixImax1<ixmax1.or.&
       ixImax2<ixmax2.or.ixImax3<ixmax3) call &
       mpistop("Error in tvdlfmf: Nonconforming input limits")

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;dxinv(3)=-qdt/dx3;
    fC=0.d0
    do idims= idimmin,idimmax
       b0i=idims

       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
       hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
       ! ixC is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3; ixCmin1=hxOmin1
       ixCmin2=hxOmin2;ixCmin3=hxOmin3;
       ! Calculate wRC=uR_{j+1/2} and wLC=uL_j+1/2
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
       jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);
       kxCmin1=ixImin1;kxCmin2=ixImin2;kxCmin3=ixImin3
       kxCmax1=ixImax1-kr(idims,1);kxCmax2=ixImax2-kr(idims,2)
       kxCmax3=ixImax3-kr(idims,3);
       kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
       kxRmin3=kxCmin3+kr(idims,3);kxRmax1=kxCmax1+kr(idims,1)
       kxRmax2=kxCmax2+kr(idims,2);kxRmax3=kxCmax3+kr(idims,3);
       ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmin3=ixCmin3;ixCRmax1=ixCmax1
       ixCRmax2=ixCmax2;ixCRmax3=ixCmax3;

       wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,&
          1:nwflux)=wCT(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,&
          1:nwflux)
       wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,&
          1:nwflux)=wCT(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,&
          1:nwflux)

       call upwindLRmf(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,ixCRmin1,&
          ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,idims,wCT,wCT,wLC,wRC,&
          x)

       ! For the high order Lax-Friedrich TVDLF scheme the limiter is based on
       ! the maximum eigenvalue, it is calculated in advance.
       ! determine mean state and store in wLC
       wmean=0.5d0*(wLC+wRC)
       call getcmaxfff(wmean,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixCmin1,&
          ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxC)

       ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each idir
       do idir=1,ndir
          call getfluxmf(wLC,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir,idims,fLC)
          call getfluxmf(wRC,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idir,idims,fRC)
          ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
          fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=half*(fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3))

          ! Add TVDLF dissipation to the flux
          if (idir==idims) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)-mf_tvdlfeps*tvdlfeps*cmaxC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3)*half*(wRC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(idir))-wLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(idir)))
          end if
          if (slab_uniform) then
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir,&
               idims)=fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3)
          else
            fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir,&
               idims)=block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3,idims)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               ixCmin3:ixCmax3)
          end if

       end do ! Next idir
    end do ! Next idims
    b0i=0

    !Now update the state:
    do idims= idimmin,idimmax
       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
       hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab_uniform) then
          fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:,&
             idims)=dxinv(idims)*fC(ixImin1:ixImax1,ixImin2:ixImax2,&
             ixImin3:ixImax3,:,idims)
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(:))=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
             mag(:)) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,:,&
             idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,:,&
             idims))
       else
          inv_volume = 1.0d0/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)
          fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:,&
             idims)=-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:,&
             idims)

          do idir = 1, ndir
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mag(idir))=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,mag(idir)) + (fC(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,idims)-fC(hxOmin1:hxOmax1,&
                hxOmin2:hxOmax2,hxOmin3:hxOmax3,idir,idims)) * inv_volume
          end do
       end if

    end do ! Next idims

    if (.not.slab) call addgeometrymf(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
       wnew,x)
    call divbclean(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,wnew,x)

  end subroutine tvdlfmf

  subroutine hancockmf(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,wCT,&
     qt,wnew,dx1,dx2,dx3,x)
    ! The non-conservative Hancock predictor for TVDLFmf
    ! on entry:
    ! input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
    ! one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
    ! on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idimmin,idimmax
    double precision, intent(in) :: qdt, qtC, qt, dx1,dx2,dx3,&
        x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), wnew(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wLC, wRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) :: fLC, fRC
    double precision :: dxinv(1:ndim)
    integer :: idims, idir, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxOmin1,&
       hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3, ixtestmin1,ixtestmin2,&
       ixtestmin3,ixtestmax1,ixtestmax2,ixtestmax3

    ! Expand limits in each idims direction in which fluxes are added
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmin3=ixOmin3;ixmax1=ixOmax1
    ixmax2=ixOmax2;ixmax3=ixOmax3;
    do idims= idimmin,idimmax
       ixmin1=ixmin1-kr(idims,1);ixmin2=ixmin2-kr(idims,2)
       ixmin3=ixmin3-kr(idims,3);ixmax1=ixmax1+kr(idims,1)
       ixmax2=ixmax2+kr(idims,2);ixmax3=ixmax3+kr(idims,3);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImin3>ixmin3.or.ixImax1<ixmax1.or.&
       ixImax2<ixmax2.or.ixImax3<ixmax3) call &
       mpistop("Error in Hancockmf: Nonconforming input limits")

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;dxinv(3)=-qdt/dx3;
    do idims= idimmin,idimmax
       b0i=idims
       ! Calculate w_j+g_j/2 and w_j-g_j/2
       ! First copy all variables, then upwind wLC and wRC.
       ! wLC is to the left of ixO, wRC is to the right of wCT.
       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
       hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);

       wRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
          1:nwflux)=wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1:nwflux)
       wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1:nwflux)=wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
          1:nwflux)

       call upwindLRmf(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
          ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,hxOmin1,hxOmin2,hxOmin3,&
          hxOmax1,hxOmax2,hxOmax3,idims,wCT,wCT,wLC,wRC,x)

       ! Advect mag(idir)
       do idir=1,ndir
          ! Calculate the fLC and fRC fluxes
          call getfluxmf(wRC,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3,idir,idims,fRC)
          call getfluxmf(wLC,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idir,idims,fLC)

          if (slab_uniform) then
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mag(idir))=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,mag(idir))+dxinv(idims)* (fLC(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)-fRC(hxOmin1:hxOmax1,&
                hxOmin2:hxOmax2,hxOmin3:hxOmax3))
          else
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mag(idir))=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,mag(idir))-qdt/block%dvolume(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3) &
                *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,idims)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3) -block%surfaceC(hxOmin1:hxOmax1,&
                hxOmin2:hxOmax2,hxOmin3:hxOmax3,idims)*fRC(hxOmin1:hxOmax1,&
                hxOmin2:hxOmax2,hxOmin3:hxOmax3))
          end if
       end do
    end do ! next idims
    b0i=0

    if (.not.slab) call addgeometrymf(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
       wnew,x)

  end subroutine hancockmf

  subroutine fdmf(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,wCT,qt,wnew,&
     fC,dx1,dx2,dx3,x)
    use mod_global_parameters
    double precision, intent(in)                                     :: qdt,&
        qtC, qt, dx1,dx2,dx3
    integer, intent(in)                                              :: &
       ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixOmin1,ixOmin2,&
       ixOmin3,ixOmax1,ixOmax2,ixOmax3, idimmin,idimmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim), intent(in)            :: x

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), intent(inout)           :: wCT, wnew
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir,1:ndim), intent(out)  :: fC

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                               :: fCT
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)                          :: fm, fp, fmR, fpL
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)                               :: v
    double precision                                                 :: &
       dxinv(1:ndim)
    integer                                                          :: idims,&
        idir, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, ixmin1,ixmin2,&
       ixmin3,ixmax1,ixmax2,ixmax3, hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,&
       hxOmax3, ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;dxinv(3)=-qdt/dx3;
    do idims= idimmin,idimmax

       ! Get fluxes for the whole grid (mesh+nghostcells)
        ixCmin1 = ixOmin1 - nghostcells * kr(idims,1)
         ixCmin2 = ixOmin2 - nghostcells * kr(idims,2)
         ixCmin3 = ixOmin3 - nghostcells * kr(idims,3)
        ixCmax1 = ixOmax1 + nghostcells * kr(idims,1)
         ixCmax2 = ixOmax2 + nghostcells * kr(idims,2)
         ixCmax3 = ixOmax3 + nghostcells * kr(idims,3)

       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
       hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);
       ! ix is centered index in the idim direction from ixOmin-1/2 to ixOmax+1/2
       ixmax1=ixOmax1;ixmax2=ixOmax2;ixmax3=ixOmax3; ixmin1=hxOmin1
       ixmin2=hxOmin2;ixmin3=hxOmin3;
       ixCRmin1=ixCmin1;ixCRmin2=ixCmin2;ixCRmin3=ixCmin3;ixCRmax1=ixCmax1
       ixCRmax2=ixCmax2;ixCRmax3=ixCmax3;

       do idir=1,ndir
          call getfluxmf(wCT,x,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,&
             ixCRmin1,ixCRmin2,ixCRmin3,ixCRmax1,ixCRmax2,ixCRmax3,idir,idims,&
             fCT)
          ! Lax-Friedrich splitting:
          fp(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
             mag(idir)) = half * (fCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,&
             ixCRmin3:ixCRmax3) + mf_tvdlfeps * tvdlfeps * cmax_global * &
             wCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
             mag(idir)))
          fm(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
             mag(idir)) = half * (fCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,&
             ixCRmin3:ixCRmax3) - mf_tvdlfeps * tvdlfeps * cmax_global * &
             wCT(ixCRmin1:ixCRmax1,ixCRmin2:ixCRmax2,ixCRmin3:ixCRmax3,&
             mag(idir)))
       end do ! iw loop

       ! now do the reconstruction of fp and fm:
       call reconstructLmf(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idims,fp,fpL)
       call reconstructRmf(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idims,fm,fmR)

       do idir=1,ndir
          if (slab_uniform) then
             fC(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idir,&
                idims) = dxinv(idims) * (fpL(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3,mag(idir)) + fmR(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3,mag(idir)))
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mag(idir))=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,mag(idir))+ (fC(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,idims)-fC(hxOmin1:hxOmax1,&
                hxOmin2:hxOmax2,hxOmin3:hxOmax3,idir,idims))
          else
             fC(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,idir,&
                idims)=-qdt*block%surfaceC(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3,idims) * (fpL(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3,mag(idir)) + fmR(ixmin1:ixmax1,ixmin2:ixmax2,&
                ixmin3:ixmax3,mag(idir)))
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mag(idir))=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,mag(idir))+ (fC(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3,idir,idims)-fC(hxOmin1:hxOmax1,&
                hxOmin2:hxOmax2,hxOmin3:hxOmax3,idir,&
                idims))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3)
          end if
       end do ! iw loop

    end do !idims loop

    if (.not.slab) call addgeometrymf(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,&
       wnew,x)
    call divbclean(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,wnew,x)

  end subroutine fdmf

  subroutine reconstructLmf(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC)
    use mod_global_parameters
    use mod_limiter

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision, intent(out)   :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision                :: ldw(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), dwC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer                         :: jxRmin1,jxRmin2,jxRmin3,jxRmax1,jxRmax2,&
       jxRmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, jxCmin1,&
       jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3, kxCmin1,kxCmin2,kxCmin3,&
       kxCmax1,kxCmax2,kxCmax3, iw

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC)
    case (limiter_weno5)
       call WENO5limiterL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC,1)
    case (limiter_wenoz5)
       call WENO5limiterL(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wLC,2)
    case default

       kxCmin1=ixImin1;kxCmin2=ixImin2;kxCmin3=ixImin3
       kxCmax1=ixImax1-kr(idims,1);kxCmax2=ixImax2-kr(idims,2)
       kxCmax3=ixImax3-kr(idims,3);

       wLC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,&
          1:nwflux) = w(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,&
          1:nwflux)

       jxRmin1=iLmin1+kr(idims,1);jxRmin2=iLmin2+kr(idims,2)
       jxRmin3=iLmin3+kr(idims,3);jxRmax1=iLmax1+kr(idims,1)
       jxRmax2=iLmax2+kr(idims,2);jxRmax3=iLmax3+kr(idims,3);

       ixCmax1=jxRmax1;ixCmax2=jxRmax2;ixCmax3=jxRmax3
       ixCmin1=iLmin1-kr(idims,1);ixCmin2=iLmin2-kr(idims,2)
       ixCmin3=iLmin3-kr(idims,3);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
       jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);

       do iw=1,nwflux
          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             jxCmin3:jxCmax3,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw)

          call dwlimiter2(dwC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,&
             type_limiter(block%level),ldw=ldw)

          wLC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw)=wLC(iLmin1:iLmax1,&
             iLmin2:iLmax2,iLmin3:iLmax3,iw)+half*ldw(iLmin1:iLmax1,&
             iLmin2:iLmax2,iLmin3:iLmax3)
       end do
    end select

  end subroutine reconstructLmf

  subroutine reconstructRmf(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wRC)
    use mod_global_parameters
    use mod_limiter

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision, intent(out)   :: wRC(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision                :: rdw(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), dwC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer                         :: jxRmin1,jxRmin2,jxRmin3,jxRmax1,jxRmax2,&
       jxRmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, jxCmin1,&
       jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3, kxCmin1,kxCmin2,kxCmin3,&
       kxCmax1,kxCmax2,kxCmax3, kxRmin1,kxRmin2,kxRmin3,kxRmax1,kxRmax2,&
       kxRmax3, iw

    select case (type_limiter(block%level))
    case (limiter_mp5)
       call MP5limiterR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iLmin1,&
          iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wRC)
    case (limiter_weno5)
       call WENO5limiterR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wRC,1)
    case (limiter_wenoz5)
       call WENO5limiterR(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          iLmin1,iLmin2,iLmin3,iLmax1,iLmax2,iLmax3,idims,w,wRC,2)
    case default

       kxCmin1=ixImin1;kxCmin2=ixImin2;kxCmin3=ixImin3
       kxCmax1=ixImax1-kr(idims,1);kxCmax2=ixImax2-kr(idims,2)
       kxCmax3=ixImax3-kr(idims,3);
       kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
       kxRmin3=kxCmin3+kr(idims,3);kxRmax1=kxCmax1+kr(idims,1)
       kxRmax2=kxCmax2+kr(idims,2);kxRmax3=kxCmax3+kr(idims,3);

       wRC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,kxCmin3:kxCmax3,&
          1:nwflux)=w(kxRmin1:kxRmax1,kxRmin2:kxRmax2,kxRmin3:kxRmax3,&
          1:nwflux)

       jxRmin1=iLmin1+kr(idims,1);jxRmin2=iLmin2+kr(idims,2)
       jxRmin3=iLmin3+kr(idims,3);jxRmax1=iLmax1+kr(idims,1)
       jxRmax2=iLmax2+kr(idims,2);jxRmax3=iLmax3+kr(idims,3);
       ixCmax1=jxRmax1;ixCmax2=jxRmax2;ixCmax3=jxRmax3
       ixCmin1=iLmin1-kr(idims,1);ixCmin2=iLmin2-kr(idims,2)
       ixCmin3=iLmin3-kr(idims,3);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
       jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);

       do iw=1,nwflux
          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3)=w(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             jxCmin3:jxCmax3,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3,iw)
          call dwlimiter2(dwC,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,&
             type_limiter(block%level),rdw=rdw)

          wRC(iLmin1:iLmax1,iLmin2:iLmax2,iLmin3:iLmax3,iw)=wRC(iLmin1:iLmax1,&
             iLmin2:iLmax2,iLmin3:iLmax3,iw)-half*rdw(jxRmin1:jxRmax1,&
             jxRmin2:jxRmax2,jxRmin3:jxRmax3)
       end do
    end select

  end subroutine reconstructRmf

  subroutine centdiff4mf(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idimmin,idimmax,qtC,wCT,&
     qt,w,fC,dx1,dx2,dx3,x)
    ! Advance the flow variables from global_time to global_time+qdt within ixO^L by
    ! fourth order centered differencing in space
    ! for the dw/dt+dF_i(w)/dx_i=S type equation.
    ! wCT contains the time centered variables at time qtC for flux and source.
    ! w is the old value at qt on input and the new value at qt+qdt on output.
    use mod_global_parameters
    use mod_comm_lib, only: mpistop

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idimmin,idimmax
    double precision, intent(in) :: qdt, qtC, qt, dx1,dx2,dx3
    double precision :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision :: fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:ndir,1:ndim)

    double precision :: v(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ndim), f(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw) :: wLC, wRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)      :: vLC, vRC,cmaxLC,cmaxRC
    double precision :: dxinv(1:ndim)
    integer :: idims, idir, idirmin,ix1,ix2,ix3
    integer :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, hxOmin1,hxOmin2,&
       hxOmin3,hxOmax1,hxOmax2,hxOmax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,&
       ixCmax2,ixCmax3, jxCmin1,jxCmin2,jxCmin3,jxCmax1,jxCmax2,jxCmax3,&
        hxCmin1,hxCmin2,hxCmin3,hxCmax1,hxCmax2,hxCmax3, kxCmin1,kxCmin2,&
       kxCmin3,kxCmax1,kxCmax2,kxCmax3, kkxCmin1,kkxCmin2,kkxCmin3,kkxCmax1,&
       kkxCmax2,kkxCmax3, kkxRmin1,kkxRmin2,kkxRmin3,kkxRmax1,kkxRmax2,&
       kkxRmax3

    ! two extra layers are needed in each direction for which fluxes are added.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmin3=ixOmin3;ixmax1=ixOmax1
    ixmax2=ixOmax2;ixmax3=ixOmax3;
    do idims= idimmin,idimmax
       ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
       ixmin3=ixmin3-2*kr(idims,3);ixmax1=ixmax1+2*kr(idims,1)
       ixmax2=ixmax2+2*kr(idims,2);ixmax3=ixmax3+2*kr(idims,3);
    end do

    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImin3>ixmin3.or.ixImax1<ixmax1.or.&
       ixImax2<ixmax2.or.ixImax3<ixmax3) then
       call mpistop("Error in evolve_CentDiff4: Non-conforming input limits")
    end if
    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;dxinv(3)=-qdt/dx3;

    ! Add fluxes to w
    do idims= idimmin,idimmax
       b0i=idims
       ixmin1=ixOmin1-2*kr(idims,1);ixmin2=ixOmin2-2*kr(idims,2)
       ixmin3=ixOmin3-2*kr(idims,3);ixmax1=ixOmax1+2*kr(idims,1)
       ixmax2=ixOmax2+2*kr(idims,2);ixmax3=ixOmax3+2*kr(idims,3);
       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmin3=ixOmin3-kr(idims,3);hxOmax1=ixOmax1-kr(idims,1)
       hxOmax2=ixOmax2-kr(idims,2);hxOmax3=ixOmax3-kr(idims,3);

       ixCmin1=hxOmin1;ixCmin2=hxOmin2;ixCmin3=hxOmin3; ixCmax1=ixOmax1
       ixCmax2=ixOmax2;ixCmax3=ixOmax3;
       hxCmin1=ixCmin1-kr(idims,1);hxCmin2=ixCmin2-kr(idims,2)
       hxCmin3=ixCmin3-kr(idims,3);hxCmax1=ixCmax1-kr(idims,1)
       hxCmax2=ixCmax2-kr(idims,2);hxCmax3=ixCmax3-kr(idims,3);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmin3=ixCmin3+kr(idims,3);jxCmax1=ixCmax1+kr(idims,1)
       jxCmax2=ixCmax2+kr(idims,2);jxCmax3=ixCmax3+kr(idims,3);
       kxCmin1=ixCmin1+2*kr(idims,1);kxCmin2=ixCmin2+2*kr(idims,2)
       kxCmin3=ixCmin3+2*kr(idims,3);kxCmax1=ixCmax1+2*kr(idims,1)
       kxCmax2=ixCmax2+2*kr(idims,2);kxCmax3=ixCmax3+2*kr(idims,3);

       kkxCmin1=ixImin1;kkxCmin2=ixImin2;kkxCmin3=ixImin3
       kkxCmax1=ixImax1-kr(idims,1);kkxCmax2=ixImax2-kr(idims,2)
       kkxCmax3=ixImax3-kr(idims,3);
       kkxRmin1=kkxCmin1+kr(idims,1);kkxRmin2=kkxCmin2+kr(idims,2)
       kkxRmin3=kkxCmin3+kr(idims,3);kkxRmax1=kkxCmax1+kr(idims,1)
       kkxRmax2=kkxCmax2+kr(idims,2);kkxRmax3=kkxCmax3+kr(idims,3);
       wRC(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,kkxCmin3:kkxCmax3,&
          1:nwflux)=wCT(kkxRmin1:kkxRmax1,kkxRmin2:kkxRmax2,kkxRmin3:kkxRmax3,&
          1:nwflux)
       wLC(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,kkxCmin3:kkxCmax3,&
          1:nwflux)=wCT(kkxCmin1:kkxCmax1,kkxCmin2:kkxCmax2,kkxCmin3:kkxCmax3,&
          1:nwflux)

       call upwindLRmf(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixCmin1,&
          ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,ixCmin1,ixCmin2,ixCmin3,&
          ixCmax1,ixCmax2,ixCmax3,idims,wCT,wCT,wLC,wRC,x)

       ! Calculate velocities from upwinded values
       call getcmaxfff(wLC,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixCmin1,&
          ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxLC)
       call getcmaxfff(wRC,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixCmin1,&
          ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3,idims,cmaxRC)
       ! now take the maximum of left and right states
       vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
          ixCmin3:ixCmax3)=max(cmaxRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
          ixCmin3:ixCmax3),cmaxLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
          ixCmin3:ixCmax3))

       do idir=1,ndir
          ! Get non-transported flux
          call getfluxmf(wCT,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
             ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,idir,idims,f)
          ! Center flux to interface
          ! f_i+1/2= (-f_(i+2) +7 f_(i+1) + 7 f_i - f_(i-1))/12
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir,&
             idims)=(-f(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
             kxCmin3:kxCmax3)+7.0d0*(f(jxCmin1:jxCmax1,jxCmin2:jxCmax2,&
             jxCmin3:jxCmax3)+f(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ixCmin3:ixCmax3))-f(hxCmin1:hxCmax1,hxCmin2:hxCmax2,&
             hxCmin3:hxCmax3))/12.0d0
          ! add rempel dissipative flux, only second order version for now
          ! one could gradually reduce the dissipative flux to improve solutions
          ! for computing steady states (Keppens et al. 2012, JCP)
          fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir,&
             idims)=fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir,&
             idims)-mf_tvdlfeps*tvdlfeps*half*vLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3) *(wRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(idir))-wLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,ixCmin3:ixCmax3,mag(idir)))

          if (slab_uniform) then
             fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir,&
                idims)=dxinv(idims)*fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                ixCmin3:ixCmax3,idir,idims)
             ! result: f_(i+1/2)-f_(i-1/2) = [-f_(i+2)+8(f_(i+1)+f_(i-1))-f_(i-2)]/12
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mag(idir))+(fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                idir,idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                idir,idims))
          else
             fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir,&
                idims)=-qdt*block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                ixCmin3:ixCmax3,idims)*fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                ixCmin3:ixCmax3,idir,idims)
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mag(idir))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                mag(idir))+ (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                ixOmin3:ixOmax3,idir,idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                hxOmin3:hxOmax3,idir,idims))/block%dvolume(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,ixOmin3:ixOmax3)
          end if
       end do    !next idir
    end do       !next idims
    b0i=0

    if (.not.slab) call addgeometrymf(qdt,ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,&
       x)
    call divbclean(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)

  end subroutine centdiff4mf

  subroutine getdtfff_courant(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew)
    ! compute CFL limited dt (for variable time stepping)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), dtnew

    double precision :: courantmax, dxinv(1:ndim)
    double precision :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       alfven(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer :: idims

    dtnew=bigdouble
    courantmax=0.d0
    dxinv(1)=1.d0/dxlevel(1);dxinv(2)=1.d0/dxlevel(2)
    dxinv(3)=1.d0/dxlevel(3);

    do idims=1,ndim
       call getcmaxfff(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,cmax)
       cmax_mype = max(cmax_mype,maxval(cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          ixOmin3:ixOmax3)))
       if (.not.slab_uniform) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)/block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3,idims)
          courantmax=max(courantmax,maxval(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)))
       else
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)=cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)*dxinv(idims)
          courantmax=max(courantmax,maxval(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             ixOmin3:ixOmax3)))
       end if
    end do
    ! courantmax='max( c/dx)'
    if (courantmax>smalldouble)  dtnew=min(dtnew,mf_cc/courantmax)

  end subroutine getdtfff_courant

  subroutine getcmaxfff(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idims,cmax)
    use mod_global_parameters

    logical :: new_cmax,needcmin
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idims
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(out) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    ! calculate alfven speed
    if(B0field) then
      cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))+block%b0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,:,0))**2,dim=ndim+1)/w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_))
    else
      cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=sqrt(sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,mag(:))**2,dim=ndim+1)/w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_))
    endif
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)=cmax(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)+abs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,mom(idims)))

  end subroutine getcmaxfff

  !> Clean divergence of magnetic field by Janhunen's and Linde's source terms
  subroutine divbclean(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim),wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw),qdt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    integer :: idims, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3, ixpmin1,&
       ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3, i1,i2,i3, iside
    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       graddivb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
       bdivb(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir)

    ! Calculate div B
    ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmin3=ixOmin3-1;ixmax1=ixOmax1+1
    ixmax2=ixOmax2+1;ixmax3=ixOmax3+1;
    call get_divb(wCT,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixmin1,&
       ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,divb)

    ixpmin1=ixOmin1;ixpmin2=ixOmin2;ixpmin3=ixOmin3;ixpmax1=ixOmax1
    ixpmax2=ixOmax2;ixpmax3=ixOmax3;

    ! Add Linde's diffusive terms
    do idims=1,ndim
       ! Calculate grad_idim(divb)
       call gradient(divb,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
          ixpmin1,ixpmin2,ixpmin3,ixpmax1,ixpmax2,ixpmax3,idims,graddivb)

       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab_uniform) then
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)=graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)*mf_cdivb/(1.0d0/dxlevel(1)**2+&
             1.0d0/dxlevel(2)**2+1.0d0/dxlevel(3)**2)
       else
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)=graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3)*mf_cdivb /(1.0d0/block%dx(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
             1)**2+1.0d0/block%dx(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
             ixpmin3:ixpmax3,2)**2+1.0d0/block%dx(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,ixpmin3:ixpmax3,3)**2)
       end if
       ! B_idim += eta*grad_idim(divb)
       ! with Janhunen's term
       !w(ixp^S,mag(idims))=w(ixp^S,mag(idims))+&
       !      graddivb(ixp^S)-qdt*w(ixp^S,mom(idims))*divb(ixp^S)
       ! without Janjunen's term
       w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
          mag(idims))=w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,ixpmin3:ixpmax3,&
          mag(idims))+graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
          ixpmin3:ixpmax3)
    end do

  end subroutine divbclean

  subroutine addgeometrymf(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,wCT,w,x)
    ! Add geometrical source terms to w
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    !.. local ..
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    integer          :: iw
    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    select case (coordinate)
    case (cylindrical)
      if(phi_>0) then
        ! s[Bphi]=(Bphi*vr-Br*vphi)/radius
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(1)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(3)))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           bphi_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           bphi_)+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)
      end if
    case (spherical)
    
      ! s[b2]=(vr*Btheta-vtheta*Br)/r
      !       + cot(theta)*psi/r
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)= wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mag(1))
      if (B0field) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(1))*block%b0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,2,0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,mom(2))*block%b0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            ixOmin3:ixOmax3,1,0)
      end if
      ! Divide by radius and add to w
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(2))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mag(2))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)
   
      if(ndir==3) then
        ! s[b3]=(vr*Bphi-vphi*Br)/r
        ! -cot(theta)*(vphi*Btheta-vtheta*Bphi)/r
        tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)=wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mag(3)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mag(1)) -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,mag(3)))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,2)) /dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,2))
        if (B0field) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)=tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(1))*block%b0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,3,0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*block%b0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,1,0) -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(3))*block%b0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2,0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,mom(2))*block%b0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,3,0))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2)) /dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              ixOmin3:ixOmax3,2))
        end if
        ! Divide by radius and add to w
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(3))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           mag(3))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)
      end if
    end select

  end subroutine addgeometrymf

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,current)
    use mod_global_parameters
    use mod_geometry

    integer :: idirmin0
    integer :: ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idirmin,&
        ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3
    double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       1:nw)
    integer :: idir

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

  !> Calculate div B within ixO
  subroutine get_divb(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
    use mod_geometry
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision                   :: divb(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision                   :: bvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)

    bvec(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,:)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,mag(:))

    select case(typediv)
    case("central")
      call divvector(bvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
    case("limited")
      call divvectorS(bvec,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
         ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,divb)
    end select
  end subroutine get_divb

end module mod_magnetofriction
