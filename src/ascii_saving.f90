!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!
!  Code VlasovPoisson : solving of the 2D Vlasov-Poisson
!    system for:
!      1) Landau damping study or
!      2) Plasma beam study
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!-----------------------------------------------------------
! file : ascii_saving.f90
!
!  used for saving the results in ASCII format
!-----------------------------------------------------------
module ascii_saving_module
  use globals, only : time_evol, &
    Phi_evol, E_evol, dens, velo, temp, &
    f1_evol, nbion_evol, Enkin_evol, Enpot_evol, &
    entropy_evol, L1_norm_evol, L2_norm_evol
  use geometry_class 

  implicit none

  !******************************
  contains
  !******************************

  !***********************************************************
  !   SAVING IN ASCII FORMAT
  !***********************************************************
  !--------------------------------------------------
  ! Saving of the simulation parameters 
  !--------------------------------------------------  
  subroutine ascii_param_write()
    use globals, only : v0, T0, epsilon, perturb, mode, Lx, vdiag_forf1Dx, &
      Phi_ext0, mode_ext, omega_ext, OhmsLaw, E0, nu_coeff

    integer     :: uout_param
    real(RKIND) :: kx, kx_ext, coef_OhmsLaw

    kx     = float(mode)*TWOPI/Lx
    kx_ext = float(mode_ext)*TWOPI/Lx

    coef_OhmsLaw = 0._RKIND
    if (OhmsLaw) then
      coef_OhmsLaw = 1._RKIND
    end if

    uout_param = 30
    open(uout_param,file='param_simu.dat',status='REPLACE',form='FORMATTED')
    write(uout_param,'(1pe12.3)') kx
    write(uout_param,'(1pe12.3)') perturb
    write(uout_param,'(1pe12.3)') v0
    write(uout_param,'(1pe12.3)') T0
    write(uout_param,'(1pe12.3)') epsilon
    write(uout_param,'(1pe12.3)') vdiag_forf1Dx
    write(uout_param,'(1pe12.3)') nu_coeff
    write(uout_param,'(1pe12.3)') Phi_ext0
    write(uout_param,'(1pe12.3)') kx_ext
    write(uout_param,'(1pe12.3)') omega_ext
    write(uout_param,'(1pe12.3)') coef_OhmsLaw
    write(uout_param,'(1pe12.3)') E0
    close(uout_param)
  end subroutine ascii_param_write


  !--------------------------------------------------
  ! Saving of the electric potential 1D : Phi(x) 
  !  (one file per diagnostic)
  !--------------------------------------------------  
  subroutine ascii_Phi1D_write(idiag_num,geom,Phi1D)
    integer                    , intent(in) :: idiag_num
    type(geometry)             , intent(in) :: geom    
    real(RKIND), dimension(0:) , intent(in) :: Phi1D
    
    integer           :: i1
    integer           :: uout_Phi
    character(LEN=50) :: Phi1D_file_name

    uout_Phi        = 40
    Phi1D_file_name = "Phi1D_    .dat"//char(0)
    write(Phi1D_file_name(7:10),'(i4.4)') idiag_num
    open(uout_Phi,file=Phi1D_file_name,status='REPLACE', form = 'FORMATTED') 
    do i1 = 0,geom%Nx
      write(uout_Phi,'((1pe12.3),(1pe20.12))') geom%xgrid(i1), Phi1D(i1)
    end do
    close(uout_Phi)
  end subroutine ascii_Phi1D_write


  !--------------------------------------------------
  ! Saving of the 1D fluid moments: dens(x), velo(x), temp(x)
  !  (one file per diagnostic)
  !--------------------------------------------------  
  subroutine ascii_fluid_moments_write(idiag_num,geom,&
    dens,velo,temp,flux)
    integer                    , intent(in) :: idiag_num
    type(geometry)             , intent(in) :: geom    
    real(RKIND), dimension(0:) , intent(in) :: dens, velo, temp, flux
    
    integer           :: i1
    integer           :: uout_dens, uout_velo, uout_temp, uout_flux
    character(LEN=50) :: dens_file_name, velo_file_name
    character(LEN=50) :: temp_file_name, flux_file_name

    uout_dens      = 50
    uout_velo      = 51
    uout_temp      = 52
    uout_flux      = 53
    dens_file_name = "dens_    .dat"//char(0)
    velo_file_name = "velo_    .dat"//char(0)
    temp_file_name = "temp_    .dat"//char(0)
    flux_file_name = "flux_    .dat"//char(0)
    write(dens_file_name(6:9),'(i4.4)') idiag_num
    write(velo_file_name(6:9),'(i4.4)') idiag_num
    write(temp_file_name(6:9),'(i4.4)') idiag_num
    write(flux_file_name(6:9),'(i4.4)') idiag_num
    open(uout_dens,file=dens_file_name,status='REPLACE', form = 'FORMATTED') 
    open(uout_velo,file=velo_file_name,status='REPLACE', form = 'FORMATTED') 
    open(uout_temp,file=temp_file_name,status='REPLACE', form = 'FORMATTED') 
    open(uout_flux,file=flux_file_name,status='REPLACE', form = 'FORMATTED') 
    do i1 = 0,geom%Nx
      write(uout_dens,'((1pe12.3),(1pe20.12))') geom%xgrid(i1), dens(i1)
      write(uout_velo,'((1pe12.3),(1pe20.12))') geom%xgrid(i1), velo(i1)
      write(uout_temp,'((1pe12.3),(1pe20.12))') geom%xgrid(i1), temp(i1)
      write(uout_flux,'((1pe12.3),(1pe20.12))') geom%xgrid(i1), flux(i1)
    end do
    close(uout_dens)
    close(uout_velo)
    close(uout_temp)
    close(uout_flux)
  end subroutine ascii_fluid_moments_write
  
  
  !--------------------------------------------------
  ! Saving of distribution 1D profile of f(x,v) 
  !  (one file per diagnostic)
  !   - f(x,v=v_diag)  for vdiag chosen
  !   - f(x=x_diag,v) for xdiag chosen
  !--------------------------------------------------  
  subroutine ascii_f1D_write(idiag_num,geom,f2D)
    use globals     , only : vdiag_forf1Dx
    use utils_module, only : locate
    integer                      , intent(in) :: idiag_num
    type(geometry)               , intent(in) :: geom    
    real(RKIND), dimension(0:,0:), intent(in) :: f2D
    
    integer           :: i1, i2
    integer           :: ixdiag, ivdiag
    integer           :: uout_f1Dx, uout_f1Dv
    character(LEN=50) :: f1Dx_filename, f1Dv_filename
    character(LEN=50), parameter :: subp = " ascii_f1D_write "//char(0)

    ixdiag = int(geom%Nx/2)
    call locate(vdiag_forf1Dx,geom%vgrid,geom%Nv,geom%dv,subp,ivdiag)

    !*** Saving of f(x,v=v_diag) ***
    uout_f1Dx     = 41
    f1Dx_filename = "f1Dx_    .dat"//char(0)
    write(f1Dx_filename(6:9),'(i4.4)') idiag_num
    open(uout_f1Dx,file=f1Dx_filename,status='REPLACE',form='FORMATTED') 
    do i1 = 0,geom%Nx
      write(uout_f1Dx,'((1pe12.3),(1pe20.12))') geom%xgrid(i1), &
        f2D(i1,ivdiag)
    end do
    close(uout_f1Dx)

    !*** Saving of f(x=x_diag,v) ***
    uout_f1Dv     = 42
    f1Dv_filename = "f1Dv_    .dat"//char(0)
    write(f1Dv_filename(6:9),'(i4.4)') idiag_num
    open(uout_f1Dv,file=f1Dv_filename,status='REPLACE',form='FORMATTED') 
    do i2 = 0,geom%Nv
      write(uout_f1Dv,'((1pe12.3),(1pe20.12))') geom%vgrid(i2), &
        f2D(ixdiag,i2)
    end do
    close(uout_f1Dv)
  end subroutine ascii_f1D_write


  !--------------------------------------------------
  ! Saving of distribution 2D : f(x,v) 
  !  (one file per diagnostic)
  !--------------------------------------------------  
  subroutine ascii_f2D_write(idiag_num,geom,f2D)
    integer                       , intent(in) :: idiag_num
    type(geometry)                , intent(in) :: geom    
    real(RKIND), dimension(0:,0:) , intent(in) :: f2D
    
    integer           :: i1, i2
    integer           :: uout_f2D
    character(LEN=50) :: f2D_file_name  

    uout_f2D      = 41
    f2D_file_name = "f2D_    .dat"//char(0)
    write(f2D_file_name(5:8),'(i4.4)') idiag_num
    open(uout_f2D,file=f2D_file_name,status='REPLACE', form = 'FORMATTED') 
    do i2 = 0,geom%Nv
      do i1 = 0,geom%Nx
        write(uout_f2D,'(2(1pe12.3),(1pe20.12))') &
          geom%xgrid(i1), geom%vgrid(i2), f2D(i1,i2)
      end do
      write(uout_f2D,*) " "
    end do
    close(uout_f2D)
  end subroutine ascii_f2D_write

 
  !---------------------------------------------- 
  ! ASCII results saving 
  !----------------------------------------------   
  subroutine ascii_results_saving(geom,nb_save)    
    type(geometry)    , intent(in) :: geom
    integer           , intent(in) :: nb_save
     
    character(LEN=50) :: resu_file_name    
    integer           :: uout_resu
    integer           :: isave

    uout_resu = 42
    write(resu_file_name,'(A)')"VlasovPoiss_res.dat"//char(0)
    open(uout_resu,file = resu_file_name, &
      status='REPLACE',form ='FORMATTED') 
    !VG!write(uout_resu,'((A12),6(A20))') 'time_diag', 'nbions', &
    !VG!  'Enkin', 'Enpot', 'entropy', 'L1_norm', 'L2_norm'
    do isave = 0,nb_save
      write(uout_resu,'((1pe12.3),6(1pe20.12))') time_evol(isave) , &
        nbion_evol(isave), Enkin_evol(isave), Enpot_evol(isave), &
        entropy_evol(isave), L1_norm_evol(isave), L2_norm_evol(isave)
    end do
    close(uout_resu)
  end subroutine ascii_results_saving
end module ascii_saving_module

