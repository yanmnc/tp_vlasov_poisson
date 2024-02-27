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
!-------------------------------------------------
! file : resol_steps.f90
! 
!  Each step of the global resolution algorithm
!-------------------------------------------------
module resol_steps_module
  use globals
  use geometry_class
  use fdistribu2d_class
  use efield_module
  implicit none
      
  !******************************
  contains
  !******************************

  !------------------------------------------------------------
  ! Solve the Vlasov equation for the distribution function f
  !  using the field E, with the shift sequence : 
  !     x/2,v,x/2
  !------------------------------------------------------------
  subroutine solve_vlasov_xvx(geom,dt,E,f)
    use globals, only : nu_coeff, diff_coeff, RHS_only
    use advec1D_SL_module
    use RHS_module
    type(geometry)                  , intent(in)    :: geom
    real(RKIND)                     , intent(in)    :: dt
    real(RKIND)      , dimension(0:), intent(in)    :: E
    type(fdistribu2d)               , intent(inout) :: f
    
    real(RKIND) :: half_dt
    
    half_dt = 0.5_RKIND*dt
    
    if (.not.RHS_only) then
      if (nu_coeff.ne.0) &
        call solve_Krook_operator(geom,f,nu_coeff,half_dt)
      if (diff_coeff.ne.0) &
        call solve_Diffusion_operator(f)
      call advec1D_x(geom,half_dt,f)
      call advec1D_v(geom,dt,E,f)
      call advec1D_x(geom,half_dt,f)
      if (diff_coeff.ne.0) &
        call solve_Diffusion_operator(f)
      if (nu_coeff.ne.0) &
        call solve_Krook_operator(geom,f,nu_coeff,half_dt)
    else
      if (diff_coeff.ne.0) &
        call solve_Diffusion_operator(f)
      if (nu_coeff.ne.0) &
        call solve_Krook_operator(geom,f,nu_coeff,dt)
      if (diff_coeff.ne.0) &
        call solve_Diffusion_operator(f)
    end if
  end subroutine solve_vlasov_xvx


  !------------------------------------------------------------
  ! Solve the Vlasov equation for the distribution function f
  !  using the field E, with the shift sequence :
  !     v/2,x,v/2 
  !------------------------------------------------------------
  subroutine solve_vlasov_vxv(geom,dt,E,f)
    use globals, only : nu_coeff, diff_coeff, RHS_only
    use advec1D_SL_module
    use RHS_module
    type(geometry)            , intent(in)    :: geom
    real(RKIND)               , intent(in)    :: dt
    real(RKIND), dimension(0:), intent(in)    :: E
    type(fdistribu2d)         , intent(inout) :: f

    real(RKIND) :: half_dt
    
    half_dt = 0.5_RKIND*dt

    if (.not.RHS_only) then
      if (nu_coeff.ne.0) &
        call solve_Krook_operator(geom,f,nu_coeff,half_dt)
      if (diff_coeff.ne.0) &
        call solve_Diffusion_operator(f)
      call advec1D_v(geom,half_dt,E,f)
      call advec1D_x(geom,dt,f)
      call advec1D_v(geom,half_dt,E,f)
      if (diff_coeff.ne.0) &
        call solve_Diffusion_operator(f)
      if (nu_coeff.ne.0) &
        call solve_Krook_operator(geom,f,nu_coeff,half_dt)
    else
      if (diff_coeff.ne.0) &
        call solve_Diffusion_operator(f)
      if (nu_coeff.ne.0) &
        call solve_Krook_operator(geom,f,nu_coeff,dt)
      if (diff_coeff.ne.0) &
        call solve_Diffusion_operator(f)
    end if

  end subroutine solve_vlasov_vxv


  !********************************************************
  !********************************************************
  !    LEAP-FROG ALGORITHM
  !********************************************************
  !********************************************************
  !---------------------------------------------------------
  ! step 0 : (step only used for the initialisation)
  !  - f1 initialisation,
  !  - computation of the corresponding Phi1 and E1,
  !  - f2 initialisation,
  !  - computation of the corresponding Phi2 and E2.
  !--------------------------------------------------------- 
  subroutine step0_LF(geom,dt,f1,Phi1,E1,f2,Phi2,E2)
    type(geometry)                  , intent(in)    :: geom
    real(RKIND)                     , intent(in)    :: dt
    type(fdistribu2d)               , intent(inout) :: f1
    real(RKIND)      , dimension(0:), intent(inout) :: Phi1
    real(RKIND)      , dimension(0:), intent(inout) :: E1
    type(fdistribu2d)               , intent(inout) :: f2
    real(RKIND)      , dimension(0:), intent(inout) :: Phi2
    real(RKIND)      , dimension(0:), intent(inout) :: E2

    print*,'**********************************************************'
    print*,'******      RESOLUTION WITH LEAP-FROG ALGORITHM     ******'
    print*,'**********************************************************'
    print*,' '

    !-> initialisation of f1 at time t=0
    call f_initialisation(f1,geom)    
    !-> computation of Phi1(t=0)
    call compute_Efield(geom,1,f1,Phi1,E1)
    !-> solve Vlasov on dt/2 by using Phi1(0) -> f2(dt/2)
    f2%values = f1%values
    call solve_vlasov_vxv(geom,0.5_RKIND*dt,E1,f2)
    !*** compute Phi2(t=dt/2) ***
    call compute_Efield(geom,1,f2,Phi2,E2)
  end subroutine step0_LF


  !------------------------------------------------------------
  ! global iteration of leap_frog
  !------------------------------------------------------------
  subroutine leap_frog_iteration(geom,iter,dt,f_tmp,f1,Phi1,E1,f2,Phi2,E2)
    type(geometry)                  , intent(in)    :: geom
    integer                         , intent(in)    :: iter
    real(RKIND)                     , intent(in)    :: dt
    type(fdistribu2d)               , intent(inout) :: f_tmp
    type(fdistribu2d)               , intent(inout) :: f1
    real(RKIND)      , dimension(0:), intent(inout) :: Phi1
    real(RKIND)      , dimension(0:), intent(inout) :: E1
    type(fdistribu2d)               , intent(inout) :: f2
    real(RKIND)      , dimension(0:), intent(inout) :: Phi2
    real(RKIND)      , dimension(0:), intent(inout) :: E2

    if (mod(iter,30).ne.0) then
      !*** f1(tn+1) computed with f1(tn) on dt using E2(tn+1/2) ***
      call solve_vlasov_xvx(geom,dt,E2,f1)
      !*** computation of E1(tn+1) ***
      call compute_Efield(geom,iter,f1,Phi1,E1)
      !*** f2(tn+3/2) computed with f2(tn+1/2) on dt using E1(tn+1) ***
      call solve_vlasov_vxv(geom,dt,E1,f2)
      !*** computation of E2(tn+3/2) ***
      call compute_Efield(geom,iter,f2,Phi2,E2)
    else
      print*,'---> leap-frog average'
      !*** average f1(tn+1)=(f2(tn+1/2)+f2(tn+3/2))/2 ***
      call solve_vlasov_xvx(geom,dt,E2,f1)
      call compute_Efield(geom,iter,f1,Phi1,E1)
      f_tmp%values = f2%values
      call solve_vlasov_vxv(geom,dt,E1,f2)
      call compute_Efield(geom,iter,f2,Phi2,E2)
      f1%values = 0.5_RKIND*(f_tmp%values + f2%values)
    end if
  end subroutine leap_frog_iteration


  !********************************************************
  !********************************************************
  !    PREDICTOR-CORRECTOR ALGORITHM
  !********************************************************
  !********************************************************
  !---------------------------------------------------------
  ! step 0 : (step only used for the initialisation)
  !  - f1 initialisation,
  !  - computation of the corresponding Phi1 and E1,
  !--------------------------------------------------------- 
  subroutine step0_PC(geom,dt,f1,Phi1,E1)
    type(geometry)                  , intent(in)    :: geom
    real(RKIND)                     , intent(in)    :: dt
    type(fdistribu2d)               , intent(inout) :: f1
    real(RKIND)      , dimension(0:), intent(inout) :: Phi1
    real(RKIND)      , dimension(0:), intent(inout) :: E1

    print*,'*************************************************************'
    print*,'******  RESOLUTION WITH PREDICTOR-CORRECTOR ALGORITHM  ******'
    print*,'*************************************************************'
    print*,' '

    !-> initialisation of f1 at time t=0
    call f_initialisation(f1,geom)    
    !-> computation of Phi1(t=0)
    call compute_Efield(geom,1,f1,Phi1,E1)
  end subroutine step0_PC


  !------------------------------------------------------------
  ! global iteration of predictor-corrector
  !------------------------------------------------------------
  subroutine predcorr_iteration(geom,iter,dt,f_tmp,f1,Phi1,E1)
    type(geometry)                  , intent(in)    :: geom
    integer                         , intent(in)    :: iter
    real(RKIND)                     , intent(in)    :: dt
    type(fdistribu2d)               , intent(inout) :: f_tmp
    type(fdistribu2d)               , intent(inout) :: f1
    real(RKIND)      , dimension(0:), intent(inout) :: Phi1
    real(RKIND)      , dimension(0:), intent(inout) :: E1

    real(RKIND) :: half_dt

    half_dt = 0.5_RKIND*dt

    !*** prediction on a dt/2 ***
    f_tmp%values = f1%values
    call solve_vlasov_xvx(geom,half_dt,E1,f_tmp)
    !--> computation of E1(tn+1/2)
    call compute_Efield(geom,iter,f_tmp,Phi1,E1)

    !*** correction on a dt ***
    call solve_vlasov_vxv(geom,dt,E1,f1)
    !--> computation of E1(tn+1)
    call compute_Efield(geom,iter,f1,Phi1,E1)
  end subroutine predcorr_iteration

end module resol_steps_module









