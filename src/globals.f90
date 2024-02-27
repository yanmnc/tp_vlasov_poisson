!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!
!  Code VlasovPoisson : solving of the 2D Vlasov-Poisson
!    system for:
!      1) Landau damping study or
!      2) Plasma beam study
 
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!-------------------------------------------------------------------
! file : globals.f90
!
! global variable definition
!-------------------------------------------------------------------
module globals
  use prec_const
  
  implicit none
  
  !**************************************************************
  !    INPUT DATA
  !**************************************************************
  !--------------------------------------------------------------
  ! GEOMETRY DATA:
  !  Variables used to define the 2D (x,v) geometry:
  !
  !    - Nx      : (Nx+1) corresponds to the number of points 
  !                 in x direction 
  !    - Nv      : (Nv+1) corresponds to the number of points 
  !                 in v direction 
  !    - Lx      : length in x direction
  !    - Lv      : length in v direction
  !    - x0      : domain in x direction is defined 
  !                 between [x0,x0+Lx]
  !    - nb_vthi : domain in v direction is defined 
  !                 between [-nb_vthi,+nb_vthi]
  !--------------------------------------------------------------
  integer    , save :: Nx      = 50
  integer    , save :: Nv      = 50
  real(RKIND), save :: Lx      = 1._RKIND
  real(RKIND), save :: Lv      = 1._RKIND
  real(RKIND), save :: x0      = 0._RKIND
  real(RKIND), save :: nb_vthi = 3._RKIND
  
  !--------------------------------------------------------------
  ! EQUILIBRIUM DATA:
  !  Variables used to define:
  !   1) The equilibrium distribution function fM as:
  !        fM(v) = f1(v) + f2(v) 
  !       with 
  !        f1(v) = (1-epsilon)/sqrt(2*PI)*exp(-v**2/2)
  !        f2(v) = epsilon/sqrt(2*PI*T0)*exp(-(v-v0)**2/(2*T0))
  ! 
  !   2) The perturbation such that:
  !       f(x,v) = fM(v)(1+perturb*cos(kx*x))
  !      with kx = 2*pi*mode/Lx
  !
  ! Rk: In the Landau damping case : epsilon=0 
  !      such that fM(v) = 1/(sqrt(2*PI))*exp(-v**2/2)  
  !--------------------------------------------------------------
  integer    , save :: BC_v_left      = 0
  integer    , save :: BC_v_right     = 0
  real(RKIND), save :: v0             = 0._RKIND
  real(RKIND), save :: T0             = 0.5_RKIND
  logical    , save :: symmetric      = .false.
  real(RKIND), save :: epsilon        = 0.5_RKIND
  real(RKIND), save :: perturb        = 0.001_RKIND  
  integer    , save :: perturb_choice = 1
  integer    , save :: mode           = 1
  real(RKIND), save :: Phi_ext0       = 0.0_RKIND  
  integer    , save :: mode_ext       = 3
  real(RKIND), save :: omega_ext      = 0.0_RKIND  
  logical    , save :: OhmsLaw        = .false.
  real(RKIND), save :: E0             = 0.0_RKIND  
  !--------------------------------------------------------------
  ! KROOK OPERATOR
  !  Used for Krook operator definition, for solving of 
  !   df/dt = nu_coeff*(f-fM)
  !--------------------------------------------------------------
  logical    , save :: RHS_only   = .false.
  real(RKIND), save :: nu_coeff   = 0._RKIND
  real(RKIND), save :: diff_coeff = 0._RKIND 

  !--------------------------------------------------------------
  ! ALGORITHM DATA
  !  Variables used to define the time scheme.
  !
  !   - time_scheme : equal to 'PC' if predictor-corrector and
  !                    'LF' if leap-frog time scheme
  !   - deltat      : time step of the numerical algorithm
  !   - nbiter      : number of iterations of the algorithm
  !--------------------------------------------------------------
  character(LEN=5), save :: time_scheme = 'LF'
  real(RKIND)     , save :: deltat      = 0.1_RKIND
  integer         , save :: nbiter      = 1

  !--------------------------------------------------------------
  ! OUTPUT DATA
  !  Variables used for the output savings.
  !
  !   - nbstep        : diagnostics are saved every 'nbstep'
  !   - vdiag_forf1Dx : position in v for 1D distribution 
  !                      function saving, i.e 
  !                      f(x,vdiag_forf1Dx) will be saved in the 
  !                      ASCII files 'f1Dx_<num_diag>.dat'
  !   - f2D_saving    : if equal to .true., the 2D distribution
  !                      function f(x,v) will be saved in the 
  !                      ASCII files 'f2D_<num_diag>.dat'
  !
  !  Rk: f(xmiddle,v) will be automatically saved in the 
  !      ASCII files 'f1Dv_<num_diag>.dat' where xmiddle 
  !      corresponds to the middle of the box in x, i.e 
  !      xmiddle = x0+Lx/2
  !--------------------------------------------------------------
  integer    , save :: nbstep        = 1
  real(RKIND), save :: vdiag_forf1Dx = 1._RKIND
  logical    , save :: f2D_saving    = .false.


  !**************************************************************
  !   OTHER USEFUL VARIABLES 
  !**************************************************************
  integer    , save :: iter_glob = 0
  integer    , save :: idiag     = 0
  real(RKIND)       :: init_time

  !--------------------------------------------------------------
  ! Array for Matlab saving
  !--------------------------------------------------------------
  real(RKIND), dimension(:)    , pointer :: time_evol
  real(RKIND), dimension(:,:)  , pointer :: Phi_evol
  real(RKIND), dimension(:,:)  , pointer :: E_evol
  real(RKIND), dimension(:)    , pointer :: dens, velo, temp, flux
  real(RKIND), dimension(:,:,:), pointer :: f1_evol 
  real(RKIND), dimension(:)    , pointer :: Enkin_evol
  real(RKIND), dimension(:)    , pointer :: Enpot_evol
  real(RKIND), dimension(:)    , pointer :: nbion_evol
  real(RKIND), dimension(:)    , pointer :: entropy_evol
  real(RKIND), dimension(:)    , pointer :: L1_norm_evol
  real(RKIND), dimension(:)    , pointer :: L2_norm_evol
end module globals



