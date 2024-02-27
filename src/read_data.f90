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

!----------------------------------------
! file : read_data.f90
!
! used for the reading of :
! - mesh datas
! - the equilibrium datas
!----------------------------------------

!------------------------------------
! mesh data reading
!------------------------------------
subroutine mesh_data
  use globals, only : Nx, Nv, Lx, &
    Lv, x0, nb_vthi 
  implicit none
  
  namelist /MESH/ Nx, Nv, Lx, Lv, x0, nb_vthi

  read(*,MESH)
  write(*,MESH)
end subroutine mesh_data


!----------------------------------------
! equilibrium data reading
!----------------------------------------
subroutine equil_data
  use prec_const
  use globals, only : BC_v_left, BC_v_right, &
      v0, T0, symmetric, epsilon, perturb, perturb_choice, mode, &
      Phi_ext0, mode_ext, omega_ext, OhmsLaw, E0
  implicit none
  
  namelist /EQUIL/ BC_v_left, BC_v_right, &
      v0, T0, symmetric, epsilon, perturb, perturb_choice, mode, &
      Phi_ext0, mode_ext, omega_ext, OhmsLaw, E0
  
  read(*,EQUIL)
  write(*,EQUIL)
end subroutine equil_data


!----------------------------------------
! RHS data
!----------------------------------------
subroutine RHS_data
  use prec_const
  use globals, only : RHS_only, nu_coeff, diff_coeff
  implicit none
  
  namelist /RHS/ RHS_only, nu_coeff, diff_coeff
  

  read(*,RHS)
  write(*,RHS)
end subroutine RHS_data


!----------------------------------------
! algorithm data reading
!----------------------------------------
subroutine algorithm_data
  use prec_const
  use globals, only : time_scheme, deltat, nbiter
  implicit none
  
  namelist /ALGORITHM/ time_scheme, deltat, nbiter
  
  read(*,ALGORITHM)
  write(*,ALGORITHM)
end subroutine algorithm_data


!----------------------------------------
! output data reading
!----------------------------------------
subroutine output_data
  use prec_const
  use globals, only : nbstep, vdiag_forf1Dx, f2D_saving
  implicit none
  
  namelist /OUTPUT/ nbstep, vdiag_forf1Dx, f2D_saving
  
  read(*,OUTPUT)
  write(*,OUTPUT)
end subroutine output_data


!----------------------------------------
!  reading of the input file .dat
!----------------------------------------
subroutine read_input
  implicit none
  
  call mesh_data
  call equil_data
  call RHS_data
  call algorithm_data
  call output_data
end subroutine read_input


!------------------------------------------
! write parameters in the text result file
!------------------------------------------
subroutine write_parameters
  use globals  
  implicit none

  namelist /MESH/ Nx, Nv, Lx, Lv, x0, nb_vthi
  namelist /EQUIL/ v0, T0, symmetric, epsilon, perturb, mode, &
    Phi_ext0, mode_ext, omega_ext, OhmsLaw, E0
  namelist /RHS/ RHS_only, nu_coeff, diff_coeff
  namelist /ALGORITHM/ time_scheme, deltat, nbiter
  namelist /OUTPUT/ nbstep, vdiag_forf1Dx, f2D_saving
  
  integer, parameter :: upar = 23
  character(LEN=50)  :: parameter_file 

  write(parameter_file,'("VlasovPoiss2D.prm")')
  open(upar,file = parameter_file,status = 'REPLACE', &
    form = 'FORMATTED')
  write(upar,MESH)
  write(upar,EQUIL)
  write(upar,RHS)
  write(upar,ALGORITHM)
  write(upar,OUTPUT)
  close(upar)
end subroutine write_parameters
