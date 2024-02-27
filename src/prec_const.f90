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
!---------------------------------------------
! file : prec_const.f90
!
! Constant and parameter initialisation
!---------------------------------------------

module prec_const
  implicit none

  !> Formatting for real simple precision
  integer, parameter :: F32 = SELECTED_REAL_KIND(p=6,r=37)
  !> Formatting for real double precision 
  integer, parameter :: F64 = SELECTED_REAL_KIND(p=13,r=200)
  !> Formatting for integer
  integer, parameter :: I32 = SELECTED_INT_KIND(9)
  !> Formatting for integer
  integer, parameter :: I64 = SELECTED_INT_KIND(18)

  !*** Precision for single ***
  integer, parameter :: SKIND = SELECTED_REAL_KIND(p=6,r=37)
  !*** Precision for real
  integer, parameter :: RKIND = SELECTED_REAL_KIND(p=13,r=200)
  !*** Precision for complex
  integer, parameter :: CKIND = RKIND
  
  !*** Some useful constants
  real(RKIND), parameter :: ZE    = 0.0_RKIND
  real(RKIND), parameter :: HF    = 0.5_RKIND
  real(RKIND), parameter :: ON    = 1.0_RKIND
  real(RKIND), parameter :: TW    = 2.0_RKIND
  real(RKIND), parameter :: TH    = 3.0_RKIND
  real(RKIND), parameter :: FO    = 4.0_RKIND
  real(RKIND), parameter :: &
    PI    = 3.141592653589793238462643383279502884197_RKIND
  real(RKIND), parameter :: &
    TWOPI = 6.283185307179586476925286766559005768394_RKIND

  !*** complex value
  complex(CKIND) :: ci = (0._RKIND,1._RKIND)
end module prec_const

