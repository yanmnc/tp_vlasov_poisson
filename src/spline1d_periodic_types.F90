!**************************************************************
!  Copyright Euratom-CEA
!
!  See AUTHORS file for the list of authors.
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!===========================================================================
!> Types for ::spline1d_periodic_class 
!>
!> \date 2000-05-15
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module spline1d_periodic_types
  use prec_const

  implicit none

  !===========================================================================
  !> Instance of type ::spline1d_periodic_class : Definition of 1D 
  !>  periodic spline 
  !---------------------------------------------------------------------------
  type :: spline1d_periodic_t

    integer   :: n    ! dimension
    real(F64) :: h    ! discretis. step  
    ! extrema of parallel region
    real(F64) :: eta1b, eta2b, eta3b, eta4b
    ! 2x2 matrix determinant
    real(F64) :: det_deltab
    ! index in the array
    integer, dimension(:)  , pointer :: ipos1d            
    ! diagonal  terms
    real(F64), dimension(:)  , pointer :: Adiag   
    ! off-diagonal terms
    real(F64), dimension(:)  , pointer :: Aodiag  
    ! A^-1*gamma
    real(F64), dimension(:,:), pointer :: Am1gamma    
    ! spline coefficients    
    real(F64), dimension(:)  , pointer :: scoef    
    ! temporary arrays
    real(F64), dimension(:)  , pointer :: xprime
    real(F64), dimension(:)  , pointer :: rhs
    ! array used for integration
    integer  , dimension(1:6) :: indx
    real(F64), dimension(1:6) :: factor_int

  end type spline1d_periodic_t
  !---------------------------------------------------------------------------

end module spline1d_periodic_types
!---------------------------------------------------------------------------

