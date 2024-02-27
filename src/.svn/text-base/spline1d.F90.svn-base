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
!> Common functions for 1D cubic spline interpolation both for 
!>  natural (non-periodic) and periodic boundary conditions
!>
!> \date 2000-05-15
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module spline1d_mod
  use prec_const

  implicit none

  private
  public :: spline1d_common_basis, &
      spline1d_common_basisderiv, spline1d_common_basis2deriv, &
      spline1d_common_integration_coef_BC

  !*** Computation of the cubic splines basis ***
  interface spline1d_common_basis
    module procedure spline1d_common_basis_point, spline1d_common_basis_knot
  end interface

  !*** Computation of the first derivate of the ***
  !***  cubic splines basis                     ***
  interface spline1d_common_basisderiv
    module procedure spline1d_common_basisderiv_point, &
        spline1d_common_basisderiv_knot
  end interface

  !*** Computation of the second derivate of the ***
  !***  cubic splines basis                      ***
  interface spline1d_common_basis2deriv
    module procedure spline1d_common_basis2deriv_point, &
        spline1d_common_basis2deriv_knot
  end interface

contains

  !===========================================================================
  !> Computation of the cubic spline basis 
  !>   forall point of the space
  !>
  !> \remark Available for periodic or natural function
  !---------------------------------------------------------------------------
  subroutine spline1d_common_basis_point( xprev, xstar, xnext, h, sbase )

    real(F64)                 , intent(in)    :: xprev,xstar
    real(F64)                 , intent(in)    :: xnext,h
    real(F64), dimension(-1:2), intent(inout) :: sbase

    real(F64) :: coeff, d_prev, d_next

    coeff     = 1._F64/(h*h*h)
    d_prev    = xstar-xprev
    d_next    = xnext-xstar
    sbase(2)  = coeff*d_prev*d_prev*d_prev
    sbase(1)  = 1._F64+coeff*(3._F64*d_prev* &
        (h*h+h*d_prev-d_prev*d_prev))
    sbase(0)  = 1._F64+coeff*(3._F64*d_next* &
        (h*h+h*d_next-d_next*d_next))
    sbase(-1) = coeff*d_next*d_next*d_next  

  end subroutine spline1d_common_basis_point
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Computation of the cubic spline basis 
  !>  for a knot of the mesh
  !>
  !> \remark Available for periodic or natural function
  !---------------------------------------------------------------------------
  subroutine spline1d_common_basis_knot( num, nb_knots, sbase )

    integer                   , intent(in)    :: num
    integer                   , intent(in)    :: nb_knots
    real(F64), dimension(-1:2), intent(INOUT) :: sbase

    if (num.lt.nb_knots) then
      sbase(-1) = 1._F64
      sbase(0)  = 4._F64
      sbase(1)  = 1._F64
      sbase(2)  = 0._F64
    else
      sbase(-1) = 0._F64
      sbase(0)  = 1._F64
      sbase(1)  = 4._F64
      sbase(2)  = 1._F64
    endif

  end subroutine spline1d_common_basis_knot
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Computation of the second derivative of the cubic spline 
  !>    basis forall point of the space
  !>
  !> \remark Available for periodic or natural function
  !---------------------------------------------------------------------------
  subroutine spline1d_common_basisderiv_point( xprev, xstar, xnext, h, sbase )

    real(F64)                 , intent(in)   :: xprev
    real(F64)                 , intent(in)   :: xstar
    real(F64)                 , intent(in)   :: xnext
    real(F64)                 , intent(in)   :: h
    real(F64), dimension(-1:2), intent(inout):: sbase

    real(F64) :: coeff, d_prev, d_next

    coeff     = 3._F64/(h*h*h)
    d_prev    = xstar-xprev
    d_next    = xnext-xstar
    sbase(2)  = coeff*d_prev*d_prev
    sbase(1)  = coeff * &
        (h*h+2._F64*h*d_prev-3._F64*d_prev*d_prev)
    sbase(0)  = -1._F64*coeff * &
        (h*h+2._F64*h*d_next-3._F64*d_next*d_next)
    sbase(-1) = -1._F64*coeff*d_next*d_next
    
  end subroutine spline1d_common_basisderiv_point
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Computation of the first derivate of 
  !>  the cubic spline basis for a knot of the mesh
  !>
  !> \remark Available for periodic or natural function
  !---------------------------------------------------------------------------
  subroutine spline1d_common_basisderiv_knot( num, nb_knots, h, sbase_prime )

    integer                   , intent(in)    :: num
    integer                   , intent(in)    :: nb_knots
    real(F64)                 , intent(in)    :: h
    real(F64), dimension(-1:2), intent(INOUT) :: sbase_prime

    if (num.lt.nb_knots) then
      sbase_prime(-1) = -3._F64/h
      sbase_prime(0)  = 0._F64
      sbase_prime(1)  = 3._F64/h
      sbase_prime(2)  = 0._F64
    else
      sbase_prime(-1) = 0._F64
      sbase_prime(0)  = -3._F64/h
      sbase_prime(1)  = 0._F64
      sbase_prime(2)  = 3._F64/h
    endif

  end subroutine spline1d_common_basisderiv_knot
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Computation of the second derivative of the cubic spline 
  !>   basis forall point of the space
  !>
  !> \remark Available for periodic or natural function
  !---------------------------------------------------------------------------
  subroutine spline1d_common_basis2deriv_point( xprev, xstar, xnext, h, &
      sbase_second )

    real(F64)                 , intent(in)   :: xprev
    real(F64)                 , intent(in)   :: xstar
    real(F64)                 , intent(in)   :: xnext
    real(F64)                 , intent(in)   :: h
    real(F64), dimension(-1:2), intent(inout):: sbase_second

    real(F64) :: coeff, d_prev, d_next

    coeff            = 6._F64/(h*h*h)
    d_prev           = xstar-xprev
    d_next           = xnext-xstar
    sbase_second(2)  = coeff*d_prev
    sbase_second(1)  = coeff*(h-3._F64*d_prev)
    sbase_second(0)  = coeff*(h-3._F64*d_next)
    sbase_second(-1) = coeff*d_next
    
  end subroutine spline1d_common_basis2deriv_point
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Computation of the second derivate of the
  !>  cubic spline basis for a knot of the mesh
  !>
  !> \remark Available for periodic or natural function
  !---------------------------------------------------------------------------
  subroutine spline1d_common_basis2deriv_knot( num, nb_knots, h, sbase_second )

    integer                   , intent(in)    :: num
    integer                   , intent(in)    :: nb_knots
    real(F64)                 , intent(in)    :: h
    real(F64), dimension(-1:2), intent(INOUT) :: sbase_second

    real(F64) :: h2

    h2 = h*h
    if (num.lt.nb_knots) then
      sbase_second(-1) = 6._F64/h2
      sbase_second(0)  = -12._F64/h2
      sbase_second(1)  = 6._F64/h2
      sbase_second(2)  = 0._F64
    else
      sbase_second(-1) = 0._F64
      sbase_second(0)  = 6._F64/h2
      sbase_second(1)  = -12._F64/h2
      sbase_second(2)  = 6._F64/h2
    endif

  end subroutine spline1d_common_basis2deriv_knot
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Computation of the coefficients (used for
  !>  the integration calculation) at the boundaries, 
  !>  because the spline -1,0,1,n-1,n,n+1 are not 
  !>  totally in the integration domain
  !>
  !> \remark Available for periodic or natural function
  !---------------------------------------------------------------------------
  subroutine spline1d_common_integration_coef_BC( n, h, indx, factor )
    
    integer                  , intent(in) :: n
    real(F64)                , intent(in) :: h
    integer  , dimension(1:6), intent(out) :: indx
    real(F64), dimension(1:6), intent(out) :: factor     

    !*** array initialisation for the boundary terms ***
    indx(1) = -1   ; factor(1) = h/4._F64
    indx(2) = 0    ; factor(2) = 3._F64*h
    indx(3) = 1    ; factor(3) = (23._F64/4._F64)*h
    indx(4) = n-1  ; factor(4) = (23._F64/4._F64)*h
    indx(5) = n    ; factor(5) = 3._F64*h
    indx(6) = n+1  ; factor(6) = h/4._F64

  end subroutine spline1d_common_integration_coef_BC
  !---------------------------------------------------------------------------

end module spline1d_mod
!---------------------------------------------------------------------------


