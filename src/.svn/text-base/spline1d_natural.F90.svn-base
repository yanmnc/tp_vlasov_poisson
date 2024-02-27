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
!> Cubic spline interpolation for 1D non-periodic
!>  boundary conditions (natural spline)
!>
!> \date 2000-05-15
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module spline1d_natural_class

  use prec_const
  use spline1d_natural_types
  use utils_module, only : locate

  implicit none

  private
  public :: spline1d_natural_new, &
      spline1d_natural_del, & 
      spline1d_natural_BC, &
      spline1d_natural_spline_coef, &
      spline1d_natural_interpol1d, &
      spline1d_natural_deriv1d

contains

  !===========================================================================
  !> 1D spline initialisation for non-periodic function,
  !>  with two possibilities for the boundary conditions : 
  !>   - dirichlet conditions 
  !>   - information on the first derivative
  !> In the case of Dirichlet conditions, an Lagrange 
  !> interpolation is used for the both extremities
  !---------------------------------------------------------------------------
  subroutine spline1d_natural_new( self, n, h )

    use spline1d_mod, only : spline1d_common_integration_coef_BC

    type(spline1d_natural_t), intent(out) :: self
    integer                 , intent(in)  :: n    ! problem dimension
    real(F64)               , intent(in)  :: h    ! discretisation step

    integer :: i, err

    self%n = n
    self%h = h

    !*** memory allocation ***
    allocate(self%xprime(0:n))
    allocate(self%rhs(0:n))
    allocate(self%ipos1d(0:n))
    allocate(self%Adiag(0:n))
    allocate(self%Aodiag(0:n-1))
    allocate(self%Am1gamma(0:n,0:1))
    allocate(self%scoef(-1:n+1))

    !*** ipos1d array initialization ***
    do i = 0,n-1
      self%ipos1d(i) = i
    enddo
    self%ipos1d(n) = n-1

    !*** A factorisation ***
    self%Adiag  = 4._F64
    self%Aodiag = 1._F64
    call DPTTRF( n+1, self%Adiag, self%Aodiag, err )
    if (err.ne.0) then
      write(6,*) 'Natural case : problem in the A factorisation'
      stop
    end if

    !*** Computation of A-1.gamma ***
    self%Am1gamma      = 0._F64
    self%Am1gamma(0,1) = 1._F64 
    self%Am1gamma(n,0) = 1._F64
    call DPTTRS( n+1, 2, self%Adiag, self%Aodiag, self%Am1gamma, n+1, err )
    if (err.ne.0) then
      write(6,*) 'Natural case : problem in A-1.gamma solving'
      stop
    end if

    !**************************************** 
    !*   2x2 matrice (deltab) computation : *
    !*    deltab = delta-lambda.A-1.gamma   *
    !*     where delta = 1  0               *
    !*                   0 -1               *
    !****************************************     
    self%eta1b = 1._F64+self%Am1gamma(n-1,0)
    self%eta2b = self%Am1gamma(n-1,1)
    self%eta3b = -self%Am1gamma(1,0)
    self%eta4b = -(1._F64+self%Am1gamma(1,1))

    ! Computation of deltab determinant
    self%det_deltab = self%eta1b*self%eta4b - &
        self%eta2b*self%eta3b
    if (self%det_deltab.eq.0._F64) then
      write(6,*) &
          'Natural case : problem with determinant equal to'
      write(6,*) self%det_deltab
      stop
    end if

    ! Initialisation of the coefficients required for integration
    call spline1d_common_integration_coef_BC( self%n, self%h, self%indx, &
        self%factor_int )

  end subroutine spline1d_natural_new
  !---------------------------------------------------------------------------


  !===========================================================================
  !> 1D natural spline destructor
  !---------------------------------------------------------------------------
  subroutine spline1d_natural_del( self )

    type(spline1d_natural_t), intent(inout) :: self

    !*** memory allocation ***
    deallocate(self%xprime)
    deallocate(self%rhs)
    deallocate(self%ipos1d)
    deallocate(self%Adiag)
    deallocate(self%Aodiag)
    deallocate(self%Am1gamma)
    deallocate(self%scoef)

  end subroutine spline1d_natural_del
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Compute boundary condition in the case of non-periodic
  !>  function, based on the approximation of the second
  !>  derivative by Lagrange polynomials
  !>   (used for the cubic spline computation)
  !> Input  : BC_left and BC_right 
  !>   for the boundary conditions specifications
  !>   = 0 if Dirichlet condition (func=0)
  !>   = 1 if Neumann (func'=0)
  !>   = 2 if Hermite conditions 
  !>     (func'=approximation of the derivate)
  !> Output : deriv_funcBC(0:1) with
  !>   - deriv_funcBC(0) = func'(n) 
  !>       = 11/6h*func(n) - 
  !>         3/h*func(n-1)+3/2h*func(n-2)-1/3h*func(n-3)  
  !>   - deriv_funcBC(1) = func'(0)
  !>       = -11/6h*func(0) + 
  !>         3/h*func(1)-3/2h*func(2)+1/3h*func(3)
  !---------------------------------------------------------------------------
  subroutine spline1d_natural_BC( func, n, h, BC_left, BC_right, deriv_funcBC )

    real(F64), dimension(0:), intent(in)    :: func
    integer                 , intent(in)    :: n
    real(F64)               , intent(in)    :: h
    integer                 , intent(in)    :: BC_left, BC_right
    real(F64), dimension(0:), intent(inout) :: deriv_funcBC

    !-> computation of the approximate derivative 
    !   at the right hand side
    select case (BC_right)
    case(0)
      deriv_funcBC(0) = (-3._F64*func(n-1) + &
          1.5_F64*func(n-2) - &
          1._F64*func(n-3)/3._F64)/h
    case(1)
      deriv_funcBC(0) = 0._F64
    case(2)
      deriv_funcBC(0) = (11._F64*func(n)/6._F64 - &
          3._F64*func(n-1) + 1.5_F64*func(n-2) - &
          1._F64*func(n-3)/3._F64)/h
    case default
      print*, 'BC_right = ', BC_right, ' must be 0, 1 or 2'
      stop
    end select
    !-> computation of the approximate derivative 
    !   at the left hand side
    select case (BC_left)
    case(0)
      deriv_funcBC(1) = (3._F64*func(1) - 1.5_F64*func(2) + &
          1._F64*func(3)/3._F64)/h
    case(1)
      deriv_funcBC(1) = 0._F64
    case(2)
      deriv_funcBC(1) = (-11._F64*func(0)/6._F64 + &
          3._F64*func(1) - 1.5_F64*func(2) + &
          1._F64*func(3)/3._F64)/h
    case default
      print*, 'BC_left = ', BC_left, ' must be 0, 1 or 2'
      stop
    end select

  end subroutine spline1d_natural_BC
  !---------------------------------------------------------------------------


  !===========================================================================
  !> 1D spline coefficient computation for 
  !>  non-periodic function, i.e solving of
  !>      Atilde| x | = | u |
  !>            | y | = | v |
  !> where x = (c0,...,cn)t
  !>       y = (cn+1,c-1)t
  !>       u = (rhs(0),...,rhs(n))t
  !> and
  !>       v = (rhs'(n),rhs'(0))t
  !>       v = (sigma1,sigma2)t
  !> where sigma1 = h/3*rhs'(n) and
  !>       sigma2 = h/3*rhs'(0)
  !> with the equivalent of a Shur complement method
  !>  Rk1 : rhs'(0) and rhs'(n) are approximated
  !>       with Lagrange polynomials 
  !>  Rk2 : The boundary conditions are given by
  !>    BC_left and BC_right which are equal to
  !>     . 0 if Dirichlet condition (func=0)
  !>     . 1 if Neumann (func'=0)
  !>     . 2 if Hermite conditions 
  !>    (func'=approximation of the derivate)
  !---------------------------------------------------------------------------
  subroutine spline1d_natural_spline_coef( self, rhs, BC_left, BC_right, &
      Sderiv_rhs )

    type(spline1d_natural_t) , intent(inout) :: self
    real(F64), dimension(0:) , intent(in)    :: rhs
    integer                  , intent(in)    :: BC_left
    integer                  , intent(in)    :: BC_right
    ! deriv_rhs(0) = rhs'(n) and deriv_rhs(1) = rhs'(0)
    real(F64)  , optional, &
               dimension(0:) , intent(in)    :: Sderiv_rhs 

    integer                   :: i, n, err
    real(F64)                 :: coeff
    real(F64), dimension(0:1) :: deriv_rhs, rhs2d, yprime 

    n = self%n

    !*** Solving of x'=A-1.u ***
    do i = 0,n
      self%xprime(i) = rhs(i)
    enddo
    !-> if dirichlet condition at r=rmin
    if (BC_left.eq.0) then
      self%xprime(0) = 0._F64
    end if
    !-> if dirichlet condition at r=rmax
    if (BC_right.eq.0) then
      self%xprime(n) = 0._F64
    end if
    call DPTTRS( n+1, 1, self%Adiag, self%Aodiag, self%xprime, n+1, err )
    if (err.ne.0) then
      write(6,*) 'Natural case : problem in x''=A-1.u solving'
      stop
    end if

    !*** computation of rhs'(0) and rhs'(n)  ***
    !-> deriv_rhs(0) = rhs'(n) and deriv_rhs(1) = rhs'(0)
    if (.not.present(Sderiv_rhs)) then
      call spline1d_natural_BC( rhs, self%n, self%h, BC_left, BC_right, &
          deriv_rhs )
    else
      deriv_rhs(0) = Sderiv_rhs(0)
      deriv_rhs(1) = Sderiv_rhs(1)
    end if

    !*** v-lamda.A-1u assembling ***
    coeff    = self%h/3._F64
    rhs2d(0) = coeff*deriv_rhs(0)+self%xprime(n-1)
    rhs2d(1) = coeff*deriv_rhs(1)-self%xprime(1)

    !*** Solving of the 2X2 system :      ***
    !***   deltab.yprime = v-lamda.A-1u   ***
    yprime(0) = (1._F64/self%det_deltab)* &
        (self%eta4b*rhs2d(0)-self%eta2b*rhs2d(1))
    yprime(1) = (1._F64/self%det_deltab)* &
        (-self%eta3b*rhs2d(0)+self%eta1b*rhs2d(1))

    !*** Computation of x = x'-A-1.gamma.y ***
    do i = 0,n
      self%scoef(i) = self%xprime(i) - &
          self%Am1gamma(i,0)*yprime(0) &
          -self%Am1gamma(i,1)*yprime(1)
    enddo
    self%scoef(-1)  = yprime(1)
    self%scoef(n+1) = yprime(0)

  end subroutine spline1d_natural_spline_coef
  !---------------------------------------------------------------------------


  !===========================================================================
  !>  natural function interpolation by using the 
  !>  cubic spline coefficients calculated 
  !---------------------------------------------------------------------------
  subroutine spline1d_natural_interpol1d( self, x, fgrid, BC_left, BC_right, &
      nbstar, xstar, finterpol )

    use spline1d_mod, only : spline1d_common_basis

    type(spline1d_natural_t)  , intent(inout) :: self
    real(F64)  , dimension(0:), intent(in)    :: x, fgrid
    integer                   , intent(in)    :: BC_left
    integer                   , intent(in)    :: BC_right
    integer                   , intent(in)    :: nbstar
    real(F64)  , dimension(0:), intent(in)    :: xstar
    real(F64)  , dimension(0:), intent(inout) :: finterpol

    integer                      :: i, ipos, k    
    real(F64)  , dimension(-1:2) :: sbase
    character(LEN=50), parameter :: subp = " spline1d_natural_interpol1d "//char(0)

    !*** spline coefficient computation *** 
    call spline1d_natural_spline_coef( self, fgrid, BC_left, BC_right )

    do i = 0, nbstar
      !*** array position location ***
      call locate( xstar(i), x, self%n, self%h, subp, ipos )

      !*** Calculation of cubic spline basis ***
      call spline1d_common_basis( x(ipos), xstar(i), x(ipos+1), self%h, sbase )

      !*** Computes f(x*) ***
      finterpol(i) = 0._F64
      do k=-1,2     
        finterpol(i) = finterpol(i)+self%scoef(ipos+k)*sbase(k)
      end do
    enddo

  end subroutine spline1d_natural_interpol1d
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Compute the first derivative of a non-periodic function
  !---------------------------------------------------------------------------
  subroutine spline1d_natural_deriv1d( self, n, h, BC_left, BC_right, func, &
      dfuncdx )

    use spline1d_mod, only : spline1d_common_basisderiv

    type(spline1d_natural_t) , intent(inout) :: self
    integer                  , intent(in)    :: n
    real(F64)                , intent(in)    :: h
    integer                  , intent(in)    :: BC_left
    integer                  , intent(in)    :: BC_right
    real(F64), dimension(0:n), intent(in)    :: func
    real(F64), dimension(0:n), intent(inout) :: dfuncdx

    integer                    :: ix, i
    real(F64), dimension(-1:2) :: sbase_prime
    real(F64)                  :: dfuncdx_tmp

    !*** computation of the spline coefficients of the function ***
    call spline1d_natural_spline_coef( self, func, BC_left, BC_right )

    !*** computation of the first derivative of the function ***
    do ix = 0,n
      call spline1d_common_basisderiv( ix, n, h, sbase_prime )
      dfuncdx_tmp = 0._F64      
      do i = -1,2
        dfuncdx_tmp = dfuncdx_tmp + & 
            sbase_prime(i) * self%scoef(self%ipos1d(ix)+i) 
      enddo
      dfuncdx(ix) = dfuncdx_tmp 
    enddo

  end subroutine spline1d_natural_deriv1d
  !---------------------------------------------------------------------------

end module spline1d_natural_class
!---------------------------------------------------------------------------


