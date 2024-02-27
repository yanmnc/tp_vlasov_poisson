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
!> Cubic spline interpolation for 1D periodic
!>  boundary conditions (natural spline)
!>
!> \date 2000-05-15
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module spline1d_periodic_class

  use prec_const
  use spline1d_periodic_types
  use utils_module, only : locate

  implicit none

  private
  public :: &
      spline1d_periodic_new, &
      spline1d_periodic_del, &
      spline1d_periodic_deriv1d, &
      spline1d_periodic_spline_coef

contains

  !===========================================================================
  !> 1D spline initialisation for periodic function
  !---------------------------------------------------------------------------
  subroutine spline1d_periodic_new( self, n, h )

    use spline1d_mod, only : spline1d_common_integration_coef_BC

    type(spline1d_periodic_t), intent(out) :: self
    real(F64)                , intent(in)  :: h
    integer                  , intent(in)  :: n

    integer   :: i, err
    real(F64) :: coeff1, coeff2

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

    !*** A factorisation ***
    self%Adiag  = 4._F64
    self%Aodiag = 1._F64
    call DPTTRF( n+1, self%Adiag, self%Aodiag, err )
    if (err.ne.0) then
      write(6,*) 'Periodic case : problem in the A factorisation'
      stop
    end if

    !*** ipos1d array initialization ***
    do i = 0,n-1
      self%ipos1d(i) = i
    enddo
    self%ipos1d(n) = n-1

    !*** Computation of A-1.gamma ***
    self%Am1gamma      = 0._F64
    self%Am1gamma(0,1) = 1._F64 
    self%Am1gamma(n,0) = 1._F64
    call DPTTRS( n+1, 2, self%Adiag, self%Aodiag, self%Am1gamma, n+1, err )
    if (err.ne.0) then
      write(6,*) 'Periodic case : problem in A-1.gamma solving'
      stop
    end if

    !**************************************** 
    !*   2x2 matrice (deltab) computation : *
    !*    deltab = delta-lambda.A-1.gamma   *
    !*     where delta = 4 0                *
    !*                   1 0                *
    !****************************************     
    coeff1 = 3._F64/h
    coeff2 = 6._F64/(h**2)
    self%eta1b = coeff1*(-1._F64-self%Am1gamma(1,0)- &
        self%Am1gamma(n-1,0))
    self%eta2b = coeff1*(-1._F64-self%Am1gamma(1,1)- &
        self%Am1gamma(n-1,1))
    self%eta3b = coeff2 * &
        (-1._F64+2._F64*self%Am1gamma(0,0) &
        -self%Am1gamma(1,0)+self%Am1gamma(n-1,0)- &
        2._F64*self%Am1gamma(n,0))
    self%eta4b = coeff2 * &
        (1._F64+2._F64*self%Am1gamma(0,1) &
        -self%Am1gamma(1,1)+self%Am1gamma(n-1,1)- &
        2._F64*self%Am1gamma(n,1))

    ! Computation of deltab determinant
    self%det_deltab = self%eta1b*self%eta4b - &
        self%eta2b*self%eta3b
    if (self%det_deltab.eq.0._F64) then
      write(6,*) &
          'Periodic case : problem with determinant equal to'
      write(6,*) self%det_deltab
      stop
    end if

    ! Initialisation of the coefficients required for integration
    call spline1d_common_integration_coef_BC( self%n, self%h, self%indx, &
        self%factor_int )

  end subroutine spline1d_periodic_new
  !---------------------------------------------------------------------------


  !===========================================================================
  !> 1D periodic spline destructor
  !---------------------------------------------------------------------------
  subroutine spline1d_periodic_del( self )

    type(spline1d_periodic_t), intent(inout) :: self

    !*** memory allocation ***
    deallocate(self%xprime)
    deallocate(self%rhs)
    deallocate(self%ipos1d)
    deallocate(self%Adiag)
    deallocate(self%Aodiag)
    deallocate(self%Am1gamma)
    deallocate(self%scoef)

  end subroutine spline1d_periodic_del
  !---------------------------------------------------------------------------


  !===========================================================================
  !> 1D spline coefficient computation for 
  !>  periodic function, i.e solving of
  !>      Atilde| x | = | u |
  !>            | y | = | v |
  !> where x = (c0,...,cn-2)t
  !>       y = cn+1
  !>       u = (rhs(0),...,rhs(n-2))t
  !>       v = rhs(n-1)
  !> with the equivalent of a Shur complement method
  !---------------------------------------------------------------------------
  subroutine spline1d_periodic_spline_coef( self, rhs )

    type(spline1d_periodic_t) , intent(inout) :: self
    real(F64)  , dimension(0:), intent(in)    :: rhs

    integer                   :: i, n, err
    real(F64), dimension(0:1) :: rhs2d, yprime  
    real(F64)                 :: coeff1, coeff2    

    n = self%n    
    !*** Solving of x'=A-1.u ***
    do i = 0,n
      self%xprime(i) = rhs(i)
    enddo
    call DPTTRS( n+1, 1, self%Adiag, self%Aodiag, self%xprime, n+1, err )
    if (err.ne.0) then
      write(6,*) 'Periodic case : problem in x''=A-1.u solving'
      stop
    end if

    !*** v-lamda.A-1u assembling ***
    coeff1   = 3._F64/self%h
    coeff2   = 6._F64/(self%h**2)
    rhs2d(0) = -coeff1*(self%xprime(1)+self%xprime(n-1))
    rhs2d(1) = -coeff2 * &
        (-2._F64*self%xprime(0) + &
        self%xprime(1)-self%xprime(n-1)+ &
        2._F64*self%xprime(n))

    !*** Solving of the 2X2 system :     ***
    !***    deltab.yprime = v-lamda.A-1u ***
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
    self%scoef(n+1) = yprime(0)
    self%scoef(-1)  = yprime(1)

  end subroutine spline1d_periodic_spline_coef
  !---------------------------------------------------------------------------


  !===========================================================================
  !>  periodic function interpolation by using the 
  !>  cubic spline coefficients calculated 
  !---------------------------------------------------------------------------
  subroutine spline1d_periodic_interpol1d( self, x, fgrid, nbstar, xstar, &
      finterpol )

    use spline1d_mod, only : spline1d_common_basis

    type(spline1d_periodic_t)     , intent(inout) :: self
    real(F64), dimension(0:self%n), intent(in)    :: x, fgrid
    integer                       , intent(in)    :: nbstar
    real(F64), dimension(0:nbstar), intent(in)    :: xstar
    real(F64), dimension(0:nbstar), intent(inout) :: finterpol

    integer                      :: i, ipos, k    
    real(F64)  , dimension(-1:2) :: sbase
    character(LEN=50), parameter :: subp = " spline1d_periodic_interpol1d "//char(0)

    !*** spline coefficient computation *** 
    call spline1d_periodic_spline_coef( self, fgrid )

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

  end subroutine spline1d_periodic_interpol1d
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Compute the first derivative of a periodic function
  !---------------------------------------------------------------------------
  subroutine spline1d_periodic_deriv1d( self, n, h, func, dfuncdx )

    use spline1d_mod, only : spline1d_common_basisderiv

    type(spline1d_periodic_t), intent(inout) :: self
    integer                  , intent(in)    :: n
    real(F64)                , intent(in)    :: h
    real(F64), dimension(0:n), intent(in)    :: func
    real(F64), dimension(0:n), intent(inout) :: dfuncdx

    integer                    :: ix, i 
    real(F64), dimension(-1:2) :: sbase_prime
    real(F64), dimension(0:1)  :: deriv
    real(F64)                  :: dfuncdx_tmp

    !*** computation of the spline coefficients of the function ***
    self%rhs(0:n) = func(0:n)
    call spline1d_periodic_spline_coef( self, self%rhs )

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

  end subroutine spline1d_periodic_deriv1d
  !---------------------------------------------------------------------------


end module spline1d_periodic_class
!---------------------------------------------------------------------------
