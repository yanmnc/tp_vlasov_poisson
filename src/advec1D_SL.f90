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

!---------------------------------------------------------
! file : advec1D_SL.f90
!
!   1D advection in x and v directions used in
!   Semi-Lagrangian method
!---------------------------------------------------------
module advec1D_SL_module
  use globals, only : Nx, Nv, OhmsLaw, E0
  use geometry_class
  use fdistribu2d_class
  use interpolation_module
  
  implicit none
      
  !******************************
  contains
  !******************************  

  !-----------------------------------------------------
  ! splitting in x direction 
  !  - Resolution of : df/dt+v*df/dx=0
  !     . dx = deltat*v
  !     . f**(x,v) = f*(x-dx,v)
  !-----------------------------------------------------   
  subroutine advec1D_x(geom,dt,f)
    type(geometry)   , intent(in)    :: geom
    real(RKIND)      , intent(in)    :: dt
    type(fdistribu2d), intent(inout) :: f

    integer     :: ix, iv
    real(RKIND) :: x0, xn, xstar, dx
    real(RKIND) :: finterpol
    logical     :: bound 
    
    real(RKIND), dimension(0:Nx)    :: rhs
    real(RKIND), dimension(-1:Nx+1) :: scoef_x

    x0 = geom%xgrid(0)
    xn = geom%xgrid(Nx)
    do iv = 0,Nv
      call boundary(geom,iv,bound)
      if (.not.bound) then     
        do ix = 0,Nx
          rhs(ix) = f%values(ix,iv)
        end do
        !*** cubic spline computation ***
        call compute_spline_x(f,rhs,scoef_x)
        !*** splitting in x direction ***
        do ix = 0,Nx-1
          dx    = dt*geom%vgrid(iv)
          xstar = geom%xgrid(ix) - dx
          if (xstar.gt.xn) then
            xstar = xstar-abs(floor(xstar/(xn-x0)))*abs(xn-x0)  
          elseif (xstar.lt.x0) then
            xstar = xstar+abs(floor(xstar/(xn-x0)))*abs(xn-x0)
          endif
          !*** f interpolation ***
          call interpol1d_x(geom,scoef_x,geom%xgrid,xstar, &
              finterpol)
          f%values(ix,iv) = finterpol
        end do
        f%values(Nx,iv) = f%values(0,iv)
      end if
    end do
  end subroutine advec1D_x


  !-----------------------------------------------------
  ! splitting in v parallel direction
  !  - Resolution of : df/dt+(dv/dt)*df/dv=0
  !     . dpar = deltat*E(x,tn)
  !     . f*(x,v) = fn(x,v-dpar)
  !-----------------------------------------------------   
  subroutine advec1D_v(geom,dt,E,f)
    type(geometry)                  , intent(in)    :: geom
    real(RKIND)                     , intent(in)    :: dt
    real(RKIND)      , dimension(0:), intent(in)    :: E
    type(fdistribu2d)               , intent(inout) :: f

    integer     :: ix, iv
    real(RKIND) :: vstar, dv
    real(RKIND) :: finterpol, coef_OhmsLaw
    logical     :: bound 
    
    real(RKIND), dimension(0:Nv)    :: rhs
    real(RKIND), dimension(0:1)     :: deriv_rhs
    real(RKIND), dimension(-1:Nv+1) :: scoef_v

    coef_OhmsLaw = 0._RKIND
    if (OhmsLaw) then
      coef_OhmsLaw = 1._RKIND
    end if

    do ix = 0,Nx-1
      do iv = 0,Nv
        rhs(iv) = f%values(ix,iv)
      end do
      !*** cubic spline computation ***
      deriv_rhs(0) = (11._RKIND*rhs(f%n2)/6._RKIND - &
        3._RKIND*rhs(f%n2-1) +  &
        1.5_RKIND*rhs(f%n2-2) - &
        1._RKIND*rhs(f%n2-3)/3._RKIND)/f%h2
      deriv_rhs(1) = (-11._RKIND*rhs(0)/6._RKIND + &
        3._RKIND*rhs(1)  - &
        1.5_RKIND*rhs(2) + &
        1._RKIND*rhs(3)/3._RKIND)/f%h2
      call compute_spline_v(f,rhs,deriv_rhs,scoef_v)
      do iv = 0,Nv
        call boundary(geom,iv,bound)
        if (.not.bound) then        
          !*** splitting in v direction ***
          dv  = dt*(E(ix)+coef_OhmsLaw*E0)
          vstar = geom%vgrid(iv) - dv       
          vstar = min(max(geom%vgrid(0),vstar), &
            geom%vgrid(Nv))   
          !*** f interpolation ***
          call interpol1d_v(geom,scoef_v, &
              geom%vgrid,vstar,finterpol)
          f%values(ix,iv) = finterpol
        end if
      end do
    end do
    f%values(Nx,0:Nv) = f%values(0,0:Nv)
  end subroutine advec1D_v
end module advec1D_SL_module










