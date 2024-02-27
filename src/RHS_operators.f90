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
! file : krook_operator.f90
! 
!  Krook operator
!-------------------------------------------------
module RHS_module
  use fdistribu2d_class
  implicit none
      
  real(RKIND), dimension(:), pointer, private :: AA, BB, CC

  !******************************
  contains
  !******************************

  !----------------------------------------------
  ! Initialisation of the diffusion coefficients
  !----------------------------------------------
  subroutine new_RHS
    use globals, only : Nx
    !--> for diffusion definition
    allocate(AA(0:Nx-1))
    allocate(BB(0:Nx-1))
    allocate(CC(0:Nx-1))
  end subroutine new_RHS


  !----------------------------------------------
  ! Delete the diffusion coefficients
  !----------------------------------------------
  subroutine del_RHS()
    deallocate(AA)
    deallocate(BB)
    deallocate(CC)
  end subroutine del_RHS


  !------------------------------------------------------------
  ! Treatment of the Krook operator
  !   df/dt = - nu (f-feq)
  ! The analytical solution of this equation is:
  !  f(t+dt) = exp(-nu*dt)f(t)+feq(1-exp(-nu*dt))
  !------------------------------------------------------------  
  subroutine solve_Krook_operator(geom,f,nu,dt)
  use globals, only : dens, velo, temp, OhmsLaw
  use physics_module, only : compute_fluid_moments
    type(geometry)   , intent(in)    :: geom
    type(fdistribu2d), intent(inout) :: f
    real(RKIND)      , intent(in)    :: nu
    real(RKIND)      , intent(in)    :: dt

    integer     :: ix, iv
    real(RKIND) :: fM, vloc, inv_sqrt_TWOPI, coef1, coef2, &
      dens_loc, temp_loc, fM0, coef_OhmsLaw

    coef_OhmsLaw = 0._RKIND
    if (OhmsLaw) then
      coef_OhmsLaw = 1._RKIND
    end if

    inv_sqrt_TWOPI = 1._RKIND/sqrt(2._RKIND*pi)
    coef1 = (0.5_RKIND+0.5_RKIND*coef_OhmsLaw)*nu*dt
    coef2 = 1._RKIND/(1._RKIND+coef1)

    ! Construct the Maxwellian
    call compute_fluid_moments(geom,f,1)

    do iv = 0,f%n2
      vloc = geom%vgrid(iv)
      do ix = 0,f%n1-1
        dens_loc = dens(ix)
        temp_loc = temp(ix)
        fM  = dens_loc*inv_sqrt_TWOPI/sqrt(temp_loc)* &
             exp(-(vloc-velo(ix))*(vloc-velo(ix))/(2._RKIND*temp_loc))
        fM0 = dens_loc*inv_sqrt_TWOPI/sqrt(temp_loc)* &
             exp(-vloc*vloc/(2._RKIND*temp_loc))
        f%values(ix,iv) = coef2*( &
          f%values(ix,iv)*(1._RKIND-coef1) + &
          (2._RKIND-coef_OhmsLaw)*coef1*(fM+coef_OhmsLaw*fM0))
      end do
    end do
  end subroutine solve_Krook_operator


  !------------------------------------------------------------
  ! initialization of diffusion matrix 
  !------------------------------------------------------------  
  subroutine init_diffusion(geom,dt)
    use geometry_class
    use globals, only : diff_coeff
    type(geometry), intent(in) :: geom     

    real(RKIND),    intent(in) :: dt  
    real(RKIND)                :: alpha
    integer                    :: ix

    alpha = diff_coeff*dt/(2._RKIND*geom%dx*geom%dx)
    
    do ix=0,geom%Nx-1
      AA(ix) = -alpha
      BB(ix) = 1._RKIND+2._RKIND*alpha
      CC(ix) = -alpha
    end do
  end subroutine init_diffusion
    

  !------------------------------------------------------------
  ! Treatment of the Diffusion operator in configuration space
  !   df/dt = D d^2(f)/dx^2
  ! Numerical method: Crank-Nicolson
  !------------------------------------------------------------  
  subroutine solve_Diffusion_operator(f)
    use utils_module
    use globals, only : Nx
    type(fdistribu2d), intent(inout)   :: f

    real(RKIND), dimension(0:Nx-1)     :: floc, RR
    integer                            :: ix, iv
    
    do iv = 0,f%n2
      floc(0:Nx-1) = f%values(0:Nx-1,iv)
      do ix = 1,f%n1-1
        RR(ix) = -AA(ix)*floc(ix-1)+(2._RKIND-BB(ix))*floc(ix)-&
          CC(ix)*floc(ix+1)
      end do
      RR(0) = -AA(0)*floc(f%n1-1)+(2._RKIND-BB(0))*floc(0)-&
        CC(0)*floc(1)
      call cyclic(AA,BB,CC,AA(0),AA(0),RR,floc,f%n1-1)
      do ix = 0,f%n1-1
        f%values(ix,iv) = floc(ix)
      end do
      f%values(f%n1,iv) = floc(0)
    end do
  end subroutine solve_Diffusion_operator
end module RHS_module
