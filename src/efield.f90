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
! file : efield.f90
!
!  . computation of the electric potential 
!   -> resolution of the Poisson equation
!      -d^2Phi\dx^2 = \int f dv - 1
!  . computation of the field
!      E = -dPhi/dx
!---------------------------------------------
module efield_module
  use prec_const
  use globals, only : Nx
  use utils_module
  use geometry_class
  use fdistribu2d_class
  
  implicit none

  !******************************
  contains
  !******************************  

  !-----------------------------------------------------
  ! used for first electric field allocation
  !-----------------------------------------------------
  subroutine new_Efield(Phi,Ex)
    use globals, only : Nx
    real(RKIND), dimension(:), pointer :: Phi
    real(RKIND), dimension(:), pointer :: Ex

    !*** electric potential allocation ***
    allocate(Phi(0:Nx))

    !*** electric field allocation ***
    allocate(Ex(0:Nx))
  end subroutine new_Efield


  !-----------------------------------------------------
  ! used for first electric field allocation
  !-----------------------------------------------------
  subroutine del_Efield(Phi,Ex)
    real(RKIND), dimension(:), pointer :: Phi
    real(RKIND), dimension(:), pointer :: Ex

    !*** electric potential deallocation ***
    deallocate(Phi)
    !*** electric field deallocation ***
    deallocate(Ex)
  end subroutine del_Efield


  !-------------------------------------------------
  ! Computation of the RHS
  !  rho = \int f dv - 1 
  ! of the quasi-neutrality equation 
  !-------------------------------------------------
  subroutine compute_rho(geom,f,rho)
    use geometry_class
    type(geometry)            , intent(in)  :: geom
    type(fdistribu2d)         , intent(in)  :: f
    real(RKIND), dimension(0:), intent(out) :: rho

    integer     :: ix, iv
    real(RKIND) :: f_tmp
    real(RKIND) :: ni
    
    !*** compute ni(x) = \int f dv  ***
    do ix = 0,f%n1-1
      ni = 0._RKIND
      do iv = 0,f%n2
        f_tmp = f%values(ix,iv)
        ni    = ni + f_tmp*geom%coeff_intdv(iv) 
      end do
      rho(ix) = ni - 1._RKIND
    end do
    rho(f%n1) = rho(0)
  end subroutine compute_rho


  !-------------------------------------------------
  ! Computation of the electric potential by solving
  !  -d^2Phi\dx^2 = rho with rho = \int f dv - 1
  ! in the Fourier space, i.e. :
  !   kx^2 Phi_kx = rho_kx 
  !-------------------------------------------------  
  subroutine compute_Phi(geom,iter,f,Phi)
    use globals, only : deltat, mode_ext, omega_ext, Phi_ext0
    use fftNRF90_module
    type(geometry)            , intent(in)    :: geom
    integer                   , intent(in)    :: iter
    type(fdistribu2d)         , intent(inout) :: f
    real(RKIND), dimension(0:), intent(out)   :: Phi
    
    complex(CKIND), dimension(:), pointer :: rhsfft_c
    real(RKIND)   , dimension(:), pointer :: kvector 

    integer     :: i, nx, nxm1
    real(RKIND) :: wkx, kx_ext, Phi_ext, time


    Phi = 0._RKIND

    !*** initialize local variables ***
    nx   = size(geom%xgrid)
    nxm1 = nx - 1 

    !*** tempory array allocation ***
    allocate(rhsfft_c(1:nxm1))
    allocate(kvector(1:nxm1))

    !*** compute the RHS of Poisson equation ***
    ! -> rho is put in  Phi
    call compute_rho(geom,f,Phi)
    do i = 1,nxm1
      rhsfft_c(i) = Phi(i-1)
    end do

    !*** perform FFT 1D in x  ***
    call fft(rhsfft_c)

    !*** compute kvector for the FFT in x direction ***
    wkx = 2._RKIND*pi/(nxm1*geom%dx)
    do i = 1,(nxm1/2)+1
      kvector(i) = (i-1) * wkx
    end do
    do i = (nxm1/2)+2,nxm1
      kvector(i) = -1._RKIND * (nxm1-(i-1)) * wkx
    end do

! ATTENTION A VERIFIER
    rhsfft_c(1) = 0._CKIND
    do i = 2,nxm1
      rhsfft_c(i) = rhsfft_c(i)/(kvector(i)**2)
    end do

    !*** Perform inverse FFT 1D  of the system solution ***
    call fft_inv(rhsfft_c)

    !*** Duplicate the result in the Phi 1D array ***
    do i = 0,nxm1-1
      Phi(i) = real(rhsfft_c(i+1), kind(Phi(i)))
    end do

    !*** Addition of an external potential ***
    kx_ext = 2._RKIND*pi*mode_ext/(nxm1*geom%dx)
    time   = (iter-1)*deltat
    do i = 0,nxm1-1
      Phi_ext = Phi_ext0*cos(kx_ext*geom%xgrid(i)-omega_ext*time)
      Phi(i)  = Phi(i) + Phi_ext
    end do
    Phi(nxm1) = Phi(0)
  end subroutine compute_Phi

  
  !-----------------------------------------------------
  !  computation of the electric field coordinates in
  !     cartesian axis : 
  !    . Ex = -dPhi/dx
  !-----------------------------------------------------   
  subroutine compute_Ex(geom,Phi,Ex)
    type(geometry)               , intent(in)    :: geom
    real(RKIND)   , dimension(0:), intent(in)    :: Phi
    real(RKIND)   , dimension(0:), intent(inout) :: Ex
    
    integer :: ix

    ! electric potential derivates
    real(RKIND), dimension(:), pointer :: dPhidx
    
    !*** dPhi/dx computation ***
    allocate(dPhidx(0:geom%Nx))
    call deriv1(Phi,dPhidx,geom%Nx,geom%dx,1)
  
    !*** Ex and Ey computation ***
    do ix = 0,geom%Nx
      Ex(ix) = -dPhidx(ix)
    end do

    !*** deallocation ***
    deallocate(dPhidx) 
  end subroutine compute_Ex


  !-----------------------------------------------------
  ! used for global electric field computation
  !-----------------------------------------------------
  subroutine compute_Efield(geom,iter,f,Phi,Ex)
    type(geometry)            , intent(in)    :: geom
    integer                   , intent(in)    :: iter
    type(fdistribu2d)         , intent(inout) :: f
    real(RKIND), dimension(0:), intent(inout) :: Phi
    real(RKIND), dimension(0:), intent(inout) :: Ex

    call compute_Phi(geom,iter,f,Phi)
    call compute_Ex(geom,Phi,Ex)
  end subroutine compute_Efield

end module efield_module



