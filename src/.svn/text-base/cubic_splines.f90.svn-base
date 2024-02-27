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
! file : cubic_splines
!
! cubic spline interpolation by using local
!  cubic spline bases
! => Developped by G. LATU and N. Crouseilles
!---------------------------------------------
module cubic_splines_module
  use prec_const
  implicit none

  public :: new_splinehh, del_splinehh 

  type, public :: splinehh
     ! number of points in x and y directions
     integer                              :: nbx, nby
     ! buffer size
     integer                              :: buf
     ! information on mesh in x
     real(RKIND)                          :: x_0, dx, invhx
     ! information on mesh in y 
     real(RKIND)                          :: y_0, dy, invhy
     ! grid in x and y directions
     real(RKIND), dimension(:)  , pointer :: xgrid, ygrid
     ! diagonal terms of L in LU decomposition
     real(RKIND), dimension(:)  , pointer :: lmatx, lmaty
     ! diagonal terms of U in LU decomposition
     real(RKIND), dimension(:)  , pointer :: umatx, umaty
     ! coefficients des splines
     real(RKIND), dimension(:,:), pointer :: coef, bcoef
     real(RKIND), dimension(:,:), pointer :: aux
     real(RKIND), dimension(:,:), pointer :: auy
  end type splinehh

  !*** public variables ***
  real(RKIND) rightconst, leftconst
  
  integer :: bufsize = 1
  integer :: stencil = 10

  !******************************
  contains
  !******************************

  !-------------------------------------------------------- 
  ! Computes the derivative of 10th order
  !--------------------------------------------------------
  function derive10(vect,dh)
    real(RKIND), dimension (-10:10), intent(in) :: vect
    real(RKIND)                    , intent(in) :: dh
    real(RKIND)                                 :: derive10

    real(RKIND) :: tpr

    tpr = &
      .2214309755e-5*vect(-10)&
      -.1771447804e-4*vect(-9)&
      +.7971515119e-4*vect(-8)&
      -.3011461267e-3*vect(-7)&
      +.1113797807e-2*vect(-6)&
      -.4145187862e-2*vect(-5)&
      +.1546473933e-1*vect(-4)&
      -.5771376946e-1*vect(-3)&
      +.2153903385*vect(-2)&
      -.8038475846*vect(-1)
    tpr = tpr &
      -.2214309755e-5*vect(10)&
      +.1771447804e-4*vect(9)&
      -.7971515119e-4*vect(8)&
      +.3011461267e-3*vect(7)&
      -.1113797807e-2*vect(6)&
      +.4145187862e-2*vect(5)&
      -.1546473933e-1*vect(4)&
      +.5771376946e-1*vect(3)&
      -.2153903385*vect(2) &
      +.8038475846*vect(1)
    derive10 = tpr/dh
  end function derive10
  

  !-------------------------------------------------------- 
  ! Computes the right contribution of the derivative of
  !  10th order
  !--------------------------------------------------------
  function rightderive10(vect,dh)
    real(RKIND), dimension(1:10), intent(in) :: vect
    real(RKIND)                 , intent(in) :: dh
    real(RKIND)                              :: rightderive10

    rightderive10 = &
      -.2214309755e-5*vect(10)&
      +.1771447804e-4*vect(9)&
      -.7971515119e-4*vect(8)&
      +.3011461267e-3*vect(7)&
      -.1113797807e-2*vect(6)&
      +.4145187862e-2*vect(5)&
      -.1546473933e-1*vect(4)&
      +.5771376946e-1*vect(3)&
      -.2153903385*vect(2) &
      +.8038475846*vect(1)
    rightderive10 = rightderive10/dh
  end function rightderive10


  !-------------------------------------------------------- 
  ! Computes the left contribution of the derivative of
  !  10th order
  !--------------------------------------------------------
  function leftderive10(vect,dh)
    real(RKIND), dimension(-10:-1), intent(in) :: vect
    real(RKIND)                   , intent(in) :: dh
    real(RKIND)                                :: leftderive10

    leftderive10 = &
      .2214309755e-5*vect(-10)&
      -.1771447804e-4*vect(-9)&
      +.7971515119e-4*vect(-8)&
      -.3011461267e-3*vect(-7)&
      +.1113797807e-2*vect(-6)&
      -.4145187862e-2*vect(-5)&
      +.1546473933e-1*vect(-4)&
      -.5771376946e-1*vect(-3)&
      +.2153903385*vect(-2)&
      -.8038475846*vect(-1)
    leftderive10 = leftderive10/dh
  end function leftderive10


  !------------------------------------------------------ 
  ! Constructor
  !------------------------------------------------------
  subroutine new_splinehh(hsthis,nbx,nby,x_0,y_0,dx,dy)
    use globals
    type(splinehh), intent(inout) :: hsthis
    ! number of points in x and y directions
    integer       , intent(in)    :: nbx, nby 
    ! coordinates of mesh origin in each direction
    real(RKIND)   , intent(in)    :: x_0, y_0
    ! cell size in each direction 
    real(RKIND)   , intent(in)    :: dx, dy

    !*** local variables ***
    real(RKIND) :: vect(-stencil:stencil) 
    integer     :: i, j, info

    ! initialization of the mesh
    hsthis%nbx = nbx
    hsthis%nby = nby
    hsthis%x_0 = x_0
    hsthis%dx = dx
    hsthis%y_0 = y_0
    hsthis%dy = dy
    allocate(hsthis%xgrid(1:nbx+3))
    do i = 1,nbx+3
      hsthis%xgrid(i)= x_0+(i-1)*dx
    end do
    allocate(hsthis%ygrid(1:nby+3))
    do i = 1,nby+3
      hsthis%ygrid(i)= y_0+(i-1)*dy
    end do

    ! memory allocation of matrices
    allocate(hsthis%coef(-2:nbx+1,-2:nby+1))
    allocate(hsthis%bcoef(-2:nbx+1,-2:nby+1))
    allocate(hsthis%aux(-2:nbx+1,-2:nby+1))
    allocate(hsthis%auy(-2:nbx+1,-2:nby+1))
    allocate(hsthis%lmatx(1:nbx+2))
    allocate(hsthis%lmaty(1:nby+2))
    allocate(hsthis%umatx(1:3))
    allocate(hsthis%umaty(1:3))

    ! Compute LU decomposition in x and y directions
    call LUHermiteSpline(hsthis%lmatx,hsthis%umatx,nbx)
    call LUHermiteSpline(hsthis%lmaty,hsthis%umaty,nby)
      
    hsthis%invhx = 1._RKIND/dx
    hsthis%invhy = 1._RKIND/dy
    vect(-stencil:stencil) = 1._RKIND
      
    rightconst = rightderive10(vect(1:stencil),1.0_RKIND)
    leftconst  = leftderive10(vect(-stencil:-1),1.0_RKIND)
  end subroutine new_splinehh


  !------------------------------------------------------ 
  ! Destructor
  !------------------------------------------------------
  subroutine del_splinehh(hsthis)
    type(splinehh), intent(inout) :: hsthis                

    deallocate(hsthis%aux)
    deallocate(hsthis%auy)
    deallocate(hsthis%xgrid)
    deallocate(hsthis%ygrid)
    deallocate(hsthis%lmatx)
    deallocate(hsthis%lmaty)
    deallocate(hsthis%umatx)
    deallocate(hsthis%umaty)
    deallocate(hsthis%bcoef)
    deallocate(hsthis%coef)
  end subroutine del_splinehh


  !*****************************************************
  ! Computation of the spline coefficients
  !*****************************************************
  !----------------------------------------------------------
  ! Compute LU decomposition for splines with 
  !  Hermite boundary conditions :
  !   -> initialization of li and di
  !      which are the inverse of diagonal terms of U matrix 
  !      in the LU decomposition
  !----------------------------------------------------------  
  subroutine LUHermiteSpline(li,di,n)
    real(RKIND), dimension(1:), intent(out) :: li,di 
    integer                   , intent(in)  :: n

    !*** local variables ***
    integer                     :: i
    real(RKIND)                 :: d
    real(RKIND), dimension(1:n) :: Dinv
   
    di    = 0._RKIND
    li    = 0._RKIND
    d     = 7._RKIND/2._RKIND
    li(1) = 0.25_RKIND
    li(2) = 2._RKIND/7._RKIND
    do i = 2,n-2
      li(i) = 1._RKIND/d
      d     = 4._RKIND-li(i)
    enddo
    di(1)   = d
    li(n-1) = 1._RKIND/di(1)
    di(2)   = 4._RKIND-li(n-1)
    li(n)   = 1._RKIND/(di(1)*di(2)) 
    di(3)   = 1._RKIND-li(n)
    li(n+1) = 1._RKIND/di(3)    
  end subroutine LUHermiteSpline


  !------------------------------------------------------------- 
  ! Compute Hermite spline coefficients
  !   rhs(:,1:ntx) contains values of function to be 
  !   interpolated on the full patch and rhs(:0) and 
  !   rhs(:,nbx+1) contains the values of the
  !   derivatives of the function at points x_0 and x_N 
  !   respectively.
  !   The transposed function is being used for better 
  !   data locality
  !------------------------------------------------------------- 
  subroutine hermite(hsthis,rhs)
    type(splinehh)                  , intent(inout) :: hsthis
    real(RKIND), &
      dimension(-bufsize:,-bufsize:), intent(inout) :: rhs
    
    !*** local variables ***
    integer     :: i, j, k 
    integer     :: nbx, nby
    real(RKIND) :: threeOverh, hOverThree, twoh
    
    !*** initializations for the x direction ***
    nbx         = hsthis%nbx
    nby         = hsthis%nby
    threeOverh = 3._RKIND/hsthis%dx
    hOverThree = hsthis%dx/3._RKIND
    twoh       = 2._RKIND*hsthis%dx

    !-> Sweep down 
    do j = -2,nby+1
      hsthis%aux(-1,j) = rhs(-1,j)
      hsthis%aux(0,j)  = rhs(0,j) + hOverThree * hsthis%aux(-1,j)
      do i = 1,nbx-1
        hsthis%aux(i,j) = rhs(i,j) - hsthis%lmatx(i) * &
          hsthis%aux(i-1,j)
      end do
      hsthis%aux(nbx,j) = rhs(nbx,j)  &
        + threeOverh * hsthis%lmatx(nbx-1) * hsthis%aux(nbx-2,j) &
        - threeOverh * hsthis%lmatx(nbx) * hsthis%aux(nbx-1,j)
    enddo

    !-> Sweep up 
    do j = -2,nby+1 
      hsthis%bcoef(nbx,j)   = twoh*hsthis%aux(nbx,j)/hsthis%umatx(3)
      hsthis%bcoef(nbx-1,j) = (6._RKIND*hsthis%aux(nbx-1,j) - &
        hsthis%bcoef(nbx,j)) / hsthis%umatx(2)
      do i = nbx-2,1,-1
        hsthis%bcoef(i,j) = (6._RKIND*hsthis%aux(i,j) - &
          hsthis%bcoef(i+1,j)) * hsthis%lmatx(i+1)
      end do
      hsthis%bcoef(0,j)    = (6._RKIND*hsthis%aux(0,j) - &
        2._RKIND*hsthis%bcoef(1,j)) * hsthis%lmatx(1)
      hsthis%bcoef(-1,j)   = hsthis%bcoef(1,j)-twoh*hsthis%aux(-1,j)
    enddo
    do j=-2,nby+1
      hsthis%bcoef(-2,j)   = 6._RKIND*rhs(-2,j) - &
        4._RKIND*hsthis%bcoef(-1,j) - hsthis%bcoef(0,j)
      hsthis%bcoef(nbx+1,j) = 6._RKIND*rhs(nbx+1,j) - &
        4._RKIND*hsthis%bcoef(nbx,j) - hsthis%bcoef(nbx-1,j)
    enddo

    !*** initializations for the y direction ***
    threeOverh = 3._RKIND/hsthis%dy
    hOverThree = hsthis%dy/3._RKIND
    twoh       = 2._RKIND*hsthis%dy

    !-> Sweep down
    do i = -2,nbx+1
      hsthis%auy(i,-1) = hsthis%bcoef(i,-1)
      hsthis%auy(i,0)  = hsthis%bcoef(i,0) + &
        hOverThree * hsthis%auy(i,-1)
      do j = 1, nby-1
        hsthis%auy(i,j) = hsthis%bcoef(i,j) - &
          hsthis%lmaty(j) * hsthis%auy(i,j-1)
      enddo
      hsthis%auy(i,nby) = hsthis%bcoef(i,nby) & 
        + threeOverh * hsthis%lmaty(nby-1) * hsthis%auy(i,nby-2) & 
        - threeOverh * hsthis%lmaty(nby) * hsthis%auy(i,nby-1)
    enddo

    !-> Sweep up
    do i = -2,nbx+1
      hsthis%coef(i,nby)   = twoh*hsthis%auy(i,nby)/hsthis%umaty(3)
      hsthis%coef(i,nby-1) = (6._RKIND*hsthis%auy(i,nby-1) - &
        hsthis%coef(i,nby)) / hsthis%umaty(2)
      do j = nby-2,1,-1
        hsthis%coef(i,j) = (6._RKIND*hsthis%auy(i,j) - &
          hsthis%coef(i,j+1)) * hsthis%lmaty(j+1)
      enddo
      hsthis%coef(i,0)    = (6._RKIND*hsthis%auy(i,0) - &
        2._RKIND*hsthis%coef(i,1)) * hsthis%lmaty(1)
      hsthis%coef(i,-1)   = hsthis%coef(i,1)-twoh*hsthis%auy(i,-1)
    enddo
    do i=-2,nbx+1
      hsthis%coef(i,-2)   = 6._RKIND*hsthis%bcoef(i,-2) - &
        4._RKIND*hsthis%coef(i,-1) - hsthis%coef(i,0)
      hsthis%coef(i,nby+1) = 6._RKIND*hsthis%bcoef(i,nby+1) - &
        4._RKIND*hsthis%coef(i,nby) - hsthis%coef(i,nby-1)
    enddo
  end subroutine hermite


  !------------------------------------------------------ 
  ! Compute spline basis in x direction
  !------------------------------------------------------  
  subroutine splinex_basis(hsthis,xprev,xstar,xnext,sbase)
    type(splinehh)              , intent(in)    :: hsthis          
    real(RKIND)                 , intent(in)    :: xprev, xstar
    real(RKIND)                 , intent(in)    :: xnext
    real(RKIND), dimension(-1:2), intent(inout) :: sbase

    real(RKIND) :: coeff, d_prev, d_next

    coeff     = hsthis%invhx*hsthis%invhx*hsthis%invhx
    d_prev    = xstar-xprev
    d_next    = xnext-xstar

    sbase(2)  = coeff*d_prev*d_prev*d_prev*0.166666666666667_RKIND
    sbase(1)  = (1._RKIND + coeff * &
      (3._RKIND*d_prev * &
      (hsthis%dx*hsthis%dx+hsthis%dx*d_prev-d_prev*d_prev))) * &
      0.166666666666667_RKIND
    sbase(0)  = (1._RKIND + coeff * &
      (3._RKIND*d_next * &
      (hsthis%dx*hsthis%dx+hsthis%dx*d_next-d_next*d_next))) * &
      0.166666666666667_RKIND
    sbase(-1) = (coeff*d_next*d_next*d_next) * &
      0.166666666666667_RKIND  
  end subroutine splinex_basis


  !------------------------------------------------------ 
  ! Compute spline basis in y direction
  !------------------------------------------------------  
  subroutine spliney_basis(hsthis,xprev,xstar,xnext,sbase)
    type(splinehh)              , intent(in)    :: hsthis
    real(RKIND)                 , intent(in)    :: xprev, xstar
    real(RKIND)                 , intent(in)    :: xnext
    real(RKIND), dimension(-1:2), intent(inout) :: sbase

    real(RKIND) :: coeff, d_prev, d_next

    coeff     = hsthis%invhy*hsthis%invhy*hsthis%invhy
    d_prev    = xstar-xprev
    d_next    = xnext-xstar

    sbase(2)  = coeff*d_prev*d_prev*d_prev*0.166666666666667_RKIND
    sbase(1)  = (1._RKIND + coeff * &
      (3._RKIND*d_prev*(hsthis%dy*hsthis%dy + &
      hsthis%dy*d_prev-d_prev*d_prev))) * &
      0.166666666666667_RKIND
    sbase(0)  = (1._RKIND + coeff * &
      (3._RKIND*d_next*(hsthis%dy*hsthis%dy + &
      hsthis%dy*d_next-d_next*d_next))) * &
      0.166666666666667_RKIND
    sbase(-1) = (coeff*d_next*d_next*d_next) * &
      0.166666666666667_RKIND
  end subroutine spliney_basis
  

  !------------------------------------------------------ 
  ! 2D interpolation
  !------------------------------------------------------  
  subroutine interpol2d(hsthis,xstar1,xstar2,finterpol)
    type(splinehh), intent(in)    :: hsthis                
    real(RKIND)   , intent(in)    :: xstar1, xstar2
    real(RKIND)   , intent(inout) :: finterpol

    integer                      :: ipos1, ipos2
    integer                      :: i, j
    real(RKIND), dimension(-1:2) :: sbase1, sbase2

    !*** array position location in r ***
    ! -> if the particle is out of the domain 
    !    is put on the boundaries
    ipos1 = 1 + hsthis%invhx*(xstar1-hsthis%xgrid(1))
    ipos1 = max(min(ipos1,hsthis%nbx+1),1)

    !*** array position location in theta ***
    ! -> if the particle is out of the domain 
    !    is put on the boundaries
    ipos2 = 1 + hsthis%invhy*(xstar2-hsthis%ygrid(1))
    ipos2 = max(min(ipos2,hsthis%nby+1),1)

    !*** calculation of cubic spline basis ***
    call splinex_basis(hsthis,hsthis%xgrid(ipos1), &
      xstar1,hsthis%xgrid(ipos1+1),sbase1)
    call spliney_basis(hsthis,hsthis%ygrid(ipos2), &
      xstar2,hsthis%ygrid(ipos2+1),sbase2)

    !*** computation of f(x*,y*) ***
    finterpol = 0._RKIND
    ipos1     = ipos1-2
    ipos2     = ipos2-2
    finterpol = finterpol + &
      hsthis%coef(ipos1-1,ipos2-1)*sbase1(-1)*sbase2(-1) + &
      hsthis%coef(ipos1,ipos2-1)*sbase1(0)*sbase2(-1) + &
      hsthis%coef(ipos1+1,ipos2-1)*sbase1(1)*sbase2(-1) + &
      hsthis%coef(ipos1+2,ipos2-1)*sbase1(2)*sbase2(-1)
    finterpol = finterpol + &
      hsthis%coef(ipos1-1,ipos2)*sbase1(-1)*sbase2(0) + &
      hsthis%coef(ipos1,ipos2)*sbase1(0)*sbase2(0) + &
      hsthis%coef(ipos1+1,ipos2)*sbase1(1)*sbase2(0) + &
      hsthis%coef(ipos1+2,ipos2)*sbase1(2)*sbase2(0)
    finterpol = finterpol + &
      hsthis%coef(ipos1-1,ipos2+1)*sbase1(-1)*sbase2(1) + &
      hsthis%coef(ipos1,ipos2+1)*sbase1(0)*sbase2(1) + &
      hsthis%coef(ipos1+1,ipos2+1)*sbase1(1)*sbase2(1) + &
      hsthis%coef(ipos1+2,ipos2+1)*sbase1(2)*sbase2(1)
    finterpol = finterpol + &
      hsthis%coef(ipos1-1,ipos2+2)*sbase1(-1)*sbase2(2) + &
      hsthis%coef(ipos1,ipos2+2)*sbase1(0)*sbase2(2) + &
      hsthis%coef(ipos1+1,ipos2+2)*sbase1(1)*sbase2(2) + &
      hsthis%coef(ipos1+2,ipos2+2)*sbase1(2)*sbase2(2)
  end subroutine interpol2d
end module cubic_splines_module
