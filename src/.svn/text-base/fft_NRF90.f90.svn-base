!-------------------------------------------------------
! file : fft_NRF90.f90 
! date : 30/04/2001
!  - fast fourier transform
! ( Numerical Recipes Version for fortran 90 )
!
! ATTENTION : to have a correct Fourier Transform, you 
!  must not take into account the last periodic point
! ex : if the function f(1:nx) is periodic, you have to
!       use the (nx-1) first points, and (nx-1) must
!       be a power of 2
!-------------------------------------------------------
module fftNRF90_module
  implicit none
  private
  public :: fft, fft_inv

  integer, parameter :: npar_arth  = 16
  integer, parameter :: npar2_arth = 8

  interface assert
     module procedure assert1
  end interface

  interface assert_eq
     module procedure assert_eq2
  end interface

  interface arth
     module procedure arth_d, arth_i
  end interface

  interface swap
     module procedure swap_zv
  end interface

  interface fourrow
     module procedure fourrow_dp
  end interface

  interface four1
     module procedure four1_dp
  end interface

  interface realft
     module procedure realft_dp
  end interface

  interface fft
     module procedure fft_r, fft_2D_r, fft_c, fft2_c, fft2_3D_c
  end interface

  interface fft_inv
     module procedure fftinv_r, fftinv_2D_r, fftinv_c, &
       fftinv2_c, fftinv2_3D_c
  end interface


  !******************************
  contains
  !******************************

  !------------------------------------------------------
  ! numerical recipes : functions extracted of nrutil.f90
  !------------------------------------------------------
  subroutine assert1(n1,string)
    character(len=*), intent(in) :: string
    logical         , intent(in) :: n1

    if (.not. n1) then
      write (*,*) 'nrerror: an assertion failed with this tag:', &
        string
      stop 'program terminated by assert1'
    end if
  end subroutine assert1

  !------------------------------------------------------
  function assert_eq2(n1,n2,string)
    character(len=*), intent(in) :: string
    integer         , intent(in) :: n1,n2

    integer :: assert_eq2

    if (n1 == n2) then
      assert_eq2 = n1
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag:', &
        string
      stop 'program terminated by assert_eq2'
    end if
  end function assert_eq2

  !------------------------------------------------------
  function arth_d(first,increment,n)
    use prec_const
    real(RKIND), intent(in) :: first, increment
    integer    , intent(in) :: n

    real(RKIND), dimension(n) :: arth_d
    integer                   :: k, k2
    real(RKIND)               :: temp

    if (n > 0) arth_d(1) = first
    if (n <= npar_arth) then
      do k = 2,n
        arth_d(k) = arth_d(k-1)+increment
      end do
    else
      do k = 2,npar2_arth
        arth_d(k) = arth_d(k-1)+increment
      end do
      temp = increment*npar2_arth
      k    = npar2_arth
      do
        if (k >= n) exit
        k2                    = k+k
        arth_d(k+1:min(k2,n)) = temp+arth_d(1:min(k,n-k))
        temp                  = temp+temp
        k                     = k2
      end do
    end if
  end function arth_d

  !------------------------------------------------------
  function arth_i(first,increment,n)
    use prec_const
    integer, intent(in) :: first, increment, n

    integer, dimension(n) :: arth_i
    integer               :: k, k2, temp

    if (n > 0) arth_i(1) = first
    if (n <= npar_arth) then
      do k = 2,n
        arth_i(k) = arth_i(k-1)+increment
      end do
    else
      do k = 2,npar2_arth
        arth_i(k) = arth_i(k-1)+increment
      end do
      temp = increment*npar2_arth
      k    = npar2_arth
      do
        if (k >= n) exit
        k2                    = k+k
        arth_i(k+1:min(k2,n)) = temp+arth_i(1:min(k,n-k))
        temp                  = temp+temp
        k                     = k2
      end do
    end if
  end function arth_i

  !------------------------------------------------------
  function zroots_unity(n,nn)
    use prec_const
    integer, intent(in) :: n, nn

    complex(CKIND), dimension(nn) :: zroots_unity
    integer                       :: k
    real(RKIND)                   :: theta

    zroots_unity(1) = 1.0_RKIND
    theta           = TWOPI/n
    k               = 1
    do
      if (k >= nn) exit
      zroots_unity(k+1)             = cmplx(cos(k*theta),sin(k*theta),CKIND)
      zroots_unity(k+2:min(2*k,nn)) = zroots_unity(k+1)*&
        zroots_unity(2:min(k,nn-k))
      k                             = 2*k
    end do
  end function zroots_unity

  !------------------------------------------------------
  subroutine swap_zv(a,b)
    use prec_const
    complex(CKIND), dimension(:), intent(inout) :: a,b

    complex(CKIND), dimension(size(a)) :: dum

    dum = a
    a   = b
    b   = dum
  end subroutine swap_zv


  !------------------------------------------------------
  ! numerical recipes : functions extracted of fourrow.f90
  !------------------------------------------------------
  subroutine fourrow_dp(data,isign)
    use prec_const
    complex(CKIND), dimension(:,:), intent(inout) :: data
    integer                       , intent(in)    :: isign

    integer        :: n, i, istep, j, m, mmax, n2
    real(RKIND)    :: theta
    complex(CKIND) :: w, wp
    complex(CKIND) :: ws
    complex(CKIND), dimension(size(data,1)) :: temp

    n = size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
    n2 = n/2
    j  = n2
    do i = 1,n-2
      if (j > i) call swap(data(:,j+1),data(:,i+1))
      m = n2
      do
        if (m < 2 .or. j < m) exit
        j = j-m
        m = m/2
      end do
      j = j+m
    end do
    mmax = 1
    do
      if (n <= mmax) exit
      istep = 2*mmax
      theta = PI/(isign*mmax)
      wp    = cmplx(-2.0_RKIND*sin(0.5_RKIND*theta)**2,sin(theta),kind=CKIND)
      w     = cmplx(1.0_RKIND,0.0_RKIND,kind=CKIND)
      do m = 1,mmax
        ws = w
        do i = m,n,istep
          j         = i+mmax
          temp      = ws*data(:,j)
          data(:,j) = data(:,i)-temp
          data(:,i) = data(:,i)+temp
        end do
        w = w*wp+w
      end do
      mmax = istep
    end do
  end subroutine fourrow_dp


  !------------------------------------------------------
  ! numerical recipes : functions extracted of four1.f90
  !------------------------------------------------------
  subroutine four1_dp(data,isign)
    use prec_const
    complex(CKIND), dimension(:), intent(inout) :: data
    integer                     , intent(in)    :: isign

    complex(CKIND), dimension(:,:), allocatable :: dat, temp
    complex(CKIND), dimension(:)  , allocatable :: w, wp
    real(RKIND)   , dimension(:)  , allocatable :: theta
    integer                                     :: n, m1, m2, j

    n = size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_dp')
    m1 = 2**ceiling(0.5_RKIND*log(real(n,RKIND))/0.693147_RKIND)
    m2 = n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat = reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta = arth(0,isign,m1)*TWOPI/n
    wp    = cmplx(-2.0_RKIND*sin(0.5_RKIND*theta)**2,sin(theta),kind=CKIND)
    w     = cmplx(1.0_RKIND,0.0_RKIND,kind=CKIND)
    do j = 2,m2
      w        = w*wp+w
      dat(:,j) = dat(:,j)*w
    end do
    temp = transpose(dat)
    call fourrow(temp,isign)
    data = reshape(temp,shape(data))
    deallocate(dat,w,wp,theta,temp)
  end subroutine four1_dp


  !------------------------------------------------------
  ! numerical recipes : functions extracted of realft.f90
  !------------------------------------------------------
  subroutine realft_dp(data,isign,zdata)
    use prec_const
    real(RKIND)   , dimension(:), intent(inout)    :: data
    integer                     , intent(in)       :: isign
    complex(CKIND), dimension(:), optional, target :: zdata

    integer        :: n, ndum, nh, nq
    complex(CKIND) :: z
    real(RKIND)    :: c1=0.5_RKIND,c2
    complex(CKIND), dimension(size(data)/4)   :: w
    complex(CKIND), dimension(size(data)/4-1) :: h1, h2
    complex(CKIND), dimension(:), pointer     :: cdata

    n = size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_dp')
    nh = n/2
    nq = n/4
    if (present(zdata)) then
      ndum  = assert_eq(n/2,size(zdata),'realft_dp')
      cdata => zdata
      if (isign == 1) cdata = cmplx(data(1:n-1:2),data(2:n:2),kind=CKIND)
    else
      allocate(cdata(n/2))
      cdata = cmplx(data(1:n-1:2),data(2:n:2),kind=CKIND)
    end if
    if (isign == 1) then
      c2 = -0.5_RKIND
      call four1(cdata,+1)
    else
      c2 = 0.5_RKIND
    end if
    w                 = zroots_unity(sign(n,isign),n/4)
    w                 = cmplx(-aimag(w),real(w),kind=CKIND)
    h1                = c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
    h2                = c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
    cdata(2:nq)       = h1+w(2:nq)*h2
    cdata(nh:nq+2:-1) = conjg(h1-w(2:nq)*h2)
    z                 = cdata(1)
    if (isign == 1) then
      cdata(1) = cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=CKIND)
    else
      cdata(1) = cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=CKIND)
      call four1(cdata,-1)
    end if
    if (present(zdata)) then
      if (isign /= 1) then
        data(1:n-1:2) = real(cdata)
        data(2:n:2)   = aimag(cdata)
      end if
    else
      data(1:n-1:2) = real(cdata)
      data(2:n:2)   = aimag(cdata)
      deallocate(cdata)
    end if
  end subroutine realft_dp


  !****************************************************
  !  FFT in ONE dimension
  !****************************************************

  !----------------------------------------------------
  ! FFT 1D computation (with real function)
  !  input/output array = 1D array
  !
  ! when 'npoints' is specified, the output data are
  !  stored in an array of size (npoints+2) to keep a 
  !  symetry in the storage of the real parts and the 
  !  imaginary parts of the (npoints/2) points
  !----------------------------------------------------  
  subroutine fft_r(data,npoints)
    use prec_const
    real(RKIND), dimension(:), intent(inout)          :: data
    integer                  , intent(in)  , optional :: npoints

    integer :: Ndata

    if (present(npoints)) then
      Ndata = size(data)
      if (Ndata.lt.(npoints+2)) then
        write(*,*) Ndata,' the size of the array must be superior to', &
          npoints+2
        stop 'program terminated by fft_r'
      end if
      call realft(data(1:npoints),1) 
      data(npoints+1) = data(2)
      data(2)         = 0._RKIND
      data(npoints+2) = 0._RKIND
    else
      call realft(data,1)
    end if
  end subroutine fft_r


  !----------------------------------------------------
  ! FFT 1D inverse computation (with real function)
  !  input/output array = 1D array
  !
  ! when 'npoints' is specified, the output data are
  !  stored in an array of size (npoints+2) to keep a 
  !  symetry in the storage of the real parts and the 
  !  imaginary parts of the (npoints/2) points
  !----------------------------------------------------  
  subroutine fftinv_r(data,npoints)
    use prec_const
    real(RKIND), dimension(:), intent(inout)        :: data
    integer                  , intent(in), optional :: npoints
    integer     :: Ndata

    Ndata = size(data)
    if (present(npoints)) then
      if (Ndata.lt.(npoints+2)) then
        write(*,*) Ndata,' the size of the array must be superior to', &
          npoints+2
        stop 'program terminated by fftinv_r'
      end if
      data(2)         = data(npoints+1)
      data(npoints+1) = 0._RKIND
      data(npoints+2) = 0._RKIND
      call realft(data(1:npoints),-1)
      data = 2._RKIND*data/npoints
    else
      call realft(data,-1)
      data = 2._RKIND*data/Ndata
    end if
  end subroutine fftinv_r


  !----------------------------------------------------
  ! FFT computation (with real function)
  !  input/output array = 1D array
  ! 
  ! the FFT is done on the first index
  !
  ! when 'npoints' is specified, the output data are
  !  stored in an array of size (npoints+2) to keep a 
  !  symetry in the storage of the real parts and the 
  !  imaginary parts of the (npoints/2) points
  !----------------------------------------------------  
  subroutine fft_2D_r(data,npoints)
    use prec_const
    real(RKIND), dimension(:,:), intent(inout)           :: data
    integer                    , intent(in)   , optional :: npoints
       
    integer :: j, Ndata1, Ndata2

    Ndata1 = size(data,1)
    Ndata2 = size(data,2)
    if (present(npoints)) then
      if (Ndata1.lt.(npoints+2)) then
        write(*,*) Ndata1,' the size of the array must be superior to', &
          npoints+2
        stop 'program terminated by fft_2D_r'
      end if
      do j = 1,Ndata2
        call realft(data(1:npoints,j),1)
      end do
      data(npoints+1,:) = data(2,:)
      data(2,:)         = 0._RKIND
      data(npoints+2,:) = 0._RKIND
    else
      do j = 1,Ndata2
        call realft(data(:,j),1)
      end do
    end if
  end subroutine fft_2D_r


  !----------------------------------------------------
  ! FFT inverse computation (with real function)
  !  input/output array = 1D array
  ! 
  ! the FFT is done on the first index
  !
  ! when 'npoints' is specified, the output data are
  !  stored in an array of size (npoints+2) to keep a 
  !  symetry in the storage of the real parts and the 
  !  imaginary parts of the (npoints/2) points
  !----------------------------------------------------  
  subroutine fftinv_2D_r(data,npoints)
    use prec_const
    real(RKIND), dimension(:,:), intent(inout)           :: data
    integer                    , intent(in)   , optional :: npoints

    integer :: j, Ndata1, Ndata2

    Ndata1 = size(data,1)
    Ndata2 = size(data,2)
    if (present(npoints)) then
      if (Ndata1.lt.(npoints+2)) then
        write(*,*) Ndata1,' the size of the array must be superior to', &
          npoints+2
        stop 'program terminated by fftinv_2D_r'
      end if
      data(2,:)         = data(npoints+1,:)
      data(npoints+1,:) = 0._RKIND
      data(npoints+2,:) = 0._RKIND
      do j = 1,Ndata2
        call realft(data(1:npoints,j),-1)
      end do
      data = 2._RKIND*data/npoints
    else
      do j = 1,Ndata2
        call realft(data(:,j),-1)
      end do
      data = 2._RKIND*data/Ndata1      
    end if
  end subroutine fftinv_2D_r


  !----------------------------------------------------
  ! FFT 1D computation (with complex function)
  !----------------------------------------------------  
  subroutine fft_c(data)
    use prec_const
    complex(CKIND), dimension(:), intent(inout) :: data
    
    integer     :: i, Ndata
    real(RKIND) :: rem

    Ndata = size(data)
    rem   = modulo(log(real(Ndata)),real(log(2._RKIND)))
    if (rem.ne.0._RKIND) then
      write(*,*) 'Warning in fft_c : ',Ndata, ' is not a power of 2'
      stop
    end if
    
    call four1(data,1)
  end subroutine fft_c


  !----------------------------------------------------
  ! FFT 1D inverse computation (with complex function)
  !----------------------------------------------------  
  subroutine fftinv_c(data)
    use prec_const
    complex(CKIND), dimension(:), intent(inout) :: data

    integer     :: i, Ndata
    real(RKIND) :: rem

    Ndata = size(data)
    rem   = modulo(log(real(Ndata)),real(log(2._RKIND)))
    if (rem.ne.0._RKIND) then
      write(*,*) 'Warning in fftinv_c : ',Ndata, ' is not a power of 2'
      stop
    end if

    call four1(data,-1)
    data = data/Ndata
  end subroutine fftinv_c


  !****************************************************
  !  FFT in TWO dimensions
  !****************************************************

  !--------------------------------------------------------
  ! Numerical Recipes : functions extracted of fourcol.f90
  !--------------------------------------------------------
  subroutine fourcol(data,isign)
    use prec_const
    implicit none
    complex(CKIND), dimension(:,:), intent(inout) :: data
    integer                       , intent(in)    :: isign

    integer        :: n, i, istep, j, m, mmax, n2
    real(RKIND)    :: theta
    complex(RKIND) :: w, wp
    complex(CKIND) :: ws
    complex(CKIND), dimension(size(data,2)) :: temp

    n = size(data,1)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourcol')
    n2 = n/2
    j  = n2
    do i = 1,n-2
      if (j > i) call swap(data(j+1,:),data(i+1,:))
      m = n2
      do
        if (m < 2 .or. j < m) exit
        j = j-m
        m = m/2
      end do
      j = j+m
    end do
    mmax = 1
    do
      if (n <= mmax) exit
      istep = 2*mmax
      theta = PI/(isign*mmax)
      wp    = cmplx(-2.0_RKIND*sin(0.5_RKIND*theta)**2,sin(theta),kind=RKIND)
      w     = cmplx(1.0_RKIND,0.0_RKIND,kind=RKIND)
      do m = 1,mmax
        ws = w
        do i = m,n,istep
          j         = i+mmax
          temp      = ws*data(j,:)
          data(j,:) = data(i,:)-temp
          data(i,:) = data(i,:)+temp
        end do
        w = w*wp+w
      end do
      mmax = istep
    end do
  end subroutine fourcol


  !-----------------------------------------------------------------------
  ! FFT 2D computation (with complex function)
  !  input/output array = 2D array
  !-----------------------------------------------------------------------
  subroutine fft2_c(f,g)
    use prec_const
    complex(CKIND), dimension (:,:), intent(inout) :: f
    complex(CKIND), dimension (:,:), intent(out)   :: g

    !***  local variables and arrays
    integer :: n1, n2, n1_f, n2_f, n1_g, n2_g
    integer :: isign

    !*** check the array sizes ***
    n1_f = size(f,1); n2_f = size(f,2)
    n1_g = size(g,1); n2_g = size(g,2)
    if ((n1_f.ne.n2_g).or.(n2_f.ne.n1_g)) then
      write(*,*) 'n1_f = ',n1_f,' n2_f = ',n2_f, &
        ' n1_g = ',n1_g,' n2_g = ',n2_g
      stop' problem with array size in fft2_c'
    end if
    n1  = n1_f
    n2  = n2_f

    !*** FFT in 2 dimensions ***
    isign = 1
    call fourcol(f(1:n1,1:n2),isign)
    g = transpose(f)
    call fourcol(g(1:n2,1:n1),isign)
  end subroutine fft2_c


  !-----------------------------------------------------------------------
  ! FFT 2D inverse computation (with complex function)
  !  input/output array = 2D array
  !-----------------------------------------------------------------------
  subroutine fftinv2_c(g,f)
    use prec_const
    complex(CKIND), dimension (:,:), intent(inout) :: g
    complex(CKIND), dimension (:,:), intent(out)   :: f

    !***  local variables and arrays
    integer :: n1, n2, n1_f, n2_f, n1_g, n2_g
    integer :: isign

    !*** check the array sizes ***
    n1_f = size(f,1); n2_f = size(f,2)
    n1_g = size(g,1); n2_g = size(g,2)
    if ((n1_f.ne.n2_g).or.(n2_f.ne.n1_g)) then
      write(*,*) 'n1_f = ',n1_f,' n2_f = ',n2_f, &
        ' n1_g = ',n1_g,' n2_g = ',n2_g
      stop' problem with array size in fftinv2_c'
    end if
    n1 = n1_f
    n2 = n2_f

    !*** FFT inv in 2 dimensions ***
    isign = -1
    call fourcol(g(1:n2,1:n1),isign)
    f = transpose(g)
    call fourcol(f(1:n1,1:n2),isign)
    f = f/(n1*n2)
  end subroutine fftinv2_c


  !-------------------------------------------------------------
  ! FFT 2D  computation (with complex function)
  !  input/output array = 3D array
  ! 
  ! the FFT 2D is done on the first and third indexes
  !
  !-------------------------------------------------------------
  subroutine fft2_3D_c(f,g)
    use prec_const
    complex(CKIND), dimension (:,:,:), intent(inout) :: f
    complex(CKIND), dimension (:,:,:), intent(out)   :: g

    !***  local variables and arrays
    integer :: n1, n2, n3, n1_f, n2_f, n3_f, n1_g, n2_g, n3_g
    integer :: isign, i, k

    !*** check the array sizes ***
    n1_f = size(f,1); n2_f = size(f,2) ; n3_f = size(f,3)
    n1_g = size(g,1); n2_g = size(g,2) ; n3_g = size(g,3)
    if ((n1_f.ne.n3_g).or.(n2_f.ne.n2_g).or.(n3_f.ne.n1_g)) then
      write(*,*) 'n1_f = ',n1_f,' n2_f = ',n2_f,' n3_f = ',n3_f, &
        ' n1_g = ',n1_g,' n2_g = ',n2_g,' n3_g = ',n3_g
      stop' problem with array size in fft2_3D_c'
    end if
    n1 = n1_f
    n2 = n2_f
    n3 = n3_f

    !*** FFT in 2 dimensions ***
    isign = 1
    do k = 1,n3
      call fourcol(f(1:n1,1:n2,k),isign)
    end do
    call transp(f,g)
    do i = 1,n1
      call fourcol(g(1:n3,1:n2,i),isign)
    end do
  end subroutine fft2_3D_c


  !-------------------------------------------------------------
  ! FFT 2D inverse computation (with complex function)
  !  input/output array = 3D array
  ! 
  ! the FFT 2D is done on the first and third indexes
  !
  !-------------------------------------------------------------
  subroutine fftinv2_3D_c(g,f)
    use prec_const
    complex(CKIND), dimension (:,:,:), intent(inout) :: g
    complex(CKIND), dimension (:,:,:), intent(out)   :: f

    !***  local variables and arrays
    integer :: n1, n2, n3, n1_f, n2_f, n3_f, n1_g, n2_g, n3_g
    integer :: isign, i, k

    !*** check the array sizes ***
    n1_f = size(f,1); n2_f = size(f,2) ; n3_f = size(f,3)
    n1_g = size(g,1); n2_g = size(g,2) ; n3_g = size(g,3)
    if ((n1_f.ne.n3_g).or.(n2_f.ne.n2_g).or.(n3_f.ne.n1_g)) then
      write(*,*) 'n1_f = ',n1_f,' n2_f = ',n2_f,' n3_f = ',n3_f, &
        ' n1_g = ',n1_g,' n2_g = ',n2_g,' n3_g = ',n3_g
      stop' problem with array size in fft2_3D_c'
    end if
    n1 = n1_f
    n2 = n2_f
    n3 = n3_f

    !*** FFT in 2 dimensions ***
    isign = -1
    do i = 1,n1
      call fourcol(g(1:n3,1:n2,i),isign)
    end do
    call transp(g,f)
    do k = 1,n3
      call fourcol(f(1:n1,1:n2,k),isign)
    end do
    f = f/(n1*n3)
  end subroutine fftinv2_3D_c


  !---------------------------------------------------------------
  ! transpose the first and the third indexes of a 3D array, i.e :
  !  - input  : f(1:n1,1:n2,1:n3)
  !  - output : g(1:n3,1:n2,1:n1) = f(1:n1,1:n2,1:n3)
  !---------------------------------------------------------------
  subroutine transp(f,g)
    use prec_const
    complex(CKIND), dimension (:,:,:), intent(in) :: f
    complex(CKIND), dimension (:,:,:), intent(out):: g
    
    integer :: j, n1, n2, n3
    
    n1 = size(f,1) ; n2 = size(f,2) ; n3 = size(f,3) ;  
    do j = 1,n2
      g(1:n3,j,1:n1) = transpose(f(1:n1,j,1:n3))
    end do
  end subroutine transp
end module fftNRF90_module



