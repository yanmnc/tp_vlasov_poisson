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
! file : physics.f90
!
!  Used for all physic output computation
!------------------------------------------------

module physics_module
  use prec_const
  use globals, only : time_evol, &
    Phi_evol, E_evol, dens, velo, temp, flux, &
    f1_evol, nbion_evol, &
    Enkin_evol, Enpot_evol, &
    entropy_evol, L1_norm_evol, L2_norm_evol

  use geometry_class
  use fdistribu2d_class

  implicit none
  
  !******************************
  contains
  !******************************

  !---------------------------------------- 
  ! memory allocation for results saving
  !---------------------------------------- 
  subroutine results_allocate(geom,nb_save)
    type(geometry), intent(in) :: geom
    integer       , intent(in) :: nb_save
    
    allocate(time_evol(0:nb_save))
    allocate(Phi_evol(0:nb_save,0:geom%Nx))
    allocate(E_evol(0:nb_save,0:geom%Nx))
    allocate(dens(0:geom%Nx))
    allocate(velo(0:geom%Nx))
    allocate(temp(0:geom%Nx))
    allocate(flux(0:geom%Nx))
    allocate(f1_evol(0:nb_save,0:geom%Nx,0:geom%Nv))
    allocate(Enkin_evol(0:nb_save))
    allocate(Enpot_evol(0:nb_save))
    allocate(nbion_evol(0:nb_save))
    allocate(entropy_evol(0:nb_save))
    allocate(L1_norm_evol(0:nb_save))
    allocate(L2_norm_evol(0:nb_save))
    call results_clearing
  end subroutine results_allocate
  

  !-------------------------------------------- 
  ! Array clearing, initialisation to 0
  !--------------------------------------------   
  subroutine results_clearing()    
    time_evol    = 0._RKIND
    Phi_evol     = 0._RKIND
    E_evol       = 0._RKIND
    dens         = 0._RKIND
    velo         = 0._RKIND
    temp         = 0._RKIND
    flux         = 0._RKIND
    f1_evol      = 0._RKIND
    Enkin_evol   = 0._RKIND
    Enpot_evol   = 0._RKIND    
    nbion_evol   = 0._RKIND
    entropy_evol = 0._RKIND
    L1_norm_evol = 0._RKIND
    L2_norm_evol = 0._RKIND
  end subroutine results_clearing
  
  
  !---------------------------------------------------- 
  ! deallocation of the memory used for results saving
  !----------------------------------------------------
  subroutine results_deallocate() 
    deallocate(time_evol)  
    deallocate(Phi_evol)
    deallocate(E_evol)
    deallocate(dens)
    deallocate(velo)
    deallocate(temp)
    deallocate(flux)
    deallocate(f1_evol)
    deallocate(Enkin_evol) 
    deallocate(Enpot_evol)
    deallocate(nbion_evol)
    deallocate(entropy_evol)
    deallocate(L1_norm_evol)
    deallocate(L2_norm_evol)
  end subroutine results_deallocate


  !-----------------------------------------------------
  ! Used to write info on the screen at each iteration
  !-----------------------------------------------------
  subroutine write_info(iter_glob_num,time_glob_num,dt_num)
    use globals, only  : nbiter
    integer    , intent(in) :: iter_glob_num
    real(RKIND), intent(in) :: time_glob_num
    real(RKIND), intent(in) :: dt_num
    
    character(len=300) :: string
    character(6)       :: iterglobnum_char
    character(17)      :: dt_char     
    character(6)       :: nbiter_char
    character(17)      :: timeglobnum_char
    character(17)      :: Phiiter_char    
    
    write(nbiter_char,'(I6)')         nbiter
    write(iterglobnum_char,'(I6)')    iter_glob_num
    write(timeglobnum_char,'(E17.8)') time_glob_num 
    write(dt_char,'(E17.8)')          dt_num
    
    write(string,*) '*** ITERATION ', &
      '('//iterglobnum_char// ' /' &
      //nbiter_char//') : time = ' //timeglobnum_char// &
      ' : dt = ' //dt_char
    write(*,*) string
  end subroutine write_info
     

  !--------------------------------------------------------
  ! Computes the 4 first fluid moments:
  !  . dens = \int f(x,v) dv
  !  . velo = \int v f(x,v) dv / dens
  !  . temp = \int (v-velo)**2 f(x,v) dv / dens
  !  . flux = \int (v-velo)**3 f(x,v) dv
  !--------------------------------------------------------
  subroutine compute_fluid_moments(geom,f,idiag_num)
    use interpolation_module, only : compute_intdv
    type(geometry)   , intent(in) :: geom
    type(fdistribu2d), intent(in) :: f
    integer          , intent(in) :: idiag_num
    
    real(RKIND), dimension(0:geom%Nx,0:geom%Nv) :: integrand
    real(RKIND), dimension(0:geom%Nx)           :: momentum, pressure
    integer                                     :: ix, iv

    !*** compute dens = \int f(x,v) dv ****
    call compute_intdv(geom,f%values(0:geom%Nx,0:geom%Nv),dens)
    !*** compute velo = \int v f(x,v) dv / dens ****
    do iv=0,geom%Nv
      do ix=0,geom%Nx
        integrand(ix,iv) = f%values(ix,iv)*geom%vgrid(iv)
      end do
    end do
    call compute_intdv(geom,integrand(0:geom%Nx,0:geom%Nv),momentum)
    velo = momentum / dens
    !*** compute temp = \int (v-velo)**2 f(x,v) dv / dens ****
    do iv=0,geom%Nv
      do ix=0,geom%Nx
        integrand(ix,iv) = f%values(ix,iv)* &
          (geom%vgrid(iv)-velo(ix))*(geom%vgrid(iv)-velo(ix))
      end do
    end do
    call compute_intdv(geom,integrand(0:geom%Nx,0:geom%Nv),pressure)
    temp = pressure / dens
    !*** compute flux = \int (v-velo)**3 f(x,v) dv ****
    do iv=0,geom%Nv
      do ix=0,geom%Nx
        integrand(ix,iv) = f%values(ix,iv)* &
          (geom%vgrid(iv)-velo(ix))*(geom%vgrid(iv)-velo(ix))*&
          (geom%vgrid(iv)-velo(ix))
      end do
    end do
    call compute_intdv(geom,integrand(0:geom%Nx,0:geom%Nv),flux)
  end subroutine compute_fluid_moments
   

  !--------------------------------------------------------
  ! Computes the total number of ions
  !  . nbions = \int(f(x,v) dx dv)
  !           = \int(ni(x) dx
  !--------------------------------------------------------
  subroutine compute_nbions(geom,f,idiag_num)
    use interpolation_module, only : compute_intdxdv
    type(geometry)   , intent(in) :: geom
    type(fdistribu2d), intent(in) :: f
    integer          , intent(in) :: idiag_num
    
    real(RKIND) :: nbions

    !*** compute nb_ions = \int(f(x,v) dx dv ****
    call compute_intdxdv(geom, &
      f%values(0:geom%Nx,0:geom%Nv),nbions)
    nbion_evol(idiag_num) = nbions
  end subroutine compute_nbions


  !--------------------------------------------------
  ! Computes the entropy
  !  H(t) = -sum_i fi(t)ln(fi(t))
  !--------------------------------------------------
  subroutine compute_entropy(geom,f,entropy)
    type(geometry)   , intent(in)  :: geom
    type(fdistribu2d), intent(in)  :: f
    real(RKIND)      , intent(out) :: entropy
    
    real(RKIND) :: f_tmp, entrop_tmp
    integer     :: ix, iv

    entropy = 0._RKIND
    do iv = 0,geom%Nv
      do ix = 0,geom%Nx
        f_tmp      = abs(f%values(ix,iv))
        entrop_tmp = -1._RKIND*f_tmp*log(f_tmp)
        entropy    = entropy + entrop_tmp * &
          geom%coeff_intdx(ix)*geom%coeff_intdv(iv)
      end do
    end do
  end subroutine compute_entropy


  !--------------------------------------------------
  ! Computes the Lp-norm
  !  Lp(t) = sum_i abs(fi(t))^p
  !--------------------------------------------------
  subroutine compute_Lpnorm(geom,f,p,Lp_norm)
    type(geometry)   , intent(in)  :: geom
    type(fdistribu2d), intent(in)  :: f
    real(RKIND)      , intent(in)  :: p
    real(RKIND)      , intent(out) :: Lp_norm

    real(RKIND) :: f_tmp, Lpnorm_tmp
    integer     :: ix, iv

    Lp_norm = 0._RKIND
    do iv = 0,geom%Nv
      do ix = 0,geom%Nx
        f_tmp      = abs(f%values(ix,iv))
        Lpnorm_tmp = f_tmp**p
        Lp_norm    = Lp_norm + Lpnorm_tmp * &
          geom%coeff_intdx(ix)*geom%coeff_intdv(iv)
      end do
    end do
  end subroutine compute_Lpnorm


  !---------------------------------------------- 
  ! Call the computing of all physic values
  !---------------------------------------------- 
  subroutine compute_physics_value(geom,f,Phi,E,idiag_num)
    use energy_module
    type(geometry)            , intent(in) :: geom
    type(fdistribu2d)         , intent(in) :: f
    real(RKIND), dimension(0:), intent(in) :: Phi
    real(RKIND), dimension(0:), intent(in) :: E
    integer                   , intent(in) :: idiag_num

    call compute_fluid_moments(geom,f,idiag_num)
    call compute_nbions(geom,f,idiag_num)
    call compute_kinetic_energy(geom,f,idiag_num)
    call compute_potential_energy(geom,f,Phi,idiag_num)
    call compute_entropy(geom,f,entropy_evol(idiag_num))
    call compute_Lpnorm(geom,f,1._RKIND,L1_norm_evol(idiag_num))
    call compute_Lpnorm(geom,f,2._RKIND,L2_norm_evol(idiag_num))
  end subroutine compute_physics_value


  !---------------------------------------------- 
  ! array saving for all the diagnostics in time 
  !  (for leap-frog)
  !---------------------------------------------- 
  subroutine diagnostics(iter_time,idiag_num,geom,f1,Phi,E)
    use globals, only : f2D_saving
    use ascii_saving_module, only : ascii_Phi1D_write, &
      ascii_fluid_moments_write, ascii_f1D_write, ascii_f2D_write
    real(RKIND)               , intent(in)    :: iter_time
    integer                   , intent(in)    :: idiag_num
    type(geometry)            , intent(in)    :: geom
    type(fdistribu2d)         , intent(inout) :: f1
    real(RKIND), dimension(0:), intent(in)    :: Phi
    real(RKIND), dimension(0:), intent(in)    :: E

    !*** Computation of all physic outputs ***
    call compute_physics_value(geom,f1,Phi,E,idiag_num)

    !*** writting of Phi(0:Nx) ***
    call ascii_Phi1D_write(idiag_num,geom,Phi)

    !*** writting of fluid moments dens(0:Nx), velo(0:Nx), temp(0:Nx) ***
    call ascii_fluid_moments_write(idiag_num,geom,dens,velo,temp,flux)

    !*** Distribution function saving ***
    call ascii_f1D_write(idiag_num,geom,f1%values)
    if (f2D_saving) & 
      call ascii_f2D_write(idiag_num,geom,f1%values)

    !*** diagnostics on the electric potential ***
    time_evol(idiag_num)   = iter_time
    Phi_evol(idiag_num,:)  = Phi(:)
    E_evol(idiag_num,:)    = E(:)
    f1_evol(idiag_num,:,:) = f1%values(:,:)
    write(*,*) 'ni = ', nbion_evol(idiag_num)
  end subroutine diagnostics
end module physics_module
