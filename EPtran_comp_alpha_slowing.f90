!---------------------------------------------------------
! EPtran_comp_alpha_slowing.f90
!
! PURPOSE:
!  compute alpha slowing down profiles
!
!---------------------------------------------------------

subroutine EPtran_comp_alpha_slowing

  use EPtran_use_parameter
  !CONTROL                 !NBI_flag
  !INPUTS	                 !nr,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
  !                        !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
  !                        !Zeff,M_DT,E_alpha,Fpeak_T,Fpeak_n,Fpeak_ei
  !OUTPUTS                 !beta_ped_percent,T_ped,beta_N_glob,arho,
  !                        !rho_hat(i),q_rho(i),Rmaj_rho(i),rmin_rho(i),
  !                        !kappa_rho(i),delta_rho(i),T_i_rho(i),T_e_rho(i),
  !                        !n_i_rho(i),n_e_rho(i),n_alpha_rho(i),T_alpha_equiv_rho(i),
  !                        !E_c_hat_rho(i)
  !OTHER                   !pi,tau_ee(i),n_alpha_ave_rho(i),S0_rho(i),tau_s_rho(i)

  !--------------------------------------
  implicit none
  integer :: i

  integer :: ii

  integer :: i_en
  integer :: n_en
  integer :: n_en_max

  real :: Numer
  real :: Denom
  real :: G_D
  real :: En
  real :: EnoTe
  real :: Convec_Factor

  real :: Z1
  
  real :: a
  !real :: tau_ee
  real :: ln_lambda

  real :: delta
  real :: rho_t
  real :: rho_p

  real, dimension(nr) :: I2
  real, dimension(nr) :: I4

  real, dimension(nr) :: rho_hat_p
  !
  character(len=80) :: comment
  !--------------------------------------


  ! rho_hat grid
  ! rho_hat(1) = 0.
  ! do i = 1,nr
  !   rho_hat(i) = real(i-1)/real(nr-1)
  !   E_alpha_grid(i) = E_alpha*1.E3 - 0.01*i
  ! enddo

  if(NBI_flag) then
    !NBI
    m_alpha = 2.0
    z_alpha = 1.0
    Z1=(5./3.) &
     *(M_DT/4.)/(M_DT/2.5)
    !Z1=1. !the correct setting, should use it after testing
  else
    !alpha particles
    m_alpha = 4.0
    z_alpha = 2.0
    Z1 = 5./3.
  endif
  
! alpha injection energy normed cross over energy E_c_hat = E_c/E_aplha 
  E_c_hat_rho(:) = 0.
  do i = 1,nr
    E_c_hat_rho(i) = (T_e_rho(i)/E_alpha)*1.E-3*(m_alpha*1836.)**(1./3.)*(3.*sqrt(pi)*Z1/4.)**(2./3.)

    a=sqrt(E_c_hat_rho(i))

    I2(i) = 1./3.*log((1.+a**3)/a**3)

    I4(i) = 1./2.-a**2*(1./6.*log((1.-a+a**2)/(1.+a)**2)+1./sqrt(3.)*(atan((2.-a)/a/sqrt(3.))+pi/6.)) 
  enddo

  if(NBI_flag) then
    ! nbi_slowing down density in 10**19 1/m**3 has been read

    do i = 1,nr
      ln_lambda = 17.
      tau_ee(i) = 1.088E-3*(T_e_rho(i))**1.5/n_e_rho(i)/ln_lambda ! in sec
      !tau_s needs be used in the transport part
      tau_s_rho(i) = 1836.*m_alpha/z_alpha**2*tau_ee(i)
      
      S0_rho(i)=n_alpha_rho(i)/tau_s_rho(i)/I2(i)  !in [10**19 1/m**3]
    enddo

  ! nbi equivalent Maxwellian temperature in keV has been read
  else
  ! alpha_slowing down density in 10**19 1/m**3

    n_alpha_rho(:) = 0.
    do i = 1,nr
      S0_rho(i)=2.5E-6*(n_i_rho(i))**2*(T_i_rho(i))**2 !in [10**19 1/m**3]/sec

      ln_lambda = 17.
      tau_ee(i) = 1.088E-3*(T_e_rho(i))**1.5/n_e_rho(i)/ln_lambda ! in sec
      !tau_s needs be used in the transport part
      tau_s_rho(i) = 1836.*m_alpha/z_alpha**2*tau_ee(i)
      
      n_alpha_rho(i)=S0_rho(i)*tau_s_rho(i)*I2(i)  !in [10**19 1/m**3]
    enddo


  endif

  ! alpha equivalent Maxwellian temperature in keV
  T_alpha_equiv_rho(:) = 0.
  do i = 1,nr
    T_alpha_equiv_rho(i)=2./3.*I4(i)/I2(i)*E_alpha*10.**3
  enddo

!calculate smeared alpha density  n_alpha_ave_rho(i)

  rho_hat_p(:)=rho_hat(:)
  do i=2,nr
    rho_t=1.02E2*(1./2.)*sqrt(4.)*sqrt(2.*T_alpha_equiv_rho(i)*1000.)/(B_t*10000.)/(rmin*100.)
    rho_p=rho_t*rmaj_rho(i)/rmin_rho(i)*q_rho(i)
 !   rho_p always greater than rho_t
    delta=rho_p
    if(rho_p.gt.rho_hat(i)) delta=rho_hat(i)  !banana half width less than radius not allowed
    if(delta.lt.rho_t) delta=rho_t   ! there is aways a minimum smearing over rho_t
    bannana_rho(i)=delta
    Numer=0.
    Denom=0.
    do ii=1,nr
      Denom=Denom+exp(-(rho_hat(i)-rho_hat_p(ii))**2/delta**2)
      Numer=Numer+exp(-(rho_hat(i)-rho_hat_p(ii))**2/delta**2)*n_alpha_rho(ii)
    enddo
    n_alpha_ave_rho(i)=Numer/Denom
  enddo
  n_alpha_ave_rho(1)=n_alpha_ave_rho(2)


! compute the increase in energy flux over simple 3/2 T_alpha convection

  n_en = 1000
!   n_en_max = 200
  n_en_max = 1000
  !print *, 'n_en'
 
  !print *, 'Convec_factor'
!  do i=1,nr
    i=5
!    i=25
!    i=45
    print *, 'I4/I2=',I4(i)/I2(i)
    print *, 'T_e_rho=',T_e_rho(i)
    print *, 'E_c_hat_rho=',E_c_hat_rho(i)
    Numer=0.0
    Denom=0.0
    do i_en = 1,n_en_max
      En = (real(i_en)-0.5)/real(n_en)
      EnoTe = En*E_alpha*10**3/T_e_rho(i)
      G_D=exp(-8.14E-5*EnoTe**4+3.77E-3*EnoTe**3-0.0553*EnoTe**2+0.036*EnoTe+0.45)
 !      if(EnoTe .le. 2.7) G_D =1.25  !Angioni-Peters Fig 2 Eq. 32  
 !test   passed G_D = 1.0 test
 !      G_D=1.0
       if(EnoTe .le. 2.7) then
         G_D=1.0
         if(EnoTe .gt. 0.3) G_D = 0.25*(EnoTe-0.3)/(2.7-0.3)+1.0
       endif
      Numer=Numer+1./real(n_en)*En*G_D*En**0.5/(E_c_hat_rho(i)**1.5+En**1.5)
      Denom=Denom+1./real(n_en)*G_D*En**0.5/(E_c_hat_rho(i)**1.5+En**1.5)
    !  print *, i_en,En,EnoTe,G_D,Numer,Denom,Numer/Denom
    enddo
    Convec_factor = Numer/Denom/(I4(i)/I2(i))
    print *, i,Numer,Denom,Convec_factor
!   enddo


end subroutine EPtran_comp_alpha_slowing
