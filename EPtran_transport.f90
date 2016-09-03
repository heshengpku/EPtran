
!---------------------------------------------------------
! EPtran_transport.f90
!
! PURPOSE:
!  Reads all transport control parameters from : EPtran_transport_control
!  compute alpha slowing down distribution transport
!
! 6.26.15
!---------------------------------------------------------

subroutine EPtran_transport

    use EPtran_use_parameter
!CONTROL    !dt,ntstep,i_new_start,i_flux_r,i_flux_E,i_Model,i_AE,
!           !i_pinch,i_thermal_pinch_off,i_BC_r1,delta1,D_factor,A_factor
!INPUTS     !nr,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
!           !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
!           !Zeff,M_DT,E_alpha,Fpeak_T,Fpeak_n,Fpeak_ei,NBI_flag
!OUTPUTS    !beta_ped_percent,T_ped,beta_N_glob,arho,
!           !rho_hat(i),q_rho(i),Rmaj_rho(i),rmin_rho(i),
!           !kappa_rho(i),delta_rho(i),T_i_rho(i),T_e_rho(i),
!           !n_i_rho(i),n_e_rho(i),n_alpha_rho(i),T_alpha_equiv_rho(i),
!           !E_c_hat_rho(i)
!OTHER      !pi,n_alpha_ave_rho(i),S0_rho(i),bannana_rho(i),tau_s_rho(i)
!           !V_prime_rho(i),chi_eff(i),D_bkg_Angioni(i),C_p_alpha(i),
!           !rg_n_alpha_th_rho(i)

    use EPtran_use_transport  !nEa,nE,dE,E_hat(j),ny,dy,y(k)
                              !Int_E,Int_y,IntE1y,IntE2y,rdiff


  !--------------------------------------
  implicit none

  real :: err = 1.e-8 !the smallest real
  real,dimension(ny) :: check_y
  real,dimension(nE) :: check_E
  real :: E_c_middle

  real :: Q_fus
  real, dimension(nr) :: V_prime_rho
  real, dimension(nr) :: chi_eff
  real, dimension(nr) :: D_bkg_Angioni  !diffusive part
  real, dimension(nr) :: C_p_alpha      !pinch part
  real, dimension(nr) :: D_bkg_Angioni_He
  real, dimension(nr) :: C_p_He
  real, dimension(nr) :: rg_n_alpha_th_rho  !AE density gradient threshold

  integer :: i,j,k,n
  integer :: ku,kd !k+1,k-1
  integer :: nt  !loop for time
  integer :: nrk !loop for Runge-Kutta
  integer :: nsav,ndt,ndata,nstart !for data save
  real :: nstart_r
  !character(len=8) :: str

  integer :: i_out,j_out

  real :: dr !radius space step
  real :: tstep

  logical :: i_isotropic ! .true.=isotropic source, .false.=beam-like source
  logical :: i_scattering ! .true.=pitch angle scattering, .false.=no
  logical :: i_start_null ! .true.=start from null distribution, .false.=from S. D.
  integer :: N_source !the injected energy component of the source

  real :: lambda0,Delta_lambda
  real :: Int_F_bpa
  real,dimension(ny) :: F_bpa !unnormalized F_bpa(lambda)

  real,dimension(nE,nr) :: nu_d

  real :: f_error_max

  real :: TeoEalpha   ! T_e_rho(i)/E_alpha
  real :: EoTe         ! E/T_e_rho(i)
  real,dimension(nr) :: c,G_D_s,c_G_D_s    !c = c_G_D_s/G_D_s
  real,dimension(nE,nr) :: G_D

  !real :: k_E
  real,dimension(nr) :: k_p
  !real,dimension(nE,nr) :: K_s  !K_s(E_hat)
  real,dimension(nE,nr) :: H_s  !H_s(E_hat)

  real :: sigma_es_p,sigma_es_p_FLR,sigma_es_t
  real :: D_es_p,D_es_p_FLR,D_es_t
  real :: epsilon_t
  real :: y_tp
  logical :: Pueschel_em_flag
  real :: sigma_em_p,sigma_em_p_FLR,sigma_em_t
  real :: D_em_p,D_em_p_FLR,D_em_t
  real :: beta_square

  logical :: i_orbit ! .true.= Orbit Width effects
  integer :: irBC_inner,irBC_outer,ir_half,ir_i,ii
  real :: ir,ir_p,lambda_tp,Orbit_factor,D_orb_add
  real,dimension(nr) :: B_theta,omega_ctheta
  real,dimension(nr) :: Delta_orb_ave_hat,rg_n_alpha_th_rho_bk
  real,dimension(nE) :: v_par0
  real,dimension(nE,nr) :: Delta_orb_hat

  real :: InPar !the total fusion particle
  real :: InPower !the total fusion power
  real,dimension(nr) :: source_rho !the alpha particle source in 10**19/sec
  real,dimension(nr) :: flux_source_rho !the source flux in 10**19/m**2/sec
  real,dimension(nr) :: flow_source_rho !the source flow in 10**19/sec
  real,dimension(nr) :: sink_rho !the alpha particle sink and He source in 10**19/sec
  real,dimension(nr) :: sinke_rho !the slowing down energy

  real,dimension(nr) :: p_alpha_rho   !slowing down pressure profile
  real,dimension(nr) :: n_tran_alpha_rho  !transported denstiy profile
  real,dimension(nr) :: T_tran_alpha_rho  !transported temperature profile
  real,dimension(nr) :: p_tran_alpha_rho  !transported pressure profile
  real,dimension(nr) :: rg_n_tran_alpha_rho  !density gradient
  real,dimension(nr) :: rg_T_tran_alpha_rho  !temperature gradient
  real,dimension(nr) :: rg_p_tran_alpha_rho  !pressure gradient

  !real,dimension(nr) :: flux
  real,dimension(nr) :: flow_rho        !the particle flow
  real,dimension(nr) :: flow_check_rho  !the particle check flow
  real,dimension(nr) :: flow_energy_rho !the energy flow

  real,dimension(nr) :: n_s_He_rho     !the untransported He ash density
  real,dimension(nr) :: n_tran_He_rho  !the transported He ash density
  real,dimension(nr) :: flux_He_rho        !the He ash particle flux
  real,dimension(nr) :: flow_He_rho        !the He ash particle flow

  real,dimension(nE,nr) :: F_s
  real,dimension(ny,nE,nr) :: f_alpha

  real,dimension(ny,nE,nr) :: D_rr,A,C_rE,C_EE

  integer :: threshold_flag !0 for density threshold, 1 for pressure threshold
  real,dimension(ny,nE,nr) :: D_ITG
  real :: D_AE
  real :: dn_th

  real,dimension(nr) :: rg_p_alpha_th_rho
  real :: dp_th

  real :: C_R
  real,dimension(nr) :: F_cor
  real,dimension(nr) :: L_s_n_alpha_rho,L_s_T_alpha_rho,L_tran_n_alpha_rho
  integer :: i_threshold

  logical :: i_energy_AE !.true. = energy-dependent critical gradient model

  real,dimension(nE,nr) :: G_AE
  real :: eps_hat
  real :: eps_hat_0
  real :: Delta_eps_hat

  real,dimension(ny,nE,nr) :: Gamma_r
  real,dimension(ny,nE,nr) :: Gamma_E

  real,dimension(ny,nE,nr) :: f_alpha_bk

  real,dimension(ny,nE,nr) :: term_csd  !term of classical slowing down
  real,dimension(ny,nE,nr) :: term_s  !term of alpha source

  integer,dimension(:),allocatable :: j_NBI !the injected energy of NBI

  real,dimension(ny,nE,nr) :: term_flux_r  !term of flux_r
  real,dimension(ny,nE,nr) :: term_flux_E  !term of flux_E
  real,dimension(ny,nE,nr) :: term_flux_y  !term of flux_y
  real,dimension(ny,nE,nr) :: term_pas  !term of pitch angle scattering

  logical :: i_distribution_out !.true.=write out the distributions in .dat files
  logical :: iexist
  !--------------------------------------------------------------
  !Control Parameters
  !--------------------------------------------------------------

  !!5.11.2015
  i_isotropic = .true. !isotropic source
  !i_isotropic = .false. !beam-like source for NBI
  !i_scattering = .true. !pitch angle scattering
  i_scattering = .false. !no pitch angle scattering
  !!5.29.2015
  i_energy_AE = .false. !energy independent
  !i_energy_AE = .true. !energy dependent
  !!9.30.2015
  i_orbit = .false. !no Orbit Width effect
  !i_orbit = .true. !Orbit Width effect
  !!6.26.2015
  i_start_null = .true. !start from null distribution
  !i_start_null = .false. !start from slowing down distribution
  !i_distribution_out = .true. !.true.=write out the distributions in .dat files
  i_distribution_out = .false.
  threshold_flag = 0 !AE density gradient threshold
  !threshold_flag = 1 !AE pressure gradient threshold
  N_source = 1 !default, one injected energy component source
  !N_source = 3 !DIII-D, three injected energy component source

  inquire(file='EPtran_transport_control',exist=iexist)
  if(iexist) then
    open(unit=1,file='EPtran_transport_control',status='old')

    read(1,*) i_new_start
    read(1,*) dt
    read(1,*) ntstep
    read(1,*) i_flux_r
    read(1,*) i_flux_E
    read(1,*) i_Model
    read(1,*) i_AE
    read(1,*) i_pinch
    read(1,*) i_thermal_pinch_off
    read(1,*) i_BC_r1
    read(1,*) delta1
    read(1,*) D_factor
    read(1,*) A_factor

    close(1)
  else
    print *, 'EPtran_transport_control file is not found'
    stop
  endif

  !--------------------------------------------------------------
  !Grids setting and check
  !--------------------------------------------------------------

  dr = rmin/real(nr-1)

! calculate the E_hat grid
  do j = 1,nE
    E_hat(j) = (j-1.)/(nE-nEa-1.)
  enddo

! calculate the lambda grid
  do k = 1,ny
    y(k) = k/(ny+1.)
  enddo

! check the integral function
  check_y = 1-y
  print *,'Integral y =',Int_y(check_y)
  print *,'Theory = 1/2'

  do j = 1,nE
    check_E(j) = 1./2./pi
    if(j .gt. nE-nEa) then
      check_E(j) = 0.
    endif
  enddo
  print *,'Integral E =',Int_E(check_E,nE)
  print *,'Theory = 2/3'

  !--------------------------------------------------------------
  !source setting, no transport
  !--------------------------------------------------------------
  if(i_thermal_pinch_off) then !turn off the thermal pinch, E_c(r) = E_c(r/a=0.5)
    E_c_middle = E_c_hat_rho((nr+1)/2)
    E_c_hat_rho(:) = E_c_middle
  endif

  ! static alpha slowing down distribution
  F_s(:,:) = 0.
  do i = 1,nr
    do j = 1,nE-nEa
      F_s(j,i) = S0_rho(i)*tau_s_rho(i)/4/pi/ &
                (E_c_hat_rho(i)**1.5+E_hat(j)**1.5)
    enddo
  enddo

  do i = 1,nr
    do j = 1,nE
      nu_d(j,i) = 1./tau_s_rho(i)*E_c_hat_rho(i)**1.5/ &
                  2./(0.17**3+E_hat(j)**1.5)
    enddo
  enddo

  !Pitch angle lambda-dependent source
  lambda0 = 0.5
  Delta_lambda = 0.2
  F_bpa(:) = exp(-(1-y(:)**2-lambda0)**2/Delta_lambda**2)
  Int_F_bpa = Int_y(F_bpa)

  !set V_prime_rho(i) in m**2
  !since  V' enters as 1/V' d[V' D dn/dr]/dr
  !size of V' or accuracy of V' not vey important
  V_prime_rho(:) = 2.*pi*kappa_rho(:)*rmin_rho(:)*2.*pi*Rmaj_rho(:)

  !The source alpha particles in 10**19/sec'
  source_rho(:) = S0_rho(:)*V_prime_rho(:)*dr

  !The source flow in 10**19/sec
  flow_source_rho(1) = 0.
  do i = 2,nr
    flow_source_rho(i) = flow_source_rho(i-1) + &
                         0.5*(source_rho(i)+source_rho(i-1))
  enddo

  !The source flux in 10**19/m**2/sec
  flux_source_rho(:) = flow_source_rho(:)/V_prime_rho(:)
  flux_source_rho(1) = 0.

  Q_fus = 10.  !default
  !!!  Q_fus = 20.  !for the 2x baseline case
  if(NBI_flag) Q_fus = 1000000. !for DIIID NBI test with effective E_alpha

  !12.4.2014 fixed T_e NOTnec= T_i and n_e NOTnec= n_i
  !plasma chi_eff from  alpha energy source flow times (1+5./Q_fus)
  !T_e NOTnec= T_i and n_e NOTnec= n_i assumed:  chi_eff in m**2/sec
  chi_eff(:) = (1.+5./Q_fus)*flux_source_rho(:)*(E_alpha*1000.)/ &
               (0.5*n_i_rho(:)*rdiff(T_i_rho(:))/dr + &
                0.5*n_e_rho(:)*rdiff(T_e_rho(:))/dr)
  chi_eff(1) = chi_eff(2)

  !--------------------------------------------------------------
  !pre-Transport setting: diffusivity
  !--------------------------------------------------------------

  !alpha ITG/TEM diffusivity from Angioni formulas normed to chi_eff
  D_bkg_Angioni(:) = chi_eff(:)* &
                     (0.02+4.5*(T_e_rho(:)/(E_alpha*1000.)) &
                     +8.0*(T_e_rho(:)/(E_alpha*1000.))**2 &
                     +350.*(T_e_rho(:)/(E_alpha*1000.))**3)

  !Angioni pinch coef C_p_alpha
  !Angi_exp = -1, corrected
  C_p_alpha(:) = (3./2.)*Rmaj_rho(:)* &
                 rdiff(T_e_rho(:))/dr/T_e_rho(:)* &
                 (1./(1.+1./E_c_hat_rho(:)**(-1.5))/ &
                 log(1.+1./E_c_hat_rho(:)**1.5)-1.)
  C_p_alpha(1) = 0.

  !He ITG/TEM diffusivity from Angioni formulas normed to chi_eff
  D_bkg_Angioni_He(:) = chi_eff(:)

  !Angioni pinch coef C_p_He = -2.0
  !Small Angioni thermal expulsion +C_T_He*n_He/L_Ti set to 0
  C_p_He(:) = -2.0
  C_p_He(1) = 0.

  if(NBI_flag) then!i_threshold = 9  !NBI DIIID caseD   SD-beam-like
    i_threshold = 0
    select case(i_threshold)
    case(1) !new GYRO gamma_AE+ITG > gamma_ITG/TEM
      do i = 1,nr
        if(rho_hat(i) .ge. 0.45) &
          rg_n_alpha_th_rho(i) = 0.0975*(1.0+(rho_hat(i)-0.45)**2/(0.1235)**2)
        if(rho_hat(i) .lt. 0.45) &
          rg_n_alpha_th_rho(i) = 0.0975*(1.0+(rho_hat(i)-0.45)**2/(0.1482)**2)
      enddo
    case(2) !TGLF n = 3, gamma_AE > 0
      do i = 1,nr
        if(rho_hat(i) .ge. 0.45) &
          rg_n_alpha_th_rho(i) = 0.0487*(1.0+(rho_hat(i)-0.45)**2/(0.0975)**2)
        if(rho_hat(i) .lt. 0.45) &
          rg_n_alpha_th_rho(i) = 0.0487*(1.0+(rho_hat(i)-0.45)**2/(0.0817)**2)
      enddo
    case(3) !TGLF worst-n, gamma_AE > 0
      do i = 1,nr
        if(rho_hat(i) .ge. 0.45) &
          rg_n_alpha_th_rho(i) = 0.0456*(1.0+(rho_hat(i)-0.45)**2/(0.1049)**2)
        if(rho_hat(i) .lt. 0.45) &
          rg_n_alpha_th_rho(i) = 0.0456*(1.0+(rho_hat(i)-0.45)**2/(0.1372)**2)
      enddo
    case default !original GYRO gamma_AE > gamma_ITG/TEM
      do i = 1,nr
        if(rho_hat(i) .ge. 0.5) &
          rg_n_alpha_th_rho(i) = 0.075*(1.0+(rho_hat(i)-0.5)**2/(0.0903)**2)
        if(rho_hat(i) .lt. 0.5) &
          rg_n_alpha_th_rho(i) = 0.075*(1.0+(rho_hat(i)-0.5)**2/(0.1035)**2)
      enddo
    end select
  else  !i_threshold = 7
    do i = 1,nr
      if(rho_hat(i) .ge. 0.4) &
        rg_n_alpha_th_rho(i) = 0.017*(1.0+(rho_hat(i)-0.4)**2/(0.2287)**2)
      if(rho_hat(i) .lt. 0.4) &
        rg_n_alpha_th_rho(i) = 0.017*(1.0+(rho_hat(i)-0.4)**2/(0.125)**2)
    enddo
  endif

  !5.29.2015, develop an energy dependent critical gradient stiff diffusion model
  if(i_energy_AE) then
    Delta_eps_hat = 1./3.
    eps_hat_0 = 3.
    !eps_hat_0 = 0.
    !eps_hat_0 = 1.
    do i = 1,nr
      do j = 1,nE
        eps_hat = E_hat(j)*E_alpha*1000./T_alpha_equiv_rho(i)
        G_AE(j,i) = Delta_eps_hat**2/((eps_hat-eps_hat_0)**2+Delta_eps_hat**2)
      enddo
    enddo
  else
    G_AE(:,:) = 1.
  endif

! calculate the D matrix
  select case(i_Model)
  case(1) !energy independent Angioni model D(r)
    do i = 1,nr
      do j = 1,nE
        do k = 1,ny
          D_rr(k,j,i) = D_factor*D_bkg_Angioni(i)
        enddo
      enddo
    enddo

    if(i_pinch) then
      do i = 1,nr
        do j = 1,nE
          do k = 1,ny
            C_rE(k,j,i) = A_factor*D_rr(k,j,i)*C_p_alpha(i)/Rmaj_rho(i)
          enddo
        enddo
      enddo
    endif

  case(2) !energy dependent Angioni model D(r,E)
    do i = 1,nr

      TeoEalpha = T_e_rho(i)/E_alpha/1.e3
      c_G_D_s(i) = 0.02 + 4.5*TeoEalpha + 8.*TeoEalpha**2 + 350.*TeoEalpha**3

      do j = 1,nE
        EoTe = E_hat(j)/TeoEalpha
        if(EoTe .le. 2.7) then
          G_D(j,i) = 1.25
        else if(EoTe .gt. 33.05) then
          G_D(j,i) = 0.
        else
          G_D(j,i) = exp(-8.14e-5*EoTe**4 + 3.77e-3*EoTe**3 &
                     - 0.0553*EoTe**2 + 0.036*EoTe + 0.45)
        endif
      enddo

      G_D_s(i) = Int_E(G_D(:,i)*F_s(:,i),nE)/n_alpha_rho(i)

      !k_E = 1.

      !K_s(i,1) = 1.5*(1./log((1./E_c_hat_rho(i)**1.5+1)*(1+E_c_hat_rho(i)**1.5))-1)
      !do j = 2,nE-nEa
      !  K_s(j,i) = 1.5*(1./log((1./E_c_hat_rho(i)**1.5+1)*(1+E_c_hat_rho(i)**1.5))-&
      !        E_c_hat_rho(i)**1.5/(E_c_hat_rho(i)**1.5+E_hat(j)**1.5))
      !enddo
      !do j = nE-nEa+1,nE
      !  K_s(j,i) = 0.
      !enddo

      if(G_D_s(i) .lt. err) then
        c(i) = 0
      else
        c(i) = c_G_D_s(i)/G_D_s(i)
      endif

      do j = 1,nE
        do k = 1,ny
          D_rr(k,j,i) = D_factor*c(i)*chi_eff(i)*G_D(j,i)
        enddo
      enddo
    enddo

    if(i_pinch)then
      do i = 1,nr

        H_s(1,i) = 0.
        do j = 2,nE-nEa
          H_s(j,i) = 1.5*(E_hat(j)**0.5/(E_c_hat_rho(i)**1.5+E_hat(j)**1.5))
        enddo
        do j = nE-nEa+1,nE
          H_s(j,i) = 0.
        enddo

        k_p(i) = 0.2*Int_E(G_D(:,i),nE)/Int_E(G_D(:,i)*H_s(:,i),nE)

        do j = 1,nE
          do k = 1,ny
            A(k,j,i) = A_factor*rmin/T_alpha_equiv_rho(i)*1.e3*E_alpha/ &
                       Rmaj_rho(i)*k_p(i)
          enddo
        enddo
      enddo
    endif

  case(3) !read dep(nE,ny,nr) from files, for DEP model
    open(unit=15,file='dep.dat',status='old',form='Unformatted')
    read(15) D_rr
    close(15)

    D_rr(:,:,:) = D_factor*D_rr(:,:,:)

    if(i_pinch) then
      !A(:,:,:) = 0. 6.24.2015
      open(unit=15,file='aep.dat',status='old',form='Unformatted')
      read(15) A
      close(15)

      A(:,:,:) = A_factor*A(:,:,:)
    endif

  case(4) !Pueschel model
    Pueschel_em_flag = .false. !only electrostatic part
    !Pueschel_em_flag = .true. !with electromagnetic counterparts
    beta_square = 0.5**2 !(beta/beta_crit)**2, 5.13.2015

    sigma_es_p = 0.292
    sigma_es_p_FLR = 0.422
    sigma_es_t = 0.527

    sigma_em_p = 0.047
    sigma_em_p_FLR = 0.125
    sigma_em_t = 0.157

    do i = 1,nr
      do j = 2,nE
        do k = 1,ny
          EoTe = E_hat(j)*E_alpha*1.e3/T_e_rho(i)
          epsilon_t = rmin_rho(i)/Rmaj_rho(i)

          D_es_p = sigma_es_p*chi_eff(i)/y(k)**2*EoTe**(-1)
          D_es_p_FLR = sigma_es_p_FLR*chi_eff(i)&
                       /(y(k)**2*(1.-y(k)**2)**0.5)*EoTe**(-3./2.)
          D_es_t = sigma_es_t*chi_eff(i)*epsilon_t**0.5&
                   /(y(k)*(1.-y(k)**2))*EoTe**(-3./2.)

          D_em_p = sigma_em_p*chi_eff(i)*beta_square
          D_em_p_FLR = sigma_em_p_FLR*chi_eff(i)*beta_square &
                       /(1.-y(k)**2)**0.5*EoTe**(-1./2.)
          D_em_t = sigma_em_t*chi_eff(i)*beta_square*epsilon_t**0.5 &
                   *y(k)/(1.-y(k)**2)*EoTe**(-1./2.)
          y_tp = sqrt(2.*epsilon_t/(1.+epsilon_t))
          if(Pueschel_em_flag) then
            if(y(k) .le. y_tp) then !trapped particles
              D_rr(k,j,i) = D_es_t + D_em_t
            else                    !passing particles
              D_rr(k,j,i) = min(D_es_p,D_es_p_FLR) + min(D_em_p,D_em_p_FLR)
            endif
          else
            if(y(k) .le. y_tp) then !trapped particles
              D_rr(k,j,i) = D_es_t
            else                    !passing particles
              D_rr(k,j,i) = min(D_es_p,D_es_p_FLR)
            endif
          endif
        enddo
      enddo
    enddo

    j = 1 !!E_hat = 0
    do i = 1,nr
      do k = 1,ny
        D_rr(k,j,i) = D_rr(k,j+1,i)
      enddo
    enddo

    D_rr = D_factor*D_rr

    if(i_pinch) then
      A = 0.1*A_factor   !!assuming a constant, 6.25.2015
    endif

  case default !const D matrix
    do i = 1,nr
      do j = 1,nE
        do k = 1,ny
          D_rr(k,j,i) = 0.1*D_factor*G_AE(j,i)    !const diffusion
        enddo
      enddo
    enddo

    if(i_pinch) then
      A = 0.1*A_factor
    endif

  end select

  ! !! Write some imformations of D_rr
  ! open(unit=7,file='EPtran_Dep.out',status='replace')
  ! write(7,*) '--------------------------------------------------------------'
  ! write(7,*) 'D_rr(r) vs r, E and lambda have been integrated with F_s'
  ! write(7,*) '--------------------------------------------------------------'
  ! do i = 1,nr
  !   if(i_isotropic) then
  !     do j = 1,nE
  !       check_E(j) = Int_y(D_rr(:,j,i))
  !     enddo
  !   else
  !     do j = 1,nE
  !       check_E(j) = Int_y(D_rr(:,j,i)*F_bpa(:))/Int_F_bpa
  !     enddo
  !   endif
  !   write(7,10) Int_E(check_E(:)*F_s(:,i),nE)/Int_E(F_s(:,i),nE)
  ! enddo

  ! write(7,*) '--------------------------------------------------------------'
  ! write(7,*) 'D_rr(r,E) vs E at different radius, lambda has been integrated'
  ! write(7,*) '--------------------------------------------------------------'
  ! do i_out = 0.1,0.9,0.2  ! r = 0.1, 0.3, 0.5, 0.7, 0.9
  !   i = nint(i_out*(nr-1)+1)
  !   write(7,*) 'r/a = ',i_out
  !   write(7,*) '--------------------------------------------------------------'
  !   if(i_isotropic) then
  !     do j = 1,nE
  !       write(7,10) Int_y(D_rr(:,j,i))
  !     enddo
  !   else
  !     do j = 1,nE
  !       write(7,10) Int_y(D_rr(:,j,i)*F_bpa(:))/Int_F_bpa
  !     enddo
  !   endif
  ! enddo
  ! close(7)

  if(.not. i_pinch) then  !no pinch
    C_rE = 0.
  endif

  !!9.30.2015, add the Orbit Width effect
  if(i_orbit) then
    !lambda = lambda0 = 0.5
    v_par0(:) = sqrt(2.*E_hat(:)*E_alpha/m_alpha/931.*(1.-lambda0))*3.e8

    !1.2.2016 B_p = B_unit*(r/Rq)|grad(r)| ~ B_0*kappa*(r/Rq)/(1+dR_0/dr)
    B_theta(:) = B_t*kappa_rho(:)*rmin_rho(:)/Rmaj_rho(:)/q_rho(:)/(1.-delR0oa)

    !1.14.2016, omega_ctheta = B_theta*Ze/m
    omega_ctheta(:) = B_theta(:)*z_alpha/m_alpha/931./1.e6*(3.e8)**2

    Orbit_factor = 1.
    !Orbit_factor = 0.5  !!10.27.2015
    do i = 2,nr
      !10.8.2015, passing and trapped particles
      epsilon_t = rmin_rho(i)/Rmaj_rho(i)
      lambda_tp = (1.-epsilon_t)/(1.+epsilon_t)
      do j = 1,nE
        !1.15.2016
        if(NBI_flag) then !For beam particles
          !10.21.2015,1/2 factor for both passing and trapped particles
          Delta_orb_hat(j,i) = Orbit_factor*0.5*v_par0(j)/omega_ctheta(i)/rmin

          !if(lambda0 .gt. lambda_tp) then !trapped
          !  Delta_orb_hat(j,i) = v_par0(j)/omega_ctheta(i)/rmin
          !else !passing
          !  Delta_orb_hat(j,i) = 0.5*v_par0(j)/omega_ctheta(i)/rmin
          !endif
        else              !For isotropic alpha
          Delta_orb_hat(j,i) = Orbit_factor*v_par0(j)/sqrt(1.-lambda0)&
                               *0.5*(1.-0.5*lambda_tp)/omega_ctheta(i)/rmin
        endif
      enddo
    enddo

    do i = 2,nr
      Delta_orb_ave_hat(i) = Int_E(F_s(:,i)*Delta_orb_hat(:,i),nE)/ &
                             Int_E(F_s(:,i),nE)
    enddo
    Delta_orb_ave_hat(1) = Delta_orb_ave_hat(2)

    do i = 1,nr
      if(Delta_orb_ave_hat(i) .lt. rho_hat(i)) exit;
      irBC_inner = i;
    enddo

    do i = nr,1,-1;
      if(Delta_orb_ave_hat(i) + rho_hat(i) .lt. 1.) exit;
      irBC_outer = i;
    enddo

    !!9.30.2015, broaden the critical gradient profile
    rg_n_alpha_th_rho_bk(:) = rg_n_alpha_th_rho(:)
    ir_half = (nr+1)/2
    !do i = 1,ir_half-1
    do i = irBC_inner,ir_half
      !ir = i + Delta_orb_ave_hat(i)*real(nr-1)
      ir = i - Delta_orb_ave_hat(i)*real(nr-1)
      if(ir .lt. 0) ir = 0.01
      ir_i = floor(ir)
      !ir_p = ir-ir_i
      do ii = ir_i+1,i-1  !!new broaden method
        rg_n_alpha_th_rho(ii) = rg_n_alpha_th_rho_bk(i)
      enddo
      !if(ir .ge. ir_half) then
      !  rg_n_alpha_th_rho(i) = rg_n_alpha_th_rho_bk(ir_half)
      !else
      !  rg_n_alpha_th_rho(i) = (1-ir_p)*rg_n_alpha_th_rho_bk(ir_i) +&
      !                         ir_p*rg_n_alpha_th_rho_bk(ir_i+1)
      !endif
    enddo

    do i = ir_half+1,nr
      ir = i - Delta_orb_ave_hat(i)*real(nr-1)
      ir_i = floor(ir)
      ir_p = ir-ir_i
      if(ir .le. ir_half) then
        rg_n_alpha_th_rho(i) = rg_n_alpha_th_rho_bk(ir_half)
      else
        rg_n_alpha_th_rho(i) = (1-ir_p)*rg_n_alpha_th_rho_bk(ir_i) + &
                               ir_p*rg_n_alpha_th_rho_bk(ir_i+1)
      endif
    enddo

    !!add a large D to simulate the new inner BC for orbit width effect
    !D_orb_add = 10.
    !D_orb_add = 20. !10.22.2015
    !print *,'D_orb_add =',D_orb_add
    !do i = 1,irBC_inner
    !  !D_rr(:,:,i) = D_rr(:,:,i) + D_orb_add
    !  D_rr(:,:,i) = D_orb_add
    !enddo
    !do i = irBC_outer,nr
    !  !D_rr(:,:,i) = D_rr(:,:,i) + D_orb_add
    !  D_rr(:,:,i) = D_orb_add
    !enddo
  else

    irBC_inner = 1
    irBC_outer = nr

  endif

  C_R = 2.1 !From TGLF simulations

  if(i_AE) then !D_ITG = D_rr for AE on
    D_ITG = D_rr

    L_s_n_alpha_rho = n_alpha_rho/(rdiff(n_alpha_rho)/dr) !08.29.16
    L_s_T_alpha_rho = T_alpha_equiv_rho/(rdiff(T_alpha_equiv_rho)/dr)
    L_s_n_alpha_rho(1) = L_s_n_alpha_rho(2) !avoid NaN at r/a = 0
    L_s_T_alpha_rho(1) = L_s_T_alpha_rho(2)

    if(threshold_flag .eq. 1) then
    !calculate the pressure gradient threshold,7.13.2014
      !(-dp_th/dr) = T_s*(-dn_th/dr)*(1+(-dlnT_s/dr)/(-dlnn_s/dr))
      rg_p_alpha_th_rho(:) = T_alpha_equiv_rho(:)*rg_n_alpha_th_rho(:)* &
                             0.16022*(1. + L_s_n_alpha_rho(:)/L_s_T_alpha_rho(:))
      !rg_p_alpha_th_rho(1) = T_alpha_equiv_rho(1)*rg_n_alpha_th_rho(1)*0.16022

    endif

    if(NBI_flag) then
      if(i_energy_AE) then
        D_AE = 25.
        !D_AE = 150.
      else
        D_AE = 10.
        !D_AE = 20.  !must be the same with Alpha_transport.f90
        !D_AE = 25.
      endif
    else
      D_AE = 0.3  !default
    endif

  endif

  if(N_source .eq. 3) then
  !the full, half, and third components of the nominal 80keV beams in DIII-D
    allocate(j_NBI(N_source))
    j_NBI(1) = nint((nE-nEa-1.)/3.)+1
    j_NBI(2) = nint((nE-nEa-1.)/2.)+1
    j_NBI(3) = nint((nE-nEa-1.)/1.)+1
  endif

  !--------------------------------------------------------------
  !Transport main code
  !--------------------------------------------------------------

  term_csd(:,:,:) = 0.
  term_s(:,:,:) = 0.
  term_pas(:,:,:) = 0.
  Gamma_r(:,:,:) = 0.
  Gamma_E(:,:,:) = 0.
  term_flux_r(:,:,:) = 0.
  term_flux_E(:,:,:) = 0.
  term_flux_y(:,:,:) = 0.

  !set the start f_alpha
  if(i_new_start) then  !new start
    nstart = 0
    if(i_start_null) then  !default from null distribution
      f_alpha = 0.
    elseif(i_isotropic) then  !from slowing down distribution
      do i = 1,nr
        do j = 1,nE
          do k = 1,ny
            f_alpha(k,j,i) = F_s(j,i)
          enddo
        enddo
      enddo
    else
      do i = 1,nr
        do j = 1,nE
          do k = 1,ny
            f_alpha(k,j,i) = F_s(j,i)*F_bpa(k)/Int_F_bpa
          enddo
        enddo
      enddo
    endif
  else  !restart
    open(unit=11,file='f_bac.dat',status='old',form='Unformatted')
    read(11) nstart_r
    nstart = nint(nstart_r)
    read(11) f_alpha
    close(11)
  endif

  nsav = 10
  ndt = (ntstep-nstart)/nsav
  ndata = 0

  !the main loop for time developing
  do nt = nstart+1,ntstep
    do nrk = 1,2

      ! 1st step of Runge-Kutta method
      if(nrk .eq. 1) then
        tstep = 0.5*dt
        !advance f_alpha_bk
        !$OMP PARALLEL DO PRIVATE(k,j,i)
        do i = 1,nr
          do j = 1,nE
            do k = 1,ny
              f_alpha_bk(k,j,i) = f_alpha(k,j,i)
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      ! 2nd step of Runge-Kutta method
      else
        tstep = dt
      endif


      !calculate the classical slowing down and source terms
      !$OMP PARALLEL DO PRIVATE(i,j,k,kd,ku)
      do i = 1,nr

        j = 1
        do k = 1,ny
          !notice here have 1/sqrt(E), 6.13.2015
          term_csd(k,j,i) = 2./tau_s_rho(i)/(dE/100.)**0.5* &
                  ((E_hat(j+1)**1.5+E_c_hat_rho(i)**1.5)*f_alpha(k,j+1,i)- &
                  (E_hat(j)**1.5+E_c_hat_rho(i)**1.5)*f_alpha(k,j,i))/dE
          !term_csd(k,j,i) = 2./tau_s_rho(i)/(E_hat(j)**0.5+E_hat(j+1)**0.5)*2.* &
          !        ((E_hat(j+1)**1.5+E_c_hat_rho(i)**1.5)*f_alpha(k,j+1,i)- &
          !        (E_hat(j)**1.5+E_c_hat_rho(i)**1.5)*f_alpha(k,j,i))/dE
        enddo

        do j = 2,nE-1
          do k = 1,ny
            term_csd(k,j,i) = 2./tau_s_rho(i)/E_hat(j)**0.5* &
                    ((E_hat(j+1)**1.5+E_c_hat_rho(i)**1.5)*f_alpha(k,j+1,i)- &
                    (E_hat(j)**1.5+E_c_hat_rho(i)**1.5)*f_alpha(k,j,i))/dE
          enddo
        enddo

        if(N_source .eq. 1) then
          j = nE-nEa
          if(i_isotropic) then !default, isotropic source
            do k = 1,ny
              term_s(k,j,i) = S0_rho(i)/2.0/pi/dE
            enddo
          else
            do k = 1,ny
              term_s(k,j,i) = S0_rho(i)/2.0/pi/dE*F_bpa(k)/Int_F_bpa
            enddo
          endif
        elseif(N_source .eq. 3) then
          do n = 1,N_source
            j = j_NBI(n)
            if(i_isotropic) then
              do k = 1,ny
                term_s(k,j,i) = S0_rho(i)/2.0/pi/(n/N_source)**0.5/dE
              enddo
            else !default, beam-like source
              do k = 1,ny
                term_s(k,j,i) = S0_rho(i)/2.0/pi/(n/N_source)**0.5/dE * &
                                F_bpa(k)/Int_F_bpa
              enddo
            endif
          enddo
        endif

        if(i_scattering) then
        !Pitch angle scattering, nu_d*d/dy((1-y**2)df/dy)
        !that is: nu_d*((1-y**2)*ddf/dy**2 - 2*y*df/dy)
          do j = 1,nE
            do k = 1,ny
              if(k .eq. 1) then
                kd = 1
              else
                kd = k-1
              endif
              if(k .eq. ny) then
                ku = ny
              else
                ku = k+1
              endif
              term_pas(k,j,i) = nu_d(j,i)*((1-y(k)**2)*&
                      (f_alpha(ku,j,i)-2*f_alpha(k,j,i)+f_alpha(kd,j,i))/dy**2&
                       - 2.*y(k)*(f_alpha(ku,j,i)-f_alpha(kd,j,i))/2./dy)
            enddo
          enddo
        endif

      enddo
      !$OMP END PARALLEL DO

      !calculate the radial flux term
      if(i_flux_r) then

        !AE on, then update the density gradient and compute the new D_rr
        if(i_AE) then !add pressure gradient threshold, 7.13.2015
          if(threshold_flag .eq. 0) then !AE density gradient threshold
            !$OMP PARALLEL DO PRIVATE(i)
            do i = 1,nr
              n_tran_alpha_rho(i) = Int_E1y(f_alpha(:,:,i),nE)
            enddo
            !$OMP END PARALLEL DO

            rg_n_tran_alpha_rho = rdiff(n_tran_alpha_rho)/dr

            !08.29.16, CDG correction factor
            L_tran_n_alpha_rho = n_tran_alpha_rho/rg_n_tran_alpha_rho
            L_tran_n_alpha_rho(1) = L_tran_n_alpha_rho(2)

            !$OMP PARALLEL DO PRIVATE(i)
            do i = 1,nr
              F_cor(i) = (1.0 + L_s_n_alpha_rho(i)*(1./L_s_T_alpha_rho(i)-C_R/Rmaj_rho(i)))/ &
                         (1.0 + L_tran_n_alpha_rho(i)*(1./L_s_T_alpha_rho(i)-C_R/Rmaj_rho(i)))
            enddo
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO PRIVATE(i,j,k,dn_th)
            do i = 1,nr
              if(F_cor(i) .gt. 0.) then
                dn_th = rg_n_tran_alpha_rho(i) - rg_n_alpha_th_rho(i) * F_cor(i) !08.29.16, CDG correction factor
              else
                dn_th = rg_n_tran_alpha_rho(i) - rg_n_alpha_th_rho(i) * 0.0
              endif
              do j = 1,nE
                do k = 1,ny
                  if(dn_th .gt. 0.) then
                    !!revised 5.29.2015
                    D_rr(k,j,i) = D_AE*rmin/n_tran_alpha_rho(i)*dn_th &
                                  *G_AE(j,i) + D_ITG(k,j,i)
                  else
                    D_rr(k,j,i) = D_ITG(k,j,i)
                  endif
                enddo
              enddo
            enddo
            !$OMP END PARALLEL DO
          else  !AE pressure gradient threshold
            !$OMP PARALLEL DO PRIVATE(i)
            do i = 1,nr
              p_tran_alpha_rho(i) = 1.e3*2./3.*E_alpha* &
                                    Int_E2y(f_alpha(:,:,i),nE)*0.16022
            enddo
            !$OMP END PARALLEL DO

            rg_p_tran_alpha_rho = rdiff(p_tran_alpha_rho)/dr

            !$OMP PARALLEL DO PRIVATE(i,j,k,dp_th)
            do i = 1,nr
              dp_th = rg_p_tran_alpha_rho(i) - rg_p_alpha_th_rho(i)
              do j = 1,nE
                do k = 1,ny
                  if(dp_th .gt. 0.) then
                    D_rr(k,j,i) = D_AE*rmin/p_tran_alpha_rho(i)*dp_th &
                                  *G_AE(j,i) + D_ITG(k,j,i)
                  else
                    D_rr(k,j,i) = D_ITG(k,j,i)
                  endif
                enddo
              enddo
            enddo
            !$OMP END PARALLEL DO
          endif
        endif

        !Update the temperature and A_alpha in D matrix
        if(i_pinch .and. (i_Model .ne. 1)) then
          !$OMP PARALLEL DO PRIVATE(i)
          do i = 1,nr
            n_tran_alpha_rho(i) = Int_E1y(f_alpha(:,:,i),nE)
            if(n_tran_alpha_rho(i) .eq. 0) then
              T_tran_alpha_rho(i) = 0.
            else
              T_tran_alpha_rho(i) = 1.e3*2./3./n_tran_alpha_rho(i)*E_alpha* &
                                    Int_E2y(f_alpha(:,:,i),nE)
            endif
          enddo
          !$OMP END PARALLEL DO

          !C_rE = T_alpha/a*A_alpha*D_rr/1.e3
          !$OMP PARALLEL DO PRIVATE(k,j,i)
          do i = 1,nr
            do j = 1,nE
              do k = 1,ny
                C_rE(k,j,i) = T_tran_alpha_rho(i)/rmin* &
                              A(k,j,i)*D_rr(k,j,i)/1.e3
              enddo
            enddo
          enddo
          !$OMP END PARALLEL DO
        endif

        if(i_Model .eq. 1) then  !for energy independent Angioni model
          !boundary condition, Gamma_r = 0 @r/a = 0.
          !$OMP PARALLEL DO PRIVATE(k,j)
          do j = 1,nE
            do k = 1,ny
              Gamma_r(k,j,1) = 0.
            enddo
          enddo
          !$OMP END PARALLEL DO

          !set Gamma_r
          !Gamma_r = -D_bkg_Angioni*df/dr+D_bkg_Angioni_C_p_alpha/R*f
          !that is, Gamma_r = -D_rr*df/dr+C_rE*f
          !C_rE = 0 for i_pinch = 0, and C_rE =/= 0 for i_pinch = 1
          !$OMP PARALLEL DO PRIVATE(k,j,i)
          do i = 2,nr
            do j = 1,nE
              do k = 1,ny
                Gamma_r(k,j,i) = -D_rr(k,j,i)* &
                                 (f_alpha(k,j,i)-f_alpha(k,j,i-1))/dr &
                                 + C_rE(k,j,i)*f_alpha(k,j,i)
              enddo
            enddo
          enddo
          !$OMP END PARALLEL DO

        else   !i_Model = 0 for const D
               !i_Model = 2 for energy dependent Angioni model
               !i_Model = 3 for DEP model, read from input files
               !i_Model = 4 for Pueschel model
          !BC Gamma_r = 0 @r/a = 0
          !$OMP PARALLEL DO PRIVATE(k,j)
          do j = 1,nE
            do k = 1,ny
              Gamma_r(k,j,1) = 0.
            enddo
          enddo
          !$OMP END PARALLEL DO

          !set Gamma_r
          !Gamma_r = -D_rr*df/dr+C_rE/E_alpha*df/dE_hat
          !$OMP PARALLEL DO PRIVATE(k,j,i)
          do i = 2,nr
            do j = 1,nE-1
              do k = 1,ny
                Gamma_r(k,j,i) = -D_rr(k,j,i)* &
                                 (f_alpha(k,j,i)-f_alpha(k,j,i-1))/dr &
                                 + C_rE(k,j,i)/E_alpha* &
                                 (f_alpha(k,j+1,i)-f_alpha(k,j,i))/dE
              enddo
            enddo
          enddo
          !$OMP END PARALLEL DO

        endif

        !set radial flux term 1/V'*d{V'*Gamma_r}/dr=1/V'*(dV'/dr)*Gamma_r + d(Gamma_r)/dr

        !Because Gamma_r(i=1)=0, the first term is 0. d(Gamma_r)/dr=(Gamma_r(k,j,i+1)-Gamma_r(k,j,i))/dr
        !Gamma_r(i) is actually the flux at i-1/2, so d(Gamma_r)/dr(i=1) = Gamma_r(k,j,2)/(dr/2)
        !BC 1/V'*d{V'*Gamma_r}/dr => 2.*Gamma_r(i=2)/dr @r/a = 0,
        !$OMP PARALLEL DO PRIVATE(j,k)
        do j = 1,nE
          do k = 1,ny
            term_flux_r(k,j,1) = 2.*Gamma_r(k,j,2)/dr
            !term_flux_r(k,j,1) = 4.*Gamma_r(k,j,2)/dr  !6.14.2015
          enddo
        enddo
        !$OMP END PARALLEL DO


        !$OMP PARALLEL DO PRIVATE(k,j,i)
        do i = 2,nr-1
          do j = 1,nE
            do k = 1,ny
              term_flux_r(k,j,i) = 1./V_prime_rho(i)&
                                   *(V_prime_rho(i+1)-V_prime_rho(i-1))/2./dr&
                                   *(Gamma_r(k,j,i)+Gamma_r(k,j,i+1))/2.&
                                   +(Gamma_r(k,j,i+1)-Gamma_r(k,j,i))/dr
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
        !boundary condition
        i = nr
        !$OMP PARALLEL DO PRIVATE(j,k)
        do j = 1,nE
          do k = 1,ny
            term_flux_r(k,j,i) = 1./V_prime_rho(i)&
                                 *(V_prime_rho(i)-V_prime_rho(i-1))/dr&
                                 *Gamma_r(k,j,i)&
                                 +(Gamma_r(k,j,i)-Gamma_r(k,j,i-1))/dr
          enddo
        enddo
        !$OMP END PARALLEL DO

      endif

      !calculate the energy flux term
      if(i_flux_E) then
        !C_EE = C_rE**2/D_rr
        !$OMP PARALLEL DO PRIVATE(k,j,i)
        do i = 1,nr
          do j = 1,nE
            do k = 1,ny
              if(D_rr(k,j,i) .eq. 0) then
                C_EE(k,j,i) = 0.
              else
                C_EE(k,j,i) = C_rE(k,j,i)**2/D_rr(k,j,i)
              endif
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

        !BC Gamma_E at r/a = 0
        i = 1
        !$OMP PARALLEL DO PRIVATE(j,k)
        do j = 1,nE-1
          do k = 1,ny
            Gamma_E(k,j,i) = -C_EE(k,j,i)/E_alpha*&
                             (f_alpha(k,j+1,i)-f_alpha(k,j,i))/dE
          enddo
        enddo
        !$OMP END PARALLEL DO

        !set Gamma_E
        !Gamma_E = -C_rE*df/dr-C_EE/E_alpha*df/dE_hat
        !$OMP PARALLEL DO PRIVATE(k,j,i)
        do i = 2,nr
          do j = 1,nE-1
            do k = 1,ny
              Gamma_E(k,j,i) = -C_rE(k,j,i)*&
                               (f_alpha(k,j,i)-f_alpha(k,j,i-1))/dr-&
                               C_EE(k,j,i)/E_alpha*&
                               (f_alpha(k,j+1,i)-f_alpha(k,j,i))/dE
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

        !set energy flux term 1/V_E*d{V_E*f}/dE, dE = dE-lambda/E*d(lambda)
        !the first term : term_flux_E = d(Gamma_E)/dE - Gamma_E/2/E
        !for V_E = 1/sqrt(E(1-lambda)), 2015.7.2
        !$OMP PARALLEL DO PRIVATE(k,j,i)
        do i = 1,nr
          !boundary condition !!!confusion here
          do k = 1,ny
            term_flux_E(k,1,i) = 0. !6.14.2015
          enddo

          do j = 2,nE-1
            do k = 1,ny
              term_flux_E(k,j,i) = (Gamma_E(k,j,i)-Gamma_E(k,j-1,i))/E_alpha/dE&
                                   - Gamma_E(k,j,i)/2./E_hat(j)/E_alpha
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

        !the second term : term_flux_y = (1-y**2)/2/E/y*(d(Gamma_E)/dy-Gamma_E/y),
        !that is : term_flux_y = ((1-y**2)/2/E*(d[(Gamma_E)/y]/dy),
        !confirmed 7.2.2015, 2.25.2016
        !$OMP PARALLEL DO PRIVATE(i,j,k,kd,ku)
        do i = 1,nr
          do k = 1,ny
            term_flux_y(k,1,i) = 0. !6.14.2015
          enddo

          do j = 2,nE-1
            do k = 1,ny
              if(k .eq. 1) then
                kd = 1
              else
                kd = k-1
              endif
              if(k .eq. ny) then
                ku = ny
              else
                ku = k+1
              endif
              !term_flux_y(k,j,i) = (1-y(k)**2)/2./E_hat(j)/E_alpha*&
              !                  (Gamma_E(ku,j,i)/y(ku)-Gamma_E(kd,j,i)/y(kd))/2./dy
              !the correct method, 7.2.2015, but a 1/2 error
              term_flux_y(k,j,i) = (1-y(k)**2)/2./E_hat(j)/E_alpha*&
                                (Gamma_E(k,j,i)/y(k)-Gamma_E(kd,j,i)/y(kd))/2./dy
              !correct 2.25.2016
              !term_flux_y(k,j,i) = (1-y(k)**2)/2./E_hat(j)/E_alpha*&
              !                  (Gamma_E(k,j,i)/y(k)-Gamma_E(kd,j,i)/y(kd))/dy
              !term_flux_y(k,j,i) = (1-y(k)**2)/2./E_hat(j)/E_alpha*&
              !                  (Gamma_E(ku,j,i)/y(ku)-Gamma_E(k,j,i)/y(k))/dy
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

      endif


      !advance f_alpha(r,E)
      !$OMP PARALLEL DO PRIVATE(k,j,i)
      do i = 1,nr
        do j = 1,nE-1
          do k = 1,ny
            f_alpha(k,j,i) = f_alpha_bk(k,j,i) + tstep * &
                             (term_csd(k,j,i) + term_s(k,j,i) &
                             + term_pas(k,j,i) - term_flux_r(k,j,i) &
                             - term_flux_E(k,j,i) - term_flux_y(k,j,i))
          enddo
        enddo

        !zero fixed BC @Energy upper boundary
        j = nE
        do k = 1,ny
          f_alpha(k,j,i) = 0.
        enddo
        !free BC @E=0, 6.25.2015
        !j = 1
        !do k = 1,ny
        !  f_alpha(k,j,i) = f_alpha(k,j+1,i)
        !enddo
      enddo
      !$OMP END PARALLEL DO

      !BC for f_alpha = delta1*F_s @r/a = 1
      if(i_BC_r1) then
        i = nr
        !$OMP PARALLEL DO PRIVATE(j,k)
        do j = 1,nE-1
          do k = 1,ny
            f_alpha(k,j,i) = delta1*F_s(j,i)
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif

      !if(irBC_outer .ne. nr) then
      !!add new outside BC for orbit width/ drift displacment
      !  !$OMP PARALLEL DO PRIVATE(j,k,i)
      !  do i = irBC_outer,nr
      !    do j = 1,nE-1
      !      do k = 1,ny
      !        f_alpha(k,j,i) = 0.
      !      enddo
      !    enddo
      !  enddo
      !  !$OMP END PARALLEL DO
      !endif

      if((mod(nt-nstart,ndt) .eq. 0) .and. (nrk .eq. 2) &
         .and. (ndata .le. nsav)) then !save data
        ndata = ndata + 1
        print *, 'ndata = ',ndata

        !!check  numerically instability, 3.20.2015
        do i = 1,nr
          do j = 1,nE
            do k = 1,ny
              if(f_alpha(k,j,i) .lt. 0.) then
                print *, 'error, minus f_alpha numerical instabilty'
                stop
              elseif(f_alpha(k,j,i) .ge. 0.) then
                !!do nothing, correct case
              else
                print *, 'error, NaN f_alpha numerical instability'
                stop
              endif
            enddo
          enddo
        enddo

        open(unit=11,file='f_bac.dat',status='replace',form='unformatted')
        write(11) real(nt)
        write(11) f_alpha
        close(11)
      endif

    enddo
  enddo
  !the main loop end

  f_error_max = maxval(abs(f_alpha-f_alpha_bk))

  !--------------------------------------------------------------
  !Outputs
  !--------------------------------------------------------------

  !The information of profiles
  open(unit=4,file='EPtran_transport.out',status='replace')

  write(4,*) 'max f_alpha error is',f_error_max
  write(4,*) '--------------------------------------------------------------'
  if(i_flux_r) then
    write(4,*) 'With radial flux'

    if(i_flux_E) then
      write(4,*) 'With energy flux'
    endif

    select case(i_Model)
    case(1)
      write(4,*) 'Using energy independent Angioni model'
    case(2)
      write(4,*) 'Using energy dependent Angioni model'
    case(3)
      write(4,*) 'Read input file, like DEP model'
    case(4)
      write(4,*) 'Using Pueschel model'
      if(Pueschel_em_flag) then
        write(4,*) 'with electromagnetic counterparts'
        write(4,*) '(beta/beta_crit)**2 =',beta_square
      else
        write(4,*) 'Only electrostatic part'
      endif
    case default
      write(4,*) 'Constant D matrix'
      if(i_energy_AE) then
        write(4,*) 'G_E =/= 1'
        write(4,*) 'eps_hat_0 =',eps_hat_0
        write(4,*) 'Delta_eps_hat =',Delta_eps_hat
      else
        write(4,*) 'G_E == 1'
      endif
    end select

    if(i_AE) then
      if(threshold_flag .eq. 0) then
        write(4,*) 'The density gradient threshold'
      else
        write(4,*) 'The pressure gradient threshold'
      endif

      if(i_energy_AE) then
        write(4,*) 'G_AE =/= 1'
        write(4,*) 'eps_hat_0 =',eps_hat_0
        write(4,*) 'Delta_eps_hat =',Delta_eps_hat
      else
        write(4,*) 'G_AE == 1'
      endif
      write(4,*) 'AE on, D_AE=',D_AE
    else
      write(4,*) 'AE off'
    endif

    if(i_pinch) then
      write(4,*) 'Energy pinch on in D matrix'
    else
      write(4,*) 'Energy pinch off in D matrix'
    endif

    if(i_thermal_pinch_off) then
      write(4,*) 'Turn off the thermal pinch'
    endif

    if(D_factor .ne. 1) then
      write(4,*) 'Increase D matrix by ',D_factor
    endif
    if(A_factor .ne. 1) then
      write(4,*) 'Increase A by ',A_factor
    endif
  else
    write(4,*) 'No transport'
  endif

  if(i_BC_r1) then
    write(4,*) 'Fixed BC f_alpha = delta1*F_s at r/a =1, delta1=',delta1
  else 
    write(4,*) 'Free BC at r/a = 1'
  endif

  write(4,*) '**************************************************************'
  write(4,*) 'The parameters'
  write(4,*) '**************************************************************'
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'dt =',dt
  write(4,*) 'dr =',dr
  write(4,*) 'dE =',dE
  write(4,*) 'dy =',dy
  write(4,*) 'ntstep =',ntstep
  write(4,*) 'nr =',nr
  write(4,*) 'nE =',nE
  write(4,*) 'ny =',ny
  write(4,*) 'NBI_flag =',NBI_flag
  if(N_source .eq. 3) then
    write(4,*) 'the full, half, and third components energy source in DIII-D'
    write(4,10) j_NBI(1),j_NBI(2),j_NBI(3)
  endif
  if(i_isotropic) then
    write(4,*) 'isotropic source'
  else
    write(4,*) 'beam-like source'
    write(4,*) 'lambda0 = ',lambda0
    write(4,*) 'Delta_lambda = ',Delta_lambda
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'F_bpa(+1,lambda) and lambda grids'
    write(4,*) '--------------------------------------------------------------'
    do k = 1,ny
      write(4,*) F_bpa(k)/Int_F_bpa,'at lambda =',1.-y(k)**2
    enddo
  endif
  if(i_scattering) then
    write(4,*) 'pitch angle scattering'
  else
    write(4,*) 'no pitch angle scattering'
  endif
  write(4,*) '--------------------------------------------------------------'

  !1.2.2016, new Fortran output menthod with array
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'tau_s(r) profile in second'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) tau_s_rho

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'S0(r) profile in 10**19/m**3/s'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) S0_rho

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'E_c/E_alpha profile'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) E_c_hat_rho

  !6.21.2015
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The thermal pinch reverse energy E_r/E_alpha profile'
  write(4,*) 'E_r = E_c*[log(1/E_c_hat**1.5+1)*(1+E_c_hat**1.5)-1]**(2/3)'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) E_c_hat_rho(:)*&
          (log(1/E_c_hat_rho(:)**1.5+1.)*(1.+E_c_hat_rho(:)**1.5)-1)**(2./3.)

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'V_prime(r) profile in m**2'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) V_prime_rho

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'chi_eff(r) profile'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) chi_eff

  if(i_Model .eq. 1)then
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'D_bkg_Angioni(r) profile'
    write(4,*) '--------------------------------------------------------------'
    write(4,10) D_bkg_Angioni

    if(i_pinch)then
      write(4,*) '--------------------------------------------------------------'
      write(4,*) 'C_p_alpha(r) profile'
      write(4,*) '--------------------------------------------------------------'
      write(4,10) C_p_alpha
    endif
  endif

  if(i_Model .eq. 2)then
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'c = c_G_D_s/G_D_s profile'
    write(4,*) '--------------------------------------------------------------'
    write(4,10) c

    if(i_pinch)then !add 'if' 1/2/2016
      write(4,*) '--------------------------------------------------------------'
      write(4,*) 'k_p profile'
      write(4,*) '--------------------------------------------------------------'
      write(4,10) k_p
    endif
  endif

  write(4,*) '**************************************************************'
  write(4,*) 'The slowing down profiles (no transport)'
  write(4,*) '**************************************************************'
  !The static slowing down profiles
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The slowing down density profile ns(r) in 10**19/m**3 analysis'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) n_alpha_rho

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The slowing down density gradient profile in 10**19/m**4'
  write(4,*) '--------------------------------------------------------------'
  !dn(r)/dr = (n(i+1)-n(i-1)/2/dr
  write(4,10) rdiff(n_alpha_rho)/dr

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The slowing down equivalent temperature profile Ts(r) in keV analysis'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) T_alpha_equiv_rho

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The slowing down temperature gradient profile in keV/m'
  write(4,*) '--------------------------------------------------------------'
  !dT(r)/dr = (T(i+1)-T(i-1))/2/dr
  write(4,10) rdiff(T_alpha_equiv_rho)/dr

  !5.13.15, revised 5.28.15, 7.13.2015
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The slowing down pressure profile ps(r) in 10**4 N/m**2'
  write(4,*) '--------------------------------------------------------------'
  p_alpha_rho(:) = n_alpha_rho(:)*T_alpha_equiv_rho(:)*0.16022
  write(4,10) p_alpha_rho

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The slowing down pressure gradient profile in [10**4 N/m**2]/m'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) rdiff(p_alpha_rho)/dr

  !the code integral static alpha profiles
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'n0(r) profile in 10**19/m**3 integral from F_s(r,E)'
  write(4,*) '--------------------------------------------------------------'
  do i = 1,nr
    write(4,10) Int_E(F_s(:,i),nE)
  enddo

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'T0(r) profile in keV integral from F_s(r,E)'
  write(4,*) '--------------------------------------------------------------'
  do i = 1,nr
    write(4,10) 1.e3*2./3./Int_E(F_s(:,i),nE)* &
                E_alpha*Int_E(E_hat(:)*F_s(:,i),nE)
  enddo

  write(4,*) '**************************************************************'
  write(4,*) 'The final stable profiles (transported)'
  write(4,*) '**************************************************************'
  !the n_tran_alpha_rho(i) and T_tran_alpha_rho(i) profiles
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The final n_alpha(r) profile in 10**19/m**3'
  write(4,*) '--------------------------------------------------------------'
  !n(r) = Int_E(f_alpha(r,E))
  do i = 1,nr
    n_tran_alpha_rho(i) = Int_E1y(f_alpha(:,:,i),nE)
  enddo
  write(4,10) n_tran_alpha_rho

  !the density gradient profiles
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The density gradient profile in 10**19/m**4'
  write(4,*) '--------------------------------------------------------------'
  !-dn(r)/dr = -(n(i+1)-n(i-1)/2/dr
  rg_n_tran_alpha_rho = rdiff(n_tran_alpha_rho)/dr
  write(4,10) rg_n_tran_alpha_rho

  if(i_AE .and. threshold_flag .eq. 0) then
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The density gradient threshold in 10**19/m**4'
    write(4,*) '--------------------------------------------------------------'
    write(4,10) rg_n_alpha_th_rho

    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The density gradient threshold correction factor'
    write(4,*) '--------------------------------------------------------------'
    L_tran_n_alpha_rho = n_tran_alpha_rho/rg_n_tran_alpha_rho
    L_tran_n_alpha_rho(1) = L_tran_n_alpha_rho(2)
    F_cor = (1.0 + L_s_n_alpha_rho*(1./L_s_T_alpha_rho-C_R/Rmaj_rho))/ &
            (1.0 + L_tran_n_alpha_rho*(1./L_s_T_alpha_rho-C_R/Rmaj_rho))
    write(4,10) F_cor

    if(i_orbit) then
      write(4,*) 'The inner BC at r/a =',real(irBC_inner-1)/real(nr-1)
      write(4,*) 'The outer BC at r/a =',real(irBC_outer-1)/real(nr-1)

      write(4,*) '--------------------------------------------------------------'
      write(4,*) 'The average normalized Orbit Width Delta_orb_ave/a'
      write(4,*) '--------------------------------------------------------------'
      write(4,10) Delta_orb_ave_hat

      write(4,*) '--------------------------------------------------------------'
      write(4,*) 'The original critical gradient profile'
      write(4,*) '--------------------------------------------------------------'
      write(4,10) rg_n_alpha_th_rho_bk
    endif

  endif

  !the temperature profiles
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The final T_alpha(r) profile in keV'
  write(4,*) '--------------------------------------------------------------'
  !T(r) = 1e3*2/3/n(r)*E_alpha*Int_E(E*f_alpha(r,E))
  do i = 1,nr
    T_tran_alpha_rho(i) = 1.e3*2./3./n_tran_alpha_rho(i)*E_alpha* &
                          Int_E2y(f_alpha(:,:,i),nE)
  enddo
  write(4,10) T_tran_alpha_rho

  !the temperature gradient profiles
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The temperature gradient profile in keV/m'
  write(4,*) '--------------------------------------------------------------'
  !-dT(r)/dr = -(T(i+1)-T(i-1))/2/dr
  rg_T_tran_alpha_rho = rdiff(T_tran_alpha_rho)/dr
  write(4,10) rg_T_tran_alpha_rho

  !the pressure profiles, 5.13.15
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The final pressure profile in 10**4 N/m**2'
  write(4,*) '--------------------------------------------------------------'
  do i = 1,nr
    p_tran_alpha_rho(i) = 1.e3*2./3.*E_alpha*Int_E2y(f_alpha(:,:,i),nE)*0.16022
  enddo
  write(4,10) p_tran_alpha_rho

  !the pressure gradient profiles, 7.13.15
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The pressure gradient profile in [10**4 N/m**2]/m'
  write(4,*) '--------------------------------------------------------------'
  rg_p_tran_alpha_rho = rdiff(p_tran_alpha_rho)/dr
  write(4,10) rg_p_tran_alpha_rho

  if(i_AE .and. threshold_flag .eq. 1) then
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The pressure gradient threshold in [10**4 N/m**2]/m'
    write(4,*) '--------------------------------------------------------------'
    write(4,10) rg_p_alpha_th_rho
  endif

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The pressure gradient check profile T*dn/dr + n*dT/dr'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) (T_tran_alpha_rho(:)*rg_n_tran_alpha_rho(:) + &
               n_tran_alpha_rho(:)*rg_T_tran_alpha_rho(:))*0.16022

  !particle and energy flows
  if(i_flux_r) then
    !the particle flow(i) profiles
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The particle flow profile in 10**19/sec'
    write(4,*) '--------------------------------------------------------------'
    !flow(r) = S(r)*flux(r) = V_prime(r)*flux(r)
    do i = 1,nr
      flow_rho(i) = V_prime_rho(i)*Int_E1y(Gamma_r(:,:,i),nE)
    enddo
    write(4,10) flow_rho

    !the particle flow_check(i) profiles
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The particle flow check profile in 10**19/sec'
    write(4,*) '--------------------------------------------------------------'
    !flow_check(r) = Int(0 to r/a) x dV(r) x [source term - "slowing down term" sink term]
    !dV(r) = V_prime(r)*dr
    flow_check_rho(1) = 0.
    do i = 2,nr
      flow_check_rho(i) = flow_check_rho(i-1)+V_prime_rho(i)* &
                          Int_E1y(term_flux_r(:,:,i),nE)*dr
    enddo
    write(4,10) flow_check_rho

    !final D_alpha, 5.13.2015
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The final D_alpha = flux(r)/(-dn/dr) in m**2/sec'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) Int_E1y(Gamma_r(:,:,i),nE)/rg_n_tran_alpha_rho(i)
    enddo

    !D_star, 2.27.2016
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'D_rr(r), E and lambda have been integrated with F_s'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      if(i_isotropic) then
        do j = 1,nE
          check_E(j) = Int_y(D_rr(:,j,i))
        enddo
      else
        do j = 1,nE
          check_E(j) = Int_y(D_rr(:,j,i)*F_bpa(:))/Int_F_bpa
        enddo
      endif
      write(4,10) Int_E(check_E(:)*F_s(:,i),nE)/Int_E(F_s(:,i),nE)
    enddo

    !the energy flow(i) profiles
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The energy flow profile in MW'
    write(4,*) '--------------------------------------------------------------'
    !flow_energy(r) = V_prime(r)*Int_E(Gamma_r(:,i)*E_hat)*1.6022*E_alpha
    do i = 1,nr
      flow_energy_rho(i) = V_prime_rho(i)* &
                           Int_E2y(Gamma_r(:,:,i),nE)*1.6022*E_alpha
    enddo
    write(4,10) flow_energy_rho

    !final D_alpha_p, 5.14.2015
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The final D_alpha_p = flux_E(r)/(-3/2*dp/dr) in m**2/sec'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      !revised, 7.13.2015
      write(4,10) Int_E2y(Gamma_r(:,:,i),nE)*E_alpha*1000. &
                  /(3./2.*rg_p_tran_alpha_rho(i)/0.16022)
    enddo

    !Convec_factor, 5.14.2015, correct T_EP, 6.21.2015
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The Convec_factor = flux_E(r)/[3/2*T(r)*flux(r)]'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) (Int_E2y(Gamma_r(:,:,i),nE)*E_alpha*1000.)/&
                  (3./2.*T_tran_alpha_rho(i)*Int_E1y(Gamma_r(:,:,i),nE))
    enddo

    !6.21.2015
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The check Convec_factor = 1 + n(r)*dT(r)/dr/(T(r)*dn(r)/dr)'
    write(4,*) '--------------------------------------------------------------'
    write(4,10) 1. + n_tran_alpha_rho(:)*rg_T_tran_alpha_rho(:)/&
                    (T_tran_alpha_rho(:)*rg_n_tran_alpha_rho(:))

    !The "cold" and "hot" profiles
    write(4,*) '--------The cold (E<1/2E0) and hot (E>1/2E0) profiles---------'
    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The no transport cold density profile in 10**19/m**3'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) Int_E(F_s(:,i),nint((nE-nEa-1)/2.))
    enddo

    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The no transport hot density profile in 10**19/m**3'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) Int_E(F_s(:,i),nE)-Int_E(F_s(:,i),nint((nE-nEa-1)/2.))
    enddo

    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The transported cold density profile in 10**19/m**3'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) Int_E1y(f_alpha(:,:,i),nint((nE-nEa-1)/2.))
    enddo

    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The transported hot density profile in 10**19/m**3'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) Int_E1y(f_alpha(:,:,i),nE)-Int_E1y(f_alpha(:,:,i),nint((nE-nEa-1)/2.))
    enddo

    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The cold particle flow profile in 10**19/sec'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) V_prime_rho(i)*Int_E1y(Gamma_r(:,:,i),nint((nE-nEa-1)/2.))
    enddo

    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The hot particle flow profile in 10**19/sec'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) V_prime_rho(i)*(Int_E1y(Gamma_r(:,:,i),nE) - &
                 Int_E1y(Gamma_r(:,:,i),nint((nE-nEa-1)/2.)))
    enddo

    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The cold energy flow profile in MW'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) V_prime_rho(i)* &
                  Int_E2y(Gamma_r(:,:,i),nint((nE-nEa-1)/2.))*1.6022*E_alpha
    enddo

    write(4,*) '--------------------------------------------------------------'
    write(4,*) 'The hot energy flow profile in MW'
    write(4,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(4,10) V_prime_rho(i)*(Int_E2y(Gamma_r(:,:,i),nE) - &
                 Int_E2y(Gamma_r(:,:,i),nint((nE-nEa-1)/2.)))*1.6022*E_alpha
    enddo

  endif

  !The source alpha particles
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The alpha particle source in 10**19/sec'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) source_rho

  !source flow, 5.13.2015
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The source flow in 10**19/sec'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) flow_source_rho

  !The total source particles, sum(S0_rho(i)*V_prime(i)*dr)
  InPar = sum(source_rho)
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The total source alpha particles in 10**19/sec,',InPar
  write(4,*) '--------------------------------------------------------------'

  !The source fusion power
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The alpha source fusion power in MW'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) source_rho*E_alpha*1.6022

  !source energy flow, 5.13.2015
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The source energy flow in MW'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) flow_source_rho*E_alpha*1.6022

  !The total fusion power, sum(S0_rho(i)*V_prime(i)*E_alpha)
  InPower = InPar*E_alpha*1.6022
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The total fusion power in MW,',InPower
  write(4,*) '--------------------------------------------------------------'

!!!11.12.14
  !The alpha particles sink and the He source
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The alpha particle sink and the He source in 10**19/sec'
  write(4,*) '--------------------------------------------------------------'
  i = 1
  sink_rho(i) = 0.
  do i = 2,nr-1
    sink_rho(i) = source_rho(i)+flow_rho(i)-flow_rho(i+1)
  enddo
  i = nr
  sink_rho(i) = source_rho(i)+flow_rho(i)
  write(4,10) sink_rho

!!!11.12.14
  !The energy sink for slowing down
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The energy sink for slowing down in MW,'
  write(4,*) '--------------------------------------------------------------'
  i = 1
  sinke_rho(i) = 0.
  do i = 2,nr-1
    sinke_rho(i) = source_rho(i)*E_alpha*1.6022 + &
                   flow_energy_rho(i)-flow_energy_rho(i+1)
  enddo
  i = nr
  sinke_rho(i) = source_rho(i)*E_alpha*1.6022+flow_energy_rho(i)
  write(4,10) sinke_rho

  !the He ash source profile
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The He ash source profile in 10**19/m**3/sec'
  write(4,*) '--------------------------------------------------------------'
  write(4,10) sink_rho/V_prime_rho/dr

  !the He ash particle flow profiles
  flux_He_rho(1)=0.
  do i = 2,nr
    flux_He_rho(i) = V_prime_rho(i-1)/V_prime_rho(i)*flux_He_rho(i-1) + &
                    0.5/V_prime_rho(i)*(sink_rho(i)+sink_rho(i-1))
  enddo

  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The He ash particle flow profile in 10**19/sec'
  write(4,*) '--------------------------------------------------------------'
  flow_He_rho(:) = V_prime_rho(:)*flux_He_rho(:)
  write(4,10) flow_He_rho

  !the He ash density profiles
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The He ash density profile (no alpha transport) in 10**19/m**3'
  write(4,*) '--------------------------------------------------------------'
  n_s_He_rho(nr) = 0.
  n_s_He_rho(nr-1) = ((flux_source_rho(nr)/D_bkg_Angioni_He(nr) &
                    +flux_source_rho(nr-1)/D_bkg_Angioni_He(nr-1))/2. &
                    *(rho_hat(nr)-rho_hat(nr-1))*rmin)/ &
           (1.+C_p_He(nr-1)*(rho_hat(nr)-rho_hat(nr-1))*rmin/Rmaj_rho(nr-1))
  do i = nr-1,2,-1
    n_s_He_rho(i-1) = n_s_He_rho(i+1)  &
                     -n_s_He_rho(i)*C_p_He(i)/Rmaj_rho(i)&
                     *(rho_hat(i+1)-rho_hat(i-1))*rmin &
                     +flux_source_rho(i)/D_bkg_Angioni_He(i)&
                     *(rho_hat(i+1)-rho_hat(i-1))*rmin
  enddo
  write(4,10) n_s_He_rho

  !the He ash density profiles
  write(4,*) '--------------------------------------------------------------'
  write(4,*) 'The He ash transported density profile in 10**19/m**3'
  write(4,*) '--------------------------------------------------------------'
  n_tran_He_rho(nr) = 0.
  n_tran_He_rho(nr-1) = ((flux_He_rho(nr)/D_bkg_Angioni_He(nr) &
                       +flux_He_rho(nr-1)/D_bkg_Angioni_He(nr-1))/2. &
                       *(rho_hat(nr)-rho_hat(nr-1))*rmin)/ &
           (1.+C_p_He(nr-1)*(rho_hat(nr)-rho_hat(nr-1))*rmin/Rmaj_rho(nr-1))
  do i = nr-1,2,-1
    n_tran_He_rho(i-1) = n_tran_He_rho(i+1)  &
                        -n_tran_He_rho(i)*C_p_He(i)/Rmaj_rho(i)&
                        *(rho_hat(i+1)-rho_hat(i-1))*rmin &
                        +flux_He_rho(i)/D_bkg_Angioni_He(i)&
                        *(rho_hat(i+1)-rho_hat(i-1))*rmin
  enddo
  write(4,10) n_tran_He_rho

  close(4)

  !The information of distribution
  open(unit=5,file='EPtran_distribution.out',status='replace')

  !write out the whole distribution
  write(5,*) '--------------------------------------------------------------'
  if(i_start_null) then
    write(5,*) 'The start distribution is f_alpha = 0.'
  else
    write(5,*) 'The start distribution is f_alpha = F_s'
  endif
  write(5,*) '--------------------------------------------------------------'

  i = nr     !The edge flow(r,E) vesus E, 5.14.2015
  write(5,*) '--------------------------------------------------------------'
  write(5,*) 'The edge flow(r,E) vesus E at r/a = ',rho_hat(i)
  write(5,*) '--------------------------------------------------------------'
  do j = 1,nE
    write(5,10) V_prime_rho(i)*Int_y(Gamma_r(:,j,i))*2*pi*sqrt(E_hat(j))
  enddo

  write(5,*) '--------------------------------------------------------------'
  write(5,*) 'sqrt(E)*f_s(r,E)/n_alpha_rho(r) vs E at different radius'
  write(5,*) '--------------------------------------------------------------'
  do i_out = 1,5  ! r = 0.1, 0.3, 0.5, 0.7, 0.9
    i = (nr-1)/5*(i_out-1)+(nr-1)/10+1
    write(5,*) 'r/a = ',rho_hat(i)
    write(5,*) '--------------------------------------------------------------'
    do j = 1,nE
      write(5,10) sqrt(E_hat(j))*F_s(j,i)/n_alpha_rho(i)
    enddo
  enddo

  write(5,*) '--------------------------------------------------------------'
  write(5,*) 'normalized F_alpha(r,E): f/n/(fs/ns) vs E at some radius'
  write(5,*) '--------------------------------------------------------------'
  do i_out = 1,5  ! r = 0.1, 0.3, 0.5, 0.7, 0.9
    i = (nr-1)/5*(i_out-1)+(nr-1)/10+1
    write(5,*) 'r/a = ',rho_hat(i)
    write(5,*) '--------------------------------------------------------------'
    do j = 1,nE
      write(5,10) Int_y(f_alpha(:,j,i))/n_tran_alpha_rho(i) &
                  /(F_s(j,i)/n_alpha_rho(i))
    enddo
  enddo

  if(ny .gt. 1) then
    write(5,*) '--------------------------------------------------------------'
    write(5,*) 'normalized F_alpha(r,E,lamda): f/n/(fs/ns) vs lambda at r/a = 0.5 and some E'
    write(5,*) '--------------------------------------------------------------'
    i = (nr-1)/2+1 ! r/a = 0.5
    do j_out = 1,6 ! E = 0.0,0.2,0.4,0.6,0.8,1.0
      j = (nE-nEa-1)/5*(j_out-1)+1
      write(5,*) 'r/a = ',rho_hat(i)
      write(5,*) 'E/E_alpha = ',E_hat(j)
      write(5,*) '--------------------------------------------------------------'
      do k = 1,ny
        write(5,10) f_alpha(k,j,i)/n_tran_alpha_rho(i)/(F_s(j,i)/n_alpha_rho(i))
      enddo
    enddo
  endif
  close(5)

  if(i_distribution_out) then

    open(unit=10,file='EPtran_f_s.dat',status='replace',form='Unformatted')
    write(10) F_s
    close(10)

    open(unit=10,file='EPtran_f_alpha.dat',status='replace',form='Unformatted')
    write(10) f_alpha
    close(10)

    open(unit=10,file='EPtran_Gamma_r.dat',status='replace',form='Unformatted')
    write(10) Gamma_r
    close(10)

    open(unit=10,file='EPtran_D_rr.dat',status='replace',form='Unformatted')
    write(10) D_rr
    close(10)

  endif

10  format(ES14.7)

end subroutine EPtran_transport
