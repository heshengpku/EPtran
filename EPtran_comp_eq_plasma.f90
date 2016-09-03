!---------------------------------------------------------
! EPtran_comp_eq_plasma.f90
!
! PURPOSE:
!  compute equalibrium and plasma profile
!
!---------------------------------------------------------

subroutine EPtran_comp_eq_plasma

    use EPtran_use_parameter
!INPUTS	                 !nr,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
!                        !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
!                        !Zeff,M_DT,E_alpha,Fpeak_T,Fpeak_n,Fpeak_ei
!OUTPUTS                 !beta_ped_percent,T_ped,beta_N_glob,arho,
!                        !rho_hat(i),q_rho(i),Rmaj_rho(i),rmin_rho(i),
!                        !kappa_rho(i),drho(i),T_i_rho(i),T_e_rho(i),
!                        !n_i_rho(i),n_e_rho(i),n_alpha_rho(i),T_alpha_equiv_rho(i),
!                        !E_c_hat_rho(i)

  !--------------------------------------
  implicit none
  integer :: i
  integer :: i_ITER_profile
  real :: avediv

!
  logical :: iexist
  real :: dummy
  character(len=80) :: comment
  !--------------------------------------


!  i_ITER_profile 
  i_ITER_profile = 0 ! default
!  i_ITER_profile = 1 ! HARDWIRED special Yang Chen et al ITER profiles
!  i_ITER_profile = 2 ! HARDWIRED  DIIID NBI effective E_alpha  test case

! rho_hat grid
  do i = 1,nr
    rho_hat(i) = real(i-1)/real(nr-1)
  enddo


  if(i_ITER_profile .eq. 0) then
  ! various predestal quantites
  ! the actual pedestal width is ignored and the top of the pedesal is at rho_hat = 1.0

    arho = sqrt(kappa_1)*rmin !in meters
    
    beta_ped_percent = beta_N_ped*I_p/B_t/rmin

    T_ped = 0.12*beta_ped_percent*(B_t)**2/n14_ped

    
  !plasma profiles
   
  ! ion temperature in keV
    do i = 1,nr
      T_i_rho(i)=T_ped*(Fpeak_T-1.)*(1.0-rho_hat(i)**2)**2+T_ped
    enddo

  ! electron temperature in keV
    do i = 1,nr
      T_e_rho(i)=T_ped*(Fpeak_T*Fpeak_ei-1.)*(1.0-rho_hat(i)**2)**2+T_ped
    enddo

  ! ion density in 10**19 1/m**3   n19 = 10.*n14     19= 10**19 1/m**3   14=10**14 1/cm**3
    do i = 1,nr
      n_i_rho(i)=10.*n14_ped*(Fpeak_n-1.)*(1.0-rho_hat(i)**2)+10.*n14_ped
    enddo
   
  ! electron desnisty in 10**19 1/m**3
    do i = 1,nr
      n_e_rho(i)=10.*n14_ped*(Fpeak_n-1.)*(1.0-rho_hat(i)**2)+10.*n14_ped
    enddo

  !equillibrium profiles

  !midplane minor radius in meters
    do i = 1,nr
      rmin_rho(i)=rmin*rho_hat(i)
    enddo

    
  !safety factor
    if(eps_q .ge. 0.0) then  !trying to avoid weak central shear
      do i = 1,nr
        if(rho_hat(i) .ge. rho_q) then
          q_rho(i)=q_0*(1+eps_q*rho_hat(i)) &
                   +(rho_hat(i)-rho_q)**2/(1.0-rho_q)**2*(q_1-(1+eps_q)*q_0) &
                   + q_0*eps_q*(rho_hat(i)-rho_q)
        else
          q_rho(i)=q_0*(1.+eps_q*rho_hat(i))
        endif
  ! This coding does not match q' at bounary
  !   if(rho_hat(i).ge.rho_q) then
  !    q_rho(i)=q_rho(i)*((1.+eps_q*rho_q)-(rho_hat(i)-rho_q)/(1.0-rho_q)*eps_q*rho_q)
  !   endif

  !  keeps s_hat > eps_q/(1+eps_q*rho_hat) for rho_hat<rho_q
      enddo
    else!include a q_min at rho_q with q_min = q_0*(1.+eps_q) !if eps_q = 0  q_rho flat inside rho_q
      do i = 1,nr
        if(rho_hat(i) .ge. rho_q) then
          q_rho(i)=q_0*(1.+eps_q) + (q_1-(q_0*(1.+eps_q)))*(rho_hat(i)-rho_q)**2/(1.0-rho_q)**2 
        else
          q_rho(i) = q_0*(1.0+eps_q) -q_0*eps_q*(rho_hat(i)-rho_q)**2/rho_q**2
        endif
      enddo
    endif

  !Shafranov shifted major radius in meters
    do i = 1,nr
      Rmaj_rho(i)=Rmaj0*(1.0+delR0oa*(rmin/Rmaj0)*(1.-rho_hat(i)))
    enddo

  !elongation
    do i = 1,nr
      kappa_rho(i)=(kappa_1-kappa_0)*rho_hat(i)+kappa_0
    enddo

  !triangularity
    do i = 1,nr
      delta_rho(i)=(delta_1-delta_0)*rho_hat(i)+delta_0
    enddo

  ! compute beta_N_glob !should be about 2-3 x beta_N_ped with ITER beta_N_glob limited to ~2.5
    beta_N_glob = 0.
    avediv=0.
    do i = 1,nr
      beta_N_glob = beta_N_glob +kappa_rho(i)*rmin_rho(i)*Rmaj_rho(i)*n_i_rho(i)*(T_i_rho(i)+T_e_rho(i))
      avediv=avediv+kappa_rho(i)*rmin_rho(i)*Rmaj_rho(i)
    enddo
      beta_N_glob = beta_N_glob/avediv  !<n*T_i+n*T_e>
      beta_N_glob = beta_N_glob/(20.*n14_ped*T_ped)*beta_N_ped

  endif !end i_ITER_profile = 0

  if(i_ITER_profile .eq. 1) then
  ! Yang Chen et al HARDWIRED
    B_t = 5.   !T   not 5.3
    I_p = 14.5 !MA  not 15 
    rmin = 2.0 !
    Rmaj0= 6.  !m   not 6.2

    delR0oa=0.16 

    q_1=3.9     !    not 3.8
    q_0=1.3     !    not 1.0
   
    kappa_1=1.8 !    not 1.75
    kappa_0=1.5 

    delta_1=0.5
    delta_0=0.
    
    n14_ped=1.0 ! not 0.9
    !density is dead flat
    
    T_ped = 0.1 !keV   L-mode like not H-mode like Kinsey et al 5.1keV at 30MW Q=10 TGLF predcition
    ! T(r_hat) = (20.-0.1)*(1.-r_hat**2)+0.1  with Te=Ti


    arho = sqrt(kappa_1)*rmin !in meters

    !beta_ped_percent=beta_N_ped*I_p/B_t/rmin

    !T_ped = 0.12*beta_ped_percent*(B_t)**2/n14_ped
   
    beta_ped_percent = T_ped*n14_ped/(0.12*B_t**2)
    beta_N_ped=beta_ped_percent/(I_p/B_t/rmin)


  !plasma profiles

  ! ion temperature in keV
    Fpeak_T = 200.  ! T(0)=20keV
    do i = 1,nr
     !T_i_rho(i)=T_ped*(Fpeak_T-1.)*(1.0-rho_hat(i)**2)**2+T_ped
      T_i_rho(i)=T_ped*(Fpeak_T-1.)*(1.0-rho_hat(i)**2)+T_ped
    enddo

  ! electron temperature in keV
    Fpeak_T = 200. 
    Fpeak_ei=1.0
    !  T(0)=20keV
    do i = 1,nr
     !T_e_rho(i)=Fpeak_ei*T_ped*(Fpeak_T-1./Fpeak_ei)*(1.0-rho_hat(i)**2)**2+T_ped
      T_e_rho(i)=Fpeak_ei*T_ped*(Fpeak_T-1./Fpeak_ei)*(1.0-rho_hat(i)**2)+T_ped
    enddo

  ! ion density in 10**19 1/m**3   n20 = 10.*n14     20= 10**20 1/m**3   14=10**14 1/cm**3
    Fpeak_n=1.
    do i = 1,nr
      n_i_rho(i)=10.*n14_ped*(Fpeak_n-1.)*(1.0-rho_hat(i)**2)+10.*n14_ped
    enddo

  ! electron desnisty in 10**19 1/m**3
    Fpeak_n=1.
    do i = 1,nr
      n_e_rho(i)=10.*n14_ped*(Fpeak_n-1.)*(1.0-rho_hat(i)**2)+10.*n14_ped
    enddo

  !equillibrium profiles

  !midplane minor radius in meters
    do i = 1,nr
      rmin_rho(i)=rmin*rho_hat(i)
    enddo


  !safety factor 
  !  do i = 1,nr
  !   q_rho(i)=q_0
  !   if(rho_hat(i).ge.rho_q) then
  !    q_rho(i)=q_0+(q_1-q_0)*(rho_hat(i)-rho_q)**2/(1.0-rho_q)**2
  !   endif
  !   if(rho_hat(i).lt.rho_q) then
  !    q_rho(i)=q_rho(i)*(1.+eps_q*rho_hat(i))
  !   endif
  !   if(rho_hat(i).ge.rho_q) then
  !    q_rho(i)=q_rho(i)*((1.+eps_q*rho_q)-(rho_hat(i)-rho_q)/(1.0-rho_q)*eps_q*rho_q)
  !   endif
  !  keeps s_hat > eps_q/(1+eps_q*rho_hat) for rho_hat<rho_q
  !  enddo
     do i = 1,nr
       q_rho(i)=q_0+(q_1-q_0)*rho_hat(i)**3
     enddo

  !Shafranov shifted major radius in meters
    do i = 1,nr
  !    Rmaj_rho(i)=Rmaj0*(1.0+delR0oa*(1.-rho_hat(i)))
  !    Rmaj_rho(i)=Rmaj0*(1.0+delR0oa*(1.-rho_hat(i)**2))
      Rmaj_rho(i)=Rmaj0*(1.0+delR0oa*(rmin/Rmaj0)*(1.-rho_hat(i)))
    enddo

  !elongation
    do i = 1,nr
  !   kappa_rho(i)=(kappa_1-kappa_0)*rho_hat(i)+kappa_0
      kappa_rho(i)=(kappa_1-kappa_0)*rho_hat(i)**4+kappa_0
    enddo

  !triangularity
    do i = 1,nr
  !   delta_rho(i)=(delta_1-delta_0)*rho_hat(i)+delta_0
      delta_rho(i)=(delta_1-delta_0)*rho_hat(i)**3+delta_0
    enddo

  ! compute beta_N_glob !should be about 2-3 x beta_N_ped with ITER beta_N_glob limited to ~2.5
    beta_N_glob = 0.
    avediv=0.
    do i = 1,nr
      beta_N_glob = beta_N_glob +kappa_rho(i)*rmin_rho(i)*Rmaj_rho(i)*n_i_rho(i)*(T_i_rho(i)+T_e_rho(i))
      avediv=avediv+kappa_rho(i)*rmin_rho(i)*Rmaj_rho(i)
    enddo
    beta_N_glob = beta_N_glob/avediv  !<n*T_i+n*T_e>
    beta_N_glob = beta_N_glob/(20.*n14_ped*T_ped)*beta_N_ped


  endif ! end i_ITER_profile=1

  if(i_ITER_profile .eq. 2) then
  !
  ! DIIID NBI effective E_alpha  test case
  !  main control parameters
  !  read in n_alpha_rho(i)  i=1,nr
  !  compute S0_rho(i) from n_alpha_rho(i)  slowing down formula
  !  guess  E_alpha (MeV) 
  !    effective NBI injection energy to match  an effective NBI T_alpha(i) temperature profile  
  !  Angioni formula dependent on  E_alpha  and T_e_rho(i) with a chi_eff(i)
  !  chi_eff(i)  used Q_fus = 1000000. and
  !
  !     chi_eff(i) =  (1.+5./Q_fus)*flux_source_rho(i)*(E_alpha*1000.)/n_i_rho(i)/ &
  !         ((-T_i_rho(i+1)+T_i_rho(i-1))/(2.*dr))
  !  
  !     flux_source_rho(i) follows from S0_rho(i) which follows from n_alpha_rho(i) input
  !         this ignores heating from ohmic and radiation loss...only NBI heating goes into chi_eff
  !  since slowing down is dtermined from T_e and n_e
  !   T_i_rho(i) set to  T_e_rho(i)
  !   n_i_rho(i) set to  n_e_rho(i)
  !


    print *, 'i_ITER_profile ==2  DIIID NBI effective E_alpha  test case'
    B_t = 2.05
    I_p = 1.0 
    rmin = 0.61 
    Rmaj0= 1.63

    delR0oa=0.13333   !Rmaj_center = 1.71

    q_1=8.6
    q_0=5.0

    kappa_1=1.67
    kappa_0=1.72

    delta_1=0.08
    delta_0=0.
   
    n14_ped=0.15

    T_ped = 0.15 
    ! T_i_rho(r_hat) = (1.5-0.15)*(1.-r_hat**2)+0.15  !F_peak_T = 10. 
    Fpeak_ei=1.75/1.5
    !  T_e_rho(r_hat) = (1.75-0.05)*(1.-r_hat**2)+0.05 


    arho = sqrt(kappa_1)*rmin !in meters

    !beta_ped_percent=beta_N_ped*I_p/B_t/rmin

    !T_ped = 0.12*beta_ped_percent*(B_t)**2/n14_ped

    beta_ped_percent = T_ped*n14_ped/(0.12*B_t**2)
    beta_N_ped=beta_ped_percent/(I_p/B_t/rmin)


  !plasma profiles

  ! ion temperature in keV
    Fpeak_T = 200.  ! T(0)=20keV
    do i = 1,nr
     !T_i_rho(i)=T_ped*(Fpeak_T-1.)*(1.0-rho_hat(i)**2)**2+T_ped
      T_i_rho(i)=T_ped*(Fpeak_T-1.)*(1.0-rho_hat(i)**2)+T_ped
    enddo

  ! electron temperature in keV
    Fpeak_T = 200.
    Fpeak_ei=1.0
    !  T(0)=20keV
    do i = 1,nr
     !T_e_rho(i)=Fpeak_ei*T_ped*(Fpeak_T-1./Fpeak_ei)*(1.0-rho_hat(i)**2)**2+T_ped
      T_e_rho(i)=Fpeak_ei*T_ped*(Fpeak_T-1./Fpeak_ei)*(1.0-rho_hat(i)**2)+T_ped
    enddo

  ! ion density in 10**19 1/m**3   n20 = 10.*n14     20= 10**20 1/m**3   14=10**14 1/cm**3
    Fpeak_n=1.
    do i = 1,nr
      n_i_rho(i)=10.*n14_ped*(Fpeak_n-1.)*(1.0-rho_hat(i)**2)+10.*n14_ped
    enddo

  ! electron desnisty in 10**19 1/m**3
    Fpeak_n=1.
    do i = 1,nr
      n_e_rho(i)=10.*n14_ped*(Fpeak_n-1.)*(1.0-rho_hat(i)**2)+10.*n14_ped
    enddo

  !equillibrium profiles

  !midplane minor radius in meters
    do i = 1,nr
      rmin_rho(i)=rmin*rho_hat(i)
    enddo


  !safety factor
  !  do i = 1,nr
  !   q_rho(i)=q_0
  !   if(rho_hat(i).ge.rho_q) then
  !    q_rho(i)=q_0+(q_1-q_0)*(rho_hat(i)-rho_q)**2/(1.0-rho_q)**2
  !   endif
  !   if(rho_hat(i).lt.rho_q) then
  !    q_rho(i)=q_rho(i)*(1.+eps_q*rho_hat(i))
  !   endif
  !   if(rho_hat(i).ge.rho_q) then
  !    q_rho(i)=q_rho(i)*((1.+eps_q*rho_q)-(rho_hat(i)-rho_q)/(1.0-rho_q)*eps_q*rho_q)
  !   endif
  !  keeps s_hat > eps_q/(1+eps_q*rho_hat) for rho_hat<rho_q
  !  enddo
     do i = 1,nr
       q_rho(i)=q_0+(q_1-q_0)*rho_hat(i)**3
     enddo

  !Shafranov shifted major radius in meters
    do i = 1,nr
  !    Rmaj_rho(i)=Rmaj0*(1.0+delR0oa*(1.-rho_hat(i)))
  !    Rmaj_rho(i)=Rmaj0*(1.0+delR0oa*(1.-rho_hat(i)**2))
      Rmaj_rho(i)=Rmaj0*(1.0+delR0oa*(rmin/Rmaj0)*(1.-rho_hat(i)))
    enddo

  !elongation
    do i = 1,nr
  !   kappa_rho(i)=(kappa_1-kappa_0)*rho_hat(i)+kappa_0
      kappa_rho(i)=(kappa_1-kappa_0)*rho_hat(i)**4+kappa_0
    enddo

  !triangularity
    do i = 1,nr
  !   delta_rho(i)=(delta_1-delta_0)*rho_hat(i)+delta_0
      delta_rho(i)=(delta_1-delta_0)*rho_hat(i)**3+delta_0
    enddo

  ! compute beta_N_glob !should be about 2-3 x beta_N_ped with ITER beta_N_glob limited to ~2.5
    beta_N_glob = 0.
    avediv=0.
    do i = 1,nr
     beta_N_glob = beta_N_glob +kappa_rho(i)*rmin_rho(i)*Rmaj_rho(i)*n_i_rho(i)*(T_i_rho(i)+T_e_rho(i))
     avediv=avediv+kappa_rho(i)*rmin_rho(i)*Rmaj_rho(i)
    enddo
    beta_N_glob = beta_N_glob/avediv  !<n*T_i+n*T_e>
    beta_N_glob = beta_N_glob/(20.*n14_ped*T_ped)*beta_N_ped


  endif ! end i_ITER_profile=2

!9.11.14
! overwrite NBI and plasma profiles
  if(NBI_flag) then
    inquire(file='out.NBIprint',exist=iexist)

    if(iexist) then
      open(unit=2,file='out.NBIprint',status='old')

    !  write(2,*) ' r_hat r_minor'
      read(2,*) 
      read(2,*)
      read(2,*) 

      do i=1,nr
        read(2,10) dummy
      enddo

    !  write(2,*) ' ni_2 10**19/m**3'
      read(2,*)
      read(2,*)
      read(2,*)

      do i=1,nr
        read(2,10) n_alpha_rho(i)  !should not be overriden if NBI_flag = 1
      enddo
    !   n_alpha_rho(i) is really n_NBI_rho(i) slowing down density

    !  write(2,*) ' -dni_2/dr 10**19/m**4'
      read(2,*)
      read(2,*)
      read(2,*)

      do i=1,nr
        read(2,10) dummy
      enddo


    !  write(2,*) ' Ti_2 keV'
      read(2,*) 
      read(2,*)
      read(2,*)
      do i=1,nr
        read(2,10) T_alpha_equiv_rho(i)   !EP or T_alpha_equiv_rho(i) is recalculated
      enddo

    !  write(2,*) ' ni 10**19/m**3'
      read(2,*)
      read(2,*)
      read(2,*)
      do i=1,nr
        read(2,10) n_i_rho(i)
      enddo

    !  write(2,*) ' Ti keV'
      read(2,*)
      read(2,*) 
      read(2,*)
      do i=1,nr
        read(2,10) T_i_rho(i)
      enddo

    !  write(2,*) ' ne 10**19/m**3'
      read(2,*)
      read(2,*)
      read(2,*)
      do i=1,nr
        read(2,10) n_e_rho(i)
      enddo

    !  write(2,*) ' Te keV'
      read(2,*)
      read(2,*) 
      read(2,*) 
      do i=1,nr
        read(2,10) T_e_rho(i)
      enddo

      close(2)
      
      beta_N_glob =2. !override
      beta_N_ped = .5 !override
      beta_ped_percent = 1.0 !override
      ! print *, 'beta_N_glob, beta_N_ped, beta_ped_percent not meaningful for NBI_flag=1'

    else

      print *, 'out.NBIprint file is not found'
      stop

    endif

  endif !NBI_flag = 1

10  format(ES14.7)

end subroutine EPtran_comp_eq_plasma