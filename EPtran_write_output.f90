!---------------------------------------------------------
! EPtran_write_output.f90
!
! PURPOSE:
!  Writes all output parameters to output. EPtran.out 
!
!---------------------------------------------------------

subroutine EPtran_write_output

  use EPtran_use_parameter
  !INPUTS	                 !nr,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
  !                        !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
  !                        !Zeff,M_DT,E_alpha,Fpeak_T,Fpeak_n,Fpeak_ei
  !OUTPUTS                 !beta_ped_percent,T_ped,beta_N_glob,arho,
  !                        !rho_hat(i),q_rho(i),Rmaj_rho(i),rmin_rho(i),
  !                        !kappa_rho(i),delta_rho(i),T_i_rho(i),T_e_rho(i),
  !                        !n_i_rho(i),n_e_rho(i),n_alpha_rho(i),T_alpha_equiv_rho(i),
  !                        !E_c_hat_rho(i)
  !OTHER                   !n_alpha_ave_rho(i),bannana_rho(i)


  !--------------------------------------
  implicit none
  !
  integer :: i
  !
  character(len=80) :: comment
  !--------------------------------------

  !----------------------------------------------------------
  ! Order and variable format in Alpha_input  must match here.
  !

  open(unit=2,file='EPtran.out',status='replace')

! write inputs
  write(2,*) 'INPUTS'
  write(2,*) '--------------------------------------------------------------'

  write(2,*) 'n_rho_grid=',nr
  write(2,*) 'beta_N_ped=',beta_N_ped
  write(2,*) 'n14_ped=',n14_ped, ' 10**14 1/cm**3'
  write(2,*) 'I_p=',I_p,' MA'
  write(2,*) 'B_t=',B_t,' Tesla'
  write(2,*) 'Rmaj0=',Rmaj0,' meters'
  write(2,*) 'delR0oa=',delR0oa
  write(2,*) 'rmin=',rmin, ' meters'
  write(2,*) 'kappa_1=',kappa_1
  write(2,*) 'kappa_0=',kappa_0
  write(2,*) 'delta_1=',delta_1
  write(2,*) 'delta_0=',delta_0
  write(2,*) 'q_1=',q_1
  write(2,*) 'q_0=',q_0
  write(2,*) 'rho_q=',rho_q
  write(2,*) 'epa_q=',eps_q
  write(2,*) 'Zeff=',Zeff
  write(2,*) 'M_DT=',M_DT
  write(2,*) 'E_apha=',E_alpha,' meV'
  write(2,*) 'Fpeak_T=',Fpeak_T
  write(2,*) 'Fpeak_n=',Fpeak_n
  write(2,*) 'Fpeak_ei=',Fpeak_ei
  write(2,'(A,I1)') 'use_T_equiv=', use_t_equiv
  write(2,*) 'NBI_flag=', NBI_flag
  write(2,*) '--------------------------------------------------------------'

!write outputs
  write(2,*) 'OUTPUTS'
  write(2,*) '--------------------------------------------------------------'
!beta_ped_percent,T_ped,beta_N_glob,arho
  if(.not. NBI_flag) then
    write(2,*) 'beta_ped_percent=',beta_ped_percent
    write(2,*) 'T_ped=',T_ped,' keV'
    write(2,*) 'beta_N_glob=',beta_N_glob
  endif
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'OUTPUT profiles versus rho_hat for plotting'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'rho_hat grid'
  do i = 1,nr
    write(2,10) rho_hat(i) 
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'plasma profile'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'ion temperature in keV'
  do i = 1,nr
    write(2,10) T_i_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'electron temperature in keV'
  do i = 1,nr
    write(2,10) T_e_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'ion density in 10**19 1/m**3'
  do i = 1,nr
    write(2,10) n_i_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'electron density in 10**19 1/m**3'
  do i = 1,nr
    write(2,10) n_e_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'equilibrium profiles'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'midplane minor radius in meters'
  do i = 1,nr
    write(2,10) rmin_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'safety factor' 
  do i = 1,nr
    write(2,10) q_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'Shafranov shifted major radius in meters'
  do i = 1,nr
    write(2,10) Rmaj_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'elongation'
  do i = 1,nr
    write(2,10) kappa_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'triangularity'
  do i = 1,nr
    write(2,10) delta_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'slowing down alpha profiles'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'alpha density in 10**19 1/m**3' 
  do i = 1,nr
    write(2,10) n_alpha_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'alpha equivalent Maxwellian temperature in keV'
  do i = 1,nr
    write(2,10) T_alpha_equiv_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'E_alpha normed E_c: E_c_hat=E_c/E_alpha   cross over energy'
  do i = 1,nr
    write(2,10) E_c_hat_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'

  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'orbit smeared slowing down alpha profiles'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'smeared alpha density in 10**19 1/m**3'
  do i = 1,nr
    write(2,10) n_alpha_ave_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'

  write(2,*) '--------------------------------------------------------'
  write(2,*) 'gradient slowing down alpha profiles'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'alpha density gradient in 10**19 1/m**3 per m'
  write(2,*) 0.0
  do i = 2,nr-1
    write(2,10) -(n_alpha_rho(i+1)-n_alpha_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
  i=nr
  write(2,10) -(n_alpha_rho(i)-n_alpha_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))
  write(2,*) '--------------------------------------------------------------'

  write(2,*) '--------------------------------------------------------'
  write(2,*) 'normed inverse temerpature gradient length'
  write(2,*) ' slowing down alpha profiles:  a/LT_alpha'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 0.0
  do i = 2,nr-1
    write(2,10) -(T_alpha_equiv_rho(i+1)-T_alpha_equiv_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))/&
                  T_alpha_equiv_rho(i)*rmin
  enddo
  i=nr
  write(2,10) -(T_alpha_equiv_rho(i)-T_alpha_equiv_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))/&
                  T_alpha_equiv_rho(i)*rmin
  write(2,*) '--------------------------------------------------------------'

  write(2,*) '--------------------------------------------------------'
  write(2,*) 'normed inverse density gradient length' 
  write(2,*) ' slowing down alpha profiles:  a/Ln_alpha'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 0.0
  do i = 2,nr-1
    write(2,10) -(n_alpha_rho(i+1)-n_alpha_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))/&
                  n_alpha_ave_rho(i)*rmin
  enddo
  i=nr
  write(2,10) -(n_alpha_rho(i)-n_alpha_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))/&
                  n_alpha_ave_rho(i)*rmin
  write(2,*) '--------------------------------------------------------------'

  write(2,*) '--------------------------------------------------------'
  write(2,*) 'gradient orbit smeared slowing down alpha profiles'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 'smeared alpha density gradient in 10**19 1/m**3 per m'
  write(2,10) 0.0
  do i = 2,nr-1
    write(2,10) -(n_alpha_ave_rho(i+1)-n_alpha_ave_rho(i-1))/rmin/(rho_hat(i+1)-rho_hat(i-1))
  enddo
  i=nr
  write(2,10) -(n_alpha_ave_rho(i)-n_alpha_ave_rho(i-1))/rmin/(rho_hat(i)-rho_hat(i-1))
  write(2,*) '--------------------------------------------------------------'

  if(.not. NBI_flag) then

    write(2,*) '--------------------------------------------------------'
    write(2,*) 'slowing down alpha beta-percent profile'
    write(2,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(2,10) beta_ped_percent*n_alpha_rho(i)*T_alpha_equiv_rho(i)/&
                  (n_i_rho(nr)*T_i_rho(nr)+n_e_rho(nr)*T_e_rho(nr))
    enddo
    write(2,*) '--------------------------------------------------------------'

    write(2,*) '--------------------------------------------------------'
    write(2,*) 'plasma beta-percent profile'
    write(2,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(2,10) beta_ped_percent*(n_i_rho(i)*T_i_rho(i)+n_e_rho(i)*T_e_rho(i))/&
                  (n_i_rho(nr)*T_i_rho(nr)+n_e_rho(nr)*T_e_rho(nr))
    enddo
    write(2,*) '--------------------------------------------------------------'

    write(2,*) '--------------------------------------------------------'
    write(2,*) 'slowing down alpha density profile used by Gorelenko(NF:2003) & Chen(PoP:2010)' 
    write(2,*) '(1.-rho_jat(i)**2)**(7./2.)'
    write(2,*) 'density  in 10**19 1/m**3'
    write(2,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(2,10)  n_alpha_rho(1)*(1.-rho_hat(i)**2)**(7./2.)
    enddo
    write(2,*) '--------------------------------------------------------------'

    write(2,*) '--------------------------------------------------------'
    write(2,*) 'slowing down alpha density gradient  profile used by Gorelenko(NF:2003) & Chen(PoP:2010)'
    write(2,*) '(1.-rho_jat(i)**2)**(7./2.)'
    write(2,*) 'density gradient in 10**19 1/m**3 per m'
    write(2,*) '--------------------------------------------------------------'
    do i = 1,nr
      write(2,10)  n_alpha_rho(1)*(1.-rho_hat(i)**2)**(5./2.)*7./2.*2.*rho_hat(i)/rmin
    enddo
    write(2,*) '--------------------------------------------------------------'
  endif !NBI_flag = 0

  write(2,*) '--------------------------------------------------------'
  write(2,*) 'normed log_e gradient slowing down alpha profiles'
  write(2,*) '--------------------------------------------------------------'
  write(2,*) 0.0
  do i = 2,nr-1
    write(2,10) -(n_alpha_rho(i+1)-n_alpha_rho(i-1))/(rho_hat(i+1)-rho_hat(i-1))/n_alpha_rho(i)
  enddo
  i=nr
  write(2,10) -(n_alpha_rho(i)-n_alpha_rho(i-1))/(rho_hat(i)-rho_hat(i-1))/n_alpha_rho(i)
  write(2,*) '--------------------------------------------------------------'

  write(2,*) '--------------------------------------------------------'
  write(2,*) ' slowing down alpha profiles normed to local electron density'
  write(2,*) '--------------------------------------------------------------'
  do i = 1,nr
    write(2,10) n_alpha_rho(i)/n_e_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'


  write(2,*) '--------------------------------------------------------'
  write(2,*) ' bannana half width relative to r_minor'
  write(2,*) '--------------------------------------------------------------'
  do i = 1,nr
    write(2,10) bannana_rho(i)
  enddo
  write(2,*) '--------------------------------------------------------------'











  close(2)

10  format(ES14.7)

end subroutine EPtran_write_output
