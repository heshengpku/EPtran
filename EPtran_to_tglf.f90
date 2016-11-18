!---------------------------------------------------------
! EPtran_to_tglf.f90
!
! PURPOSE:
!  compute the needed parameters for DEP and TGLF
!
! 6.11.14
!---------------------------------------------------------
module EPtran_to_tglf
  
  implicit none
  
  integer,parameter :: nn = 51    !must be same as nr

  !species
  real :: mi_tglf,zeff_tglf,m_alpha_tglf,z_alpha_tglf
  !geometry
  real,dimension(nn) :: r_hat_tglf_r,rmaj_hat_tglf_r,drmajdr_tglf_r,&
                        kappa_tglf_r,s_kappa_tglf_r,delta_tglf_r,s_delta_tglf_r,&
                        q_tglf_r,s_hat_tglf_r,q_prime_tglf_r,p_prime_tglf_r
  !profiles
  real,dimension(nn) :: ni_hat_tglf_r,Ti_hat_tglf_r,n_EP_hat_r,T_EP_hat_r,&
                        aoLn_e_tglf_r,aoLT_e_tglf_r,aoLn_i_tglf_r,aoLT_i_tglf_r,&
                        aoLn_EP_r,aoLT_EP_r,betae_unit_tglf_r,xnue_tglf_r!,&
                        !vexb_shear_tglf_r,vpar_shear_e_tglf_r,vpar_shear_i_tglf_r

  !not TGLF input parameters
  real,dimension(nn) :: kymark,B_unit,beta_unit_r,p_r,dlnpdr_r,c_s,rho_s,rho_star,chi_gB
  real,dimension(nn) :: omega_TAE,omega_plus,omega_minus
  
end module EPtran_to_tglf



subroutine get_tglf_parameters
  
  use EPtran_use_parameter
!INPUTS                  !nr,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
!                        !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
!                        !Zeff,M_DT
!OUTPUTS                 !beta_ped_percent,T_ped,beta_N_glob,arho,
!                        !rho_hat(i),q_rho(i),Rmaj_rho(i),rminnho(i),
!                        !kappa_rho(i),delta_rho(i),T_i_rho(i),T_e_rho(i),
!                        !n_i_rho(i),n_e_rho(i),n_alpha_rho(i),T_alpha_equiv_rho(i),
!OTHERS                  !pi,tau_ee(i)

  use EPtran_to_tglf
  
  implicit none
  
!local parameter
  integer :: i
  integer :: n = 3 !toroidal number
  real :: dr
  logical :: dump_flag = .true.
  
  dr = rmin/real(nr-1)
  
  !ion mass in uint of proton mass
  mi_tglf = M_DT/2.0

  zeff_tglf = Zeff

  m_alpha_tglf = m_alpha/2.0
  z_alpha_tglf = z_alpha

  !minor radius normalized by rmin
  r_hat_tglf_r(:) = rho_hat(:)
  
  !normalized centroid major radius
  rmaj_hat_tglf_r(:) = Rmaj_rho(:)/rmin

  !drmaj/dr
  i = 1
  drmajdr_tglf_r(i) = (Rmaj_rho(i+1)-Rmaj_rho(i))/dr
  do i = 2,nr-1
    drmajdr_tglf_r(i) = (Rmaj_rho(i+1)-Rmaj_rho(i-1))/2./dr
  enddo
  i = nr
  drmajdr_tglf_r(i) = (Rmaj_rho(i)-Rmaj_rho(i-1))/dr

  !elongation of flux surface
  kappa_tglf_r(:) = kappa_rho(:)

  !shear in elongation : r/kappa*dkappa/dr
  i = 1
  s_kappa_tglf_r(i) = rmin*rho_hat(i)/kappa_rho(i)*(kappa_rho(i+1)-kappa_rho(i))/dr
  do i = 2,nr-1
    s_kappa_tglf_r(i) = rmin*rho_hat(i)/kappa_rho(i)*(kappa_rho(i+1)-kappa_rho(i-1))/2./dr
  enddo
  i = nr
  s_kappa_tglf_r(i) = rmin*rho_hat(i)/kappa_rho(i)*(kappa_rho(i)-kappa_rho(i-1))/dr

  !triangularity of flux surface
  delta_tglf_r(:) = delta_rho(:)

  !shear in trangularity : r/sqrt(1-delta**2)*ddelta/dr
  i = 1
  s_delta_tglf_r(i) = rmin*rho_hat(i)/sqrt(1.-delta_rho(i)**2)*(delta_rho(i+1)-delta_rho(i))/dr
  do i = 2,nr-1
    !s_delta_tglf_r(i) = rho_hat(i)/sqrt(1.-delta_rho(i)**2)*(delta_1-delta_0)
    s_delta_tglf_r(i) = rmin*rho_hat(i)/sqrt(1.-delta_rho(i)**2)*(delta_rho(i+1)-delta_rho(i-1))/2./dr
  enddo
  i = nr
  s_delta_tglf_r(i) = rmin*rho_hat(i)/sqrt(1.-delta_rho(i)**2)*(delta_rho(i)-delta_rho(i-1))/dr

  !safety factor
  q_tglf_r(:) = q_rho(:)

  !magnetic shear : r/q*dq/dr
  !s_hat_tglf_r(1) = 0.
  i = 1
  s_hat_tglf_r(i) = rmin*rho_hat(i)/q_rho(i)*(q_rho(i+1)-q_rho(i))/dr
  do i = 2,nr-1
    s_hat_tglf_r(i) = rmin*rho_hat(i)/q_rho(i)*(q_rho(i+1)-q_rho(i-1))/2./dr
  enddo
  i = nr
  s_hat_tglf_r(i) = rmin*rho_hat(i)/q_rho(i)*(q_rho(i)-q_rho(i-1))/dr
  
  !pressure profile : p = n_e*T_e + n_i*T_i + n_EP*T_EP
  !the total pressure including EP species 7.18.2016
  p_r(:) = n_e_rho(:)*T_e_rho(:) + n_i_rho(:)*T_i_rho(:) + n_alpha_rho(:)*T_alpha_equiv_rho(:)

  !pressure gradient profile : dlnp/dr = dp/dr/p, it is usually negative 7.18.2016
  dlnpdr_r(1) = 0.
  do i = 2,nr-1
    dlnpdr_r(i) = (p_r(i+1)-p_r(i-1))/2./dr/p_r(i)
  enddo
  i = nr
  dlnpdr_r(i) = (p_r(i)-p_r(i-1))/dr/p_r(i)

  !B_unit .= B_t*kappa
  B_unit(:) = B_t*kappa_rho(:)
  
  !beta profile : beta_unit(:)  = 8*pi*p(:)/B_unit**2 in cgs
  !                             = 2*mu0*p(:)/B_unit**2 in SI
  !mu0 = 4*pi*1e-7, e = 1.6022e-19 in SI
  !p(:) unit [10**19/m**3]*[keV] --> SI = [J/m**3]
  beta_unit_r(:) = 2.*(4.*pi*1.e-7)*1.6022*1.e3*p_r(:)/B_unit(:)**2
  
  !normalized dq/dpsi : q/r*dq/dr*rmin**2 = (q/r_hat)**2*s_hat
  q_prime_tglf_r(:) = (q_rho(:)/rho_hat(:))**2*s_hat_tglf_r(:)
  
  !normalized dp/dpsi : q/r*dp/dr*rmin**2*8*pi/B_unit**2 = q/r_hat*dlnp/dr*rmin*beta_N
  !looking for tgyro_tglf_map.f90
  p_prime_tglf_r(:) = q_rho(:)/rho_hat(:)*(beta_unit_r(:)/(8*pi))*(rmin*dlnpdr_r(:))
  
  !ion density normalized by electron density
  ni_hat_tglf_r(:) = n_i_rho(:)/n_e_rho(:)

  !ion temperature normalized by electron temperature
  Ti_hat_tglf_r(:) = T_i_rho(:)/T_e_rho(:)
  
  !normalized electron density gradient : -rmin/ne*dne/dr
  aoLn_e_tglf_r(:) = rmin/n_e_rho(:)*rdiff(n_e_rho(:))/dr

  !normalized electron temperature gradient : -rmin/Te*dTe/dr
  aoLT_e_tglf_r(:) = rmin/T_e_rho(:)*rdiff(T_e_rho(:))/dr

  !normalized ion density gradient : -rmin/ni*dni/dr
  aoLn_i_tglf_r(:) = rmin/n_i_rho(:)*rdiff(n_i_rho(:))/dr

  !normalized ion temperature gradient : -rmin/Ti*dTi/dr
  aoLT_i_tglf_r(:) = rmin/T_i_rho(:)*rdiff(T_i_rho(:))/dr

  !electron beta profile
  betae_unit_tglf_r(:) = beta_unit_r(:)*n_e_rho(:)*T_e_rho(:)/p_r(:)
  
  !cs0 = sqrt(T_e/m_D) = sqrt(T_e[keV]/(2.*938MeV/c**2))
  c_s(:) = 3.0e5*sqrt(T_e_rho(:)/2./0.938) !in m/s

  !gamma_eb(:) = -rho_hat(:)*rmin/q_rho(:)*w0p(:)

  !gamma_p(:)  = -Rmaj_rho(:)*w0p(:)
  
  !electron collision frequency nuei/(cs0/rmin), 4.20.2016, tau_ei = m_i/m_e*tau_ee
  xnue_tglf_r(:) = 1./tau_ee(:)/M_DT/1836.*rmin/c_s(:)

  rho_s(:) = c_s(:)/(9.577e7/2.*B_unit(:))

  rho_star(:) = rho_s(:)/rmin  !4.19.2016

  !chi_gB in m**2/sec
  chi_gB(:) = rho_s(:)**2*c_s(:)/rmin

  !normalized energitic particle temperature for DEP
  T_EP_hat_r(:) = T_alpha_equiv_rho(:)/T_e_rho(:)

  !add, 4.19.2016
  !normalized EP density
  n_EP_hat_r(:) = n_alpha_rho(:)/n_e_rho(:)

  !normalized EP density gradient: -rmin/n_EP*dn_EP/dr
  aoLn_EP_r(:) = rmin/n_alpha_rho(:)*rdiff(n_alpha_rho(:))/dr

  !normalized EP temperature gradient : -rmin/T_EP*dT_EP/dr
  aoLT_EP_r(:) = rmin/T_alpha_equiv_rho(:)*rdiff(T_alpha_equiv_rho(:))/dr

  !ky = (nq/r)rho_s = (nq/(r/a))rho_star
  kymark = n*q_tglf_r/r_hat_tglf_r*rho_star

  !omega_TAE/(c_s/a) = sqrt(2/beta_e/ni_hat)/2/q/(R/a)
  omega_TAE = (2./betae_unit_tglf_r/ni_hat_tglf_r)**0.5/2./q_tglf_r/rmaj_hat_tglf_r

  !omega_plus[minus] = omega_TAE/sqrt(1-[+]2r/R)
  omega_plus  = omega_TAE / (1. - 2.*r_hat_tglf_r/rmaj_hat_tglf_r)**0.5
  omega_minus = omega_TAE / (1. + 2.*r_hat_tglf_r/rmaj_hat_tglf_r)**0.5

  if(dump_flag) then
    open(unit=2,file='EPtran_to_tglf.out',status='replace')
    !just as TGLFEP input.profile format, 11.17.2016
    write(2,*)-1.0,'  SIGN_BT'
    write(2,*)-1.0,'  SIGN_IT'
    write(2,*)nr-1,'  NR'
    write(2,*)3,'  NS'
    write(2,*)1,'  GEOMETRY_FLAG'
    write(2,*)'--------------------------------------------------------------'
    write(2,*)'# electron species '
    write(2,*)-1.0,'  ZS'
    write(2,*)1./1836./2.,'  MASS m/m_D'
    write(2,*)'# normalized density gradients: rlns'
    write(2,10)aoLn_e_tglf_r(2:nr)
    write(2,*)'# normalized temperature gradients: rlts'
    write(2,10)aoLT_e_tglf_r(2:nr)
    write(2,*)'--------------------------------------------------------------'
    write(2,*)'# ion species ',1
    write(2,*)1.0,'  ZS' 
    write(2,*)mi_tglf,'  MASS m/m_D'
    write(2,*)'# normalized density: as'
    write(2,10)ni_hat_tglf_r(2:nr)
    write(2,*)'# normalized temperature: taus'
    write(2,10)Ti_hat_tglf_r(2:nr)
    write(2,*)'# normalized density gradients: rlns'
    write(2,10)aoLn_i_tglf_r(2:nr)
    write(2,*)'# normalized temperature gradients: rlts'
    write(2,10)aoLT_i_tglf_r(2:nr)
    write(2,*)'--------------------------------------------------------------'
    write(2,*)'# ion species ',2
    write(2,*)z_alpha_tglf,'  ZS' 
    write(2,*)m_alpha_tglf,'  MASS m/m_D'
    write(2,*)'# normalized density: as'
    write(2,10)n_EP_hat_r(2:nr)
    write(2,*)'# normalized temperature: taus'
    write(2,10)T_EP_hat_r(2:nr)
    write(2,*)'# normalized density gradients: rlns'
    write(2,10)aoLn_EP_r(2:nr)
    write(2,*)'# normalized temperature gradients: rlts'
    write(2,10)aoLT_EP_r(2:nr)
    write(2,*)'--------------------------------------------------------------'
    write(2,*)'# Geometry '
    write(2,*)'# minor radius: rmin'
    write(2,10)r_hat_tglf_r(2:nr)
    write(2,*)'# major radius: rmaj'
    write(2,10)rmaj_hat_tglf_r(2:nr)
    write(2,*)'# safety factor: q'
    write(2,10)q_tglf_r(2:nr)
    write(2,*)'# magnetic shear: shear'
    write(2,10)s_hat_tglf_r(2:nr)
    write(2,*)'# q_prime'
    write(2,10)q_prime_tglf_r(2:nr)
    write(2,*)'# p_prime'
    write(2,10)p_prime_tglf_r(2:nr)
    write(2,*)'# shift'
    write(2,10)drmajdr_tglf_r(2:nr)
    write(2,*)'# elogation: kappa'
    write(2,10)kappa_tglf_r(2:nr)
    write(2,*)'# shear in elogation: s_kappa'
    write(2,10)s_kappa_tglf_r(2:nr)
    write(2,*)'# triangularity: delta'
    write(2,10)delta_tglf_r(2:nr)
    write(2,*)'# shear in triangularity: s_delta'
    write(2,10)s_delta_tglf_r(2:nr)
    write(2,*)'# squareness: zeta'
    do i = 2,nr
      write(2,10)0.0
    enddo
    write(2,*)'# shear in squareness: s_zeta'
    do i = 2,nr
      write(2,10)0.0
    enddo
    write(2,*)'--------------------------------------------------------------'
    write(2,*)'# effective ion charge: zeff'
    do i = 2,nr
      write(2,10)zeff_tglf
    enddo
    write(2,*)'# betae'
    write(2,10)betae_unit_tglf_r(2:nr)
    write(2,*)'# rho_star = rho_s/a'
    write(2,10)rho_star(2:nr)
    write(2,*)'# omega_TAE / (c_s/a)'
    write(2,10)omega_TAE(2:nr)

    ! do i = 1,nr
    !   write(2,*)'--------------------------------------------------------------'
    !   write(2,20)'r_hat',r_hat_tglf_r(i)
    !   write(2,20)'rmaj_hat',rmaj_hat_tglf_r(i)
    !   write(2,20)'drmajdr',drmajdr_tglf_r(i)
    !   write(2,20)'kappa',kappa_tglf_r(i)
    !   write(2,20)'s_kappa',s_kappa_tglf_r(i)
    !   write(2,20)'delta',delta_tglf_r(i)
    !   write(2,20)'s_delta',s_delta_tglf_r(i) 
    !   write(2,20)'q',q_tglf_r(i)
    !   write(2,20)'s_hat',s_hat_tglf_r(i) 
    !   write(2,20)'q_prime',q_prime_tglf_r(i)  
    !   write(2,20)'p_prime',p_prime_tglf_r(i)
    !   write(2,20)'aoLn_e',aoLn_e_tglf_r(i)
    !   write(2,20)'aoLT_e',aoLT_e_tglf_r(i)
    !   write(2,20)'ni_hat',ni_hat_tglf_r(i)
    !   write(2,20)'Ti_hat',Ti_hat_tglf_r(i)
    !   write(2,20)'aoLn_i',aoLn_i_tglf_r(i)
    !   write(2,20)'aoLT_i',aoLT_i_tglf_r(i)
    !   write(2,20)'n_EP_hat',n_EP_hat_r(i)
    !   write(2,20)'T_EP_hat',T_EP_hat_r(i)
    !   write(2,20)'aoLn_EP',aoLn_EP_r(i)
    !   write(2,20)'aoLT_EP',aoLT_EP_r(i)
    !   write(2,20)'betae_unit',betae_unit_tglf_r(i)
    !   write(2,20)'xnue',xnue_tglf_r(i)
    !   write(2,20)'rho_star',rho_star(i)
    !   write(2,20)'ky_mark',kymark(i)
    !   write(2,20)'omega_TAE/(c_s/a)',omega_TAE(i)
    !   write(2,20)'omega_plus',omega_plus(i)
    !   write(2,20)'omega_minus',omega_minus(i)
    !   write(2,20)'chi_gB in m**2/sec',chi_gB(i)
    !   write(2,*)'--------------------------------------------------------------'
    ! enddo

    close(2)

  endif !dump

10  format(F11.6)
20  format(T2,A,T23,':' ,F11.6)

end subroutine get_tglf_parameters
