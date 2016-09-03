!-------------------------------------------------------------------------
! EPtran_use_parameter.f90
!
! PURPOSE:
!   controls, inputs, outputs and other parameters
!
!-------------------------------------------------------------------------

module EPtran_use_parameter

  implicit none

    !!CONTROL PARAMETERS
    real :: dt = 1.e-5       !time step
    integer :: ntstep = 4e5 !number of time steps
    logical :: i_new_start = .true.  !flag for restart

    logical :: NBI_flag = .false. !NBI_flag

    logical :: i_flux_r = .true.
    logical :: i_flux_E = .false.
    integer :: i_Model = 2
    logical :: i_AE = .false.
    logical :: i_pinch = .false.
    logical :: i_thermal_pinch_off = .false.
    logical :: i_BC_r1 = .true.
    real :: delta1 = 0

    real :: D_factor = 1.0
    real :: A_factor = 1.0

    real :: m_alpha = 4.0  !mass of the energetic particle
    real :: z_alpha = 2.0  !charge of the energetic particle

    !!INPUT PARAMETERS
    integer, parameter :: nr = 51     !the grid number of rho

    !defaults
    !Kinsey et al Standard case for ITER:
    !J.E. Kinsey, G.M. Staebler, R.E. Waltz, and R. V. Budny, Nucl. Fusion 51 92012) 083001
    !see Figs 7,11,13
    !defaults
    real :: beta_N_ped = 0.92
    real :: n14_ped = 0.9
    real :: I_p = 15
    real :: B_t = 5.3
    real :: Rmaj0 = 6.2
    real :: delR0oa = 0.3
    real :: rmin = 2.0
    real :: kappa_1 = 1.75
    real :: kappa_0 = 1.0
    real :: delta_1 = 0.2
    real :: delta_0 = 0.0
    real :: q_1 = 3.8
    real :: q_0 = 1.0
    real :: rho_q = 0.5
    real :: eps_q = 0.1
    real :: Zeff = 1.7
    real :: M_DT = 2.5
    real :: E_alpha = 3.5
    real :: Fpeak_T = 3.0
    real :: Fpeak_n = 1.1
    real :: Fpeak_ei = 1.1
    integer :: use_t_equiv = 0

    !!OUTPUT PARAMETERS
    real :: beta_ped_percent
    real :: T_ped
    real :: beta_N_glob
    real :: arho

    real :: rho_hat(nr)
    real :: q_rho(nr)
    real :: Rmaj_rho(nr)
    real :: rmin_rho(nr)
    real :: kappa_rho(nr)
    real :: delta_rho(nr)
    real :: T_i_rho(nr)
    real :: T_e_rho(nr)
    real :: n_i_rho(nr)
    real :: n_e_rho(nr)
    real :: n_alpha_rho(nr)
    real :: T_alpha_equiv_rho(nr)
    real :: E_c_hat_rho(nr)
    real :: E_alpha_grid(nr)

    !!OTHER PARAMETERS
    real, parameter :: pi = 3.141592565

    real, dimension(nr) ::tau_ee

    real, dimension(nr) :: n_alpha_ave_rho
    real, dimension(nr) :: S0_rho
    real, dimension(nr) :: bannana_rho
    real, dimension(nr) :: tau_s_rho

    contains

    function rdiff(rho) !the radial diff function
      !diff for radial, rd(i) = -d() = -(rho(i+1)-rho(i-1))/2.

      implicit none

      real,dimension(nr) :: rho,rdiff

      rdiff(1) = 0.
      rdiff(2:nr-1) = -(rho(3:nr) - rho(1:nr-2))/2.
      rdiff(nr) = -(rho(nr) - rho(nr-1))

    end function rdiff

end module EPtran_use_parameter
