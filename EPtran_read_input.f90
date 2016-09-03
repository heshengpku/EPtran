!---------------------------------------------------------
! EPtran_read_input.f90
!
! PURPOSE:
!  Reads all input parameters from input: EPtran_input
!
!---------------------------------------------------------

subroutine EPtran_read_input

    use EPtran_use_parameter
!INPUTS	                 !nr,beta_N_ped,n14_ped,I_p,B_t,Rmaj0,delR0oa,rmin,
!                        !kappa_1,kappa_0,delta_1,delta_0,q_1,q_0,rho_q,eps_q,
!                        !Zeff,M_DT,E_alpha,Fpeak_T,Fpeak_n,Fpeak_ei,NBI_flag

  !--------------------------------------
  implicit none
  !
  logical :: i_read_default,iexist
  !
  character(len=80) :: comment
  !--------------------------------------

  i_read_default = .true.

  if(i_read_default) then

    inquire(file='EPtran_input',exist=iexist)

    if(iexist) then
      open(unit=1,file='EPtran_input',status='old')

      !----------------------------------------------------------
      ! Order and variable format in EPtran_input  must match here.
      !

      !read(1,*) nr
      read(1,*) beta_N_ped
      read(1,*) n14_ped
      read(1,*) I_p
      read(1,*) B_t
      read(1,*) Rmaj0
      read(1,*) delR0oa
      read(1,*) rmin
      read(1,*) kappa_1
      read(1,*) kappa_0
      read(1,*) delta_1
      read(1,*) delta_0
      read(1,*) q_1
      read(1,*) q_0
      read(1,*) rho_q
      read(1,*) eps_q
      read(1,*) Zeff
      read(1,*) M_DT
      read(1,*) E_alpha
      read(1,*) Fpeak_T
      read(1,*) Fpeak_n
      read(1,*) Fpeak_ei
      read(1,'(I1)') use_t_equiv
      read(1,*) NBI_flag  !false=Alpha   true=NBI

      close(1)
    
    else

      print *, 'EPtran_input file is not found'
      stop

    endif

  endif


end subroutine EPtran_read_input
