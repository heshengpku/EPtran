!-------------------------------------------------------------------------
! EPtran_use_transport.f90
!
! PURPOSE:
!   The velocity space grids for EPtran_transport.f90 
!
!-------------------------------------------------------------------------

module EPtran_use_transport


  implicit none

  integer, parameter :: nEa = 10    !add grid number for E_hat > 1
  integer, parameter :: nE = 1+500+nEa  !number of E_hat grid
  
  real, parameter :: dE = 1./(nE-nEa-1)
  real, dimension(nE) :: E_hat

  integer, parameter :: ny = 1 !y=sqrt(1-lambda) or lambda=1-y**2
  real, parameter :: dy = 1./(ny+1)
  real, dimension(ny) :: y

  contains
  !The integral functions (for E and lambda)

  real function Int_y(f)
    !Integral for y, Int f(y) from 0 to 1

    implicit none

    real,dimension(ny),intent(in) :: f

    integer :: k

    Int_y = 0.5*f(1)*dy
    do k = 1,ny
       Int_y = Int_y+f(k)*dy
    enddo
    Int_y = Int_y+0.5*f(ny)*dy

  end function Int_y

  real function Int_E(f,N_grid)
    !Integral for E, Int(f*2*pi*E_hat**0.5)dE_hat from 0 to E_c at N_grid for energy

    implicit none

    real,dimension(nE),intent(in) :: f
    integer,intent(in) :: N_grid

    integer :: j
    real :: pi = 3.141592565

    Int_E = pi*E_hat(1)**0.5*f(1)*dE
    do j = 2,N_grid-1
      Int_E = Int_E+2*pi*E_hat(j)**0.5*f(j)*dE
    enddo
    j = N_grid
    Int_E = Int_E+pi*E_hat(j)**0.5*f(j)*dE

  end function Int_E

  real function Int_E1y(f,N_grid)
    !Integral for E, Int(f*2*pi*E_hat**0.5)dE_hat from 0 to N_grid for energy
    !& Integral y from 0 to 1

    implicit none

    real,dimension(ny,nE),intent(in) :: f
    integer,intent(in) :: N_grid

    integer :: j
    real :: pi = 3.141592565

    Int_E1y = pi*E_hat(1)**0.5*Int_y(f(:,1))*dE
    do j = 2,N_grid-1
      Int_E1y = Int_E1y+2*pi*E_hat(j)**0.5*Int_y(f(:,j))*dE
    enddo
    j = N_grid
    Int_E1y = Int_E1y+pi*E_hat(j)**0.5*Int_y(f(:,j))*dE

  end function Int_E1y

  real function Int_E2y(f,N_grid)
    !Integral for E, Int(f*2*pi*E_hat**1.5)dE_hat from 0 to N_grid for energy
    !& Integral y from 0 to 1

    implicit none

    real,dimension(ny,nE),intent(in) :: f
    integer,intent(in) :: N_grid

    integer :: j
    real :: pi = 3.141592565

    Int_E2y = pi*E_hat(1)**1.5*Int_y(f(:,1))*dE
    do j = 2,N_grid-1
      Int_E2y = Int_E2y+2*pi*E_hat(j)**1.5*Int_y(f(:,j))*dE
    enddo
    j = N_grid
    Int_E2y = Int_E2y+pi*E_hat(j)**1.5*Int_y(f(:,j))*dE

  end function Int_E2y

end module EPtran_use_transport
