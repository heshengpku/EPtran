!------------------------------------------------------------------
! EPtran_driver.f90
!
! PURPOSE:
!  driver calls EPtran_mainsub
!------------------------------------------------------------------

program EPtran_driver

  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------
   
    print *, '----------------------------------------'
    print *, 'EPtran_driver calling EPtran_mainsub'
    print *, '----------------------------------------'
   call EPtran_mainsub 

   call EPtran_write_output

    print *, '----------------------------------------'
    print *, 'EPtran_driver calling Transport'
    print *, '----------------------------------------'
    
   call EPtran_transport

    print *, '----------------------------------------'
    print *, 'EPtran_driver done'
    print *, '----------------------------------------'

end program EPtran_driver
