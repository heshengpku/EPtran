subroutine EPtran_mainsub

  !---------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------
    
    call EPtran_read_input

    call EPtran_comp_eq_plasma

    call EPtran_comp_alpha_slowing

    call get_tglf_parameters
    
    print *, '----------------------------------------'
    print *,'EPtran_mainsub done'
    print *, '----------------------------------------'
    
end subroutine EPtran_mainsub
    

