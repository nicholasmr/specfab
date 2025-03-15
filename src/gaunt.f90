! Nicholas Rathmann <rathmann@nbi.ku.dk> 2020-

! Update Gaunt coefs by running python3 make_gaunt_coefs.py L 

module gaunt

    implicit none 
    integer, parameter, private :: dp = 8 ! Default precision
    include "include/gaunt__head.f90"
    include "include/gaunt__head_new.f90"

contains

    subroutine set_gaunts()
        include "include/gaunt__body.f90"
        include "include/gaunt__body_new.f90"
    end
    
end module gaunt
