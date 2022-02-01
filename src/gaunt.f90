! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

! To allow for larger "L" truncation, update the below gaunt coef matrices by running: python3 make_gaunt_coefs.py L 

module gaunt

    implicit none 
    integer, parameter, private :: dp = 8 ! Default precision
    include "include/gaunt__head.f90"        

contains

    subroutine set_gaunts()
        include "include/gaunt__body.f90"
    end
    
end module gaunt
