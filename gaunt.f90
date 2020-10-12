! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

module gaunt

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    include "include/gaunt__head_L60.f90"        
        
contains

subroutine set_gaunts()
    include "include/gaunt__body_L60.f90"
end
    
end module gaunt
