! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

module moments 

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real, parameter, private    :: Pi = 3.1415927
!    integer, private :: ii

contains      

    ! Moments (structure tensors) <c^k> for k = 0,2,4,6,8

    function f_ev_c0(n00) result(ev)
        ! Integral over orientation distribution = total number of c-axes (grains)
        complex(kind=dp), intent(in) :: n00
        real(kind=dp) :: ev
        ev = REAL(sqrt(4*Pi)*n00)
    end

    function f_ev_c2(n00,n2m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2)
        real(kind=dp) :: k, ev(3,3)
        include "include/ev_c2__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c4(n00,n2m,n4m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4)
        real(kind=dp) :: k, ev(3,3, 3,3)
        include "include/ev_c4__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c6(n00,n2m,n4m,n6m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        real(kind=dp) :: k, ev(3,3, 3,3, 3,3)
!        k = 0
!        ev = reshape([(0, ii=1,3**6)], [3,3, 3,3, 3,3])
        include "include/ev_c6__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c8(n00,n2m,n4m,n6m,n8m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6), n8m(-8:8)
        real(kind=dp) :: k, ev(3,3, 3,3, 3,3, 3,3)
!        k = 0
!        ev = reshape([(0, ii=1,3**8)], [3,3, 3,3, 3,3, 3,3])
        include "include/ev_c8__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end

end module moments
