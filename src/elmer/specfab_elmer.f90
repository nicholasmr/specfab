include "elmer/include/IBOF.f90"
    
function a4_IBOF(a2) 

    implicit none
    
    real(kind=dp), intent(in) :: a2(3,3)
    real(kind=dp)             :: a4_IBOF(3,3,3,3)
    real(kind=dp)             :: ae2(6), ae4_IBOF(9)
    
    ! IBOF routine assumes dimension ae2(6), but ae2-methods assume ae2(5)
    ae2(1)=a2(1,1)
    ae2(2)=a2(2,2)
    ae2(3)=a2(3,3)
    ae2(4)=a2(1,2)
    ae2(5)=a2(2,3)
    ae2(6)=a2(1,3)
    
    call IBOF(ae2, ae4_IBOF)
    a4_IBOF = ae4_to_a4(a2_to_ae2(a2), ae4_IBOF)
end

function a2_to_ae2(a2) result(ae2)

    implicit none
    
    real(kind=dp), intent(in) :: a2(3,3)
    real(kind=dp)             :: ae2(5)
    
    ae2 = 0.0
    
    ae2(1) = a2(1, 1)
    ae2(2) = a2(2, 2)
    ae2(3) = a2(1, 2)
    ae2(4) = a2(2, 3)
    ae2(5) = a2(1, 3)
end

function a4_to_ae4(a4) result(ae4)

    !* 1111, 2222, 1122, 1123, 2231, 1131, 1112, 2223, 2212
    
    implicit none
    
    real(kind=dp), intent(in) :: a4(3,3,3,3)
    real(kind=dp)             :: ae4(9)
    
    ae4 = 0.0 ! init
    
    ae4(1) = a4(1,1,1,1)
    ae4(2) = a4(2,2,2,2)
    ae4(3) = a4(1,1,2,2)
    ae4(4) = a4(1,1,2,3)
    ae4(5) = a4(1,2,2,3)
    ae4(6) = a4(1,1,1,3)
    ae4(7) = a4(1,1,1,2)
    ae4(8) = a4(2,2,2,3)
    ae4(9) = a4(1,2,2,2)
end


function ae2_to_a2(ae2) result(a2)

    implicit none

    real(kind=dp), intent(in) :: ae2(5)
    real(kind=dp)             :: a2(3,3)
    
    a2 = 0.0
    include "elmer/include/ae2_to_a2__body.f90"
end

function ae4_to_a4(ae2, ae4) result(a4)

    implicit none

    real(kind=dp), intent(in) :: ae2(5), ae4(9)
    Real(kind=dp)             :: a4(3,3,3,3)
    integer :: i,j,k,l

    a4 = 0.0 ! init
    include "elmer/include/ae4_to_a4__body.f90"
end

