! N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

module mandel 

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real(kind=dp), parameter, private    :: s = sqrt(2.0d0) 
        
contains      

    function mat_to_vec(M) result (v)

        implicit none

        real(kind=dp), intent(in) :: M(3,3) 
        real(kind=dp)             :: v(6) 
        
        v(:) = [M(1,1), M(2,2), M(3,3), s*M(2,3), s*M(1,3), s*M(1,2)]
    end

    function vec_to_mat(v) result (M)
        
        implicit none
        real(kind=dp), intent(in) :: v(6) 
        real(kind=dp)             :: M(3,3) 
        
        M(1,:) = [v(1),   v(6)/s, v(5)/s]
        M(2,:) = [v(6)/s, v(2),   v(4)/s]
        M(3,:) = [v(5)/s, v(4)/s, v(3)]
    end

    function a4_to_mat(A) result (M)

        implicit none
        real(kind=dp), intent(in) :: A(3,3,3,3) 
        real(kind=dp)             :: M(6,6)

        M(1,:) = [  A(1,1,1,1),   A(1,1,2,2),   A(1,1,3,3), s*A(1,1,2,3), s*A(1,1,1,3), s*A(1,1,1,2)]
        M(2,:) = [  A(2,2,1,1),   A(2,2,2,2),   A(2,2,3,3), s*A(2,2,2,3), s*A(2,2,1,3), s*A(2,2,1,2)]
        M(3,:) = [  A(3,3,1,1),   A(3,3,2,2),   A(3,3,3,3), s*A(3,3,2,3), s*A(3,3,1,3), s*A(3,3,1,2)]
        M(4,:) = [s*A(2,3,1,1), s*A(2,3,2,2), s*A(2,3,3,3), 2*A(2,3,2,3), 2*A(2,3,1,3), 2*A(2,3,1,2)]
        M(5,:) = [s*A(1,3,1,1), s*A(1,3,2,2), s*A(1,3,3,3), 2*A(1,3,2,3), 2*A(1,3,1,3), 2*A(1,3,1,2)]
        M(6,:) = [s*A(1,2,1,1), s*A(1,2,2,2), s*A(1,2,3,3), 2*A(1,2,2,3), 2*A(1,2,1,3), 2*A(1,2,1,2)]
    end

end module mandel
