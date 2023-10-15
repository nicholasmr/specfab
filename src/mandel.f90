! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2022

! Mandel's notation (like the Voigt notation) provides a way to represent symmetric tensors in terms of matrices and vectors.
! In this compact notation (which exclusively accounts for the independent tensor components), inner products can conveniently be expressed as matrix-vector multiplications.

module mandel 

    implicit none 

    integer, parameter, private       :: dp = 8 ! Default precision
    real(kind=dp), parameter, private :: s = sqrt(2.0d0) 
        
contains      

    function mat_to_vec(M) result (v)

        ! 3x3 symmetric matrix to 6x1 Mandel vector

        implicit none
        real(kind=dp), intent(in) :: M(3,3) 
        real(kind=dp)             :: v(6) 
        
        v(:) = [M(1,1), M(2,2), M(3,3), s*M(2,3), s*M(1,3), s*M(1,2)]
    end

    function vec_to_mat(v) result (M)
        
        ! 6x1 Mandel vector to symmetric 3x3 matrix
        
        implicit none
        real(kind=dp), intent(in) :: v(6) 
        real(kind=dp)             :: M(3,3) 
        
        M(1,:) = [v(1),   v(6)/s, v(5)/s]
        M(2,:) = [v(6)/s, v(2),   v(4)/s]
        M(3,:) = [v(5)/s, v(4)/s, v(3)]
    end

    function vec_to_mat_voigt(v) result (M)
        
        ! 6x1 Voigt vector to symmetric 3x3 matrix
        
        implicit none
        real(kind=dp), intent(in) :: v(6) 
        real(kind=dp)             :: M(3,3) 
        
        M(1,:) = [v(1),   v(6)/1, v(5)/1]
        M(2,:) = [v(6)/1, v(2),   v(4)/1]
        M(3,:) = [v(5)/1, v(4)/1, v(3)]
    end

    function a4_to_mat(A) result (M)
    
        ! 3x3x3x3 symmetric tensor (e.g. a4) to 6x6 Mandel matrix

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

    function mat3333_to_mat99(Q) result(P)
        ! 3x3x3x3 arbitrary tensor to 9x9 matrix
        implicit none
        real(kind=dp), intent(in)  :: Q(3,3,3,3) 
        real(kind=dp)              :: P(9,9) 
        include "include/P3333_to_P99.f90"
    end
    
    function mat33_to_vec9(Q) result(V)
        ! 3x3 symmetric tensor to 9x1 vector
        implicit none
        real(kind=dp), intent(in)  :: Q(3,3) 
        real(kind=dp)              :: V(9) 
        V = [Q(1,1),Q(1,2),Q(1,3), & 
             Q(1,2),Q(2,2),Q(2,3), &
             Q(1,3),Q(2,3),Q(3,3)]
    end

    function vec9_to_mat33(V) result(Q)
        ! 9x1 vector to 3x3 tensor
        implicit none
        real(kind=dp), intent(in)  :: V(9) 
        real(kind=dp)              :: Q(3,3) 
        Q(1,:) = [V(1),V(2),V(3)]
        Q(2,:) = [V(4),V(5),V(6)]
        Q(3,:) = [V(7),V(8),V(9)]
    end

end module mandel
