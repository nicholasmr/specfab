! N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

! *** IN DEVELOPMENT ***

module wavepropagation  

    use tensorproducts
    use mandel
    use moments

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    integer, private            :: nn ! loop index
!    integer, parameter :: identity(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    
       
contains      

    function vj(nlm, theta_n,phi_n, omega, rho, lam,mu,eta1,eta2,eta3) 
    
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: theta_n(:), phi_n(:) ! arrays of theta and phi values to calculate phase velocities (vj) along
        complex(kind=dp)             :: nlm4(15), nlm4_rot(15) ! nlm truncated at L=4, and its rotated form
        real(kind=dp), intent(in)    :: omega, rho, lam,mu,eta1,eta2,eta3
        real(kind=dp)                :: ri(0:6)
        complex(kind=dp)             :: qjn(3) ! z-comp of wavector (k=[0,0,q]) for qS1, qS2, qP waves
        real(kind=dp)                :: vj(3,size(theta_n)) ! qS1, qS2, qP phase velocities
        
        nlm4 = nlm(:(I_l6-1)) ! l=0,2,4 coef array
        
        do nn = 1,size(theta_n)
            nlm4_rot = rotate_nlm4(nlm4, -theta_n(nn), -phi_n(nn)) ! negative angles because we are are rotating the specified direction (back) into the vertical orientation
            ri = detQ_coefs(nlm, 0.0d0, omega, rho, lam,mu,eta1,eta2,eta3)
            qjn = qj(ri) 
            vj(:,nn) = omega/abs(qjn(:))
        end do
    end

    function rotate_nlm4(nlm4, theta,phi) result (nlm4_rot)
    
        implicit none

        complex(kind=dp), intent(in) :: nlm4(15) ! nlm truncated at L=4
        real(kind=dp), intent(in)    :: theta, phi 
        complex(kind=dp)             :: nlm4_rot(15) !, n2m_rot(5), n4m_rot(9)
        complex(kind=dp)             :: D2mn(-2:2,-2:2), D4mn(-4:4,-4:4)
    
        include "include/Dlmn.f90"
        nlm4_rot(I_l0) = nlm4(I_l0)
        nlm4_rot(I_l2:(I_l4-1)) = matmul(D2mn, nlm4(I_l2:(I_l4-1)))
        nlm4_rot(I_l4:(I_l6-1)) = matmul(D4mn, nlm4(I_l4:(I_l6-1)))
    end

    function qj(ri) 

        implicit none
        
!        complex(kind=dp), intent(in) :: nlm(:)
!        real(kind=dp), intent(in)    :: omega, rho, lam,mu,eta1,eta2,eta3
        real(kind=dp), intent(in)    :: ri(0:6)
        complex(kind=dp)             :: qj(3) ! z-comp of wavector (k=[0,0,q]) for qS1, qS2, qP waves
        complex(kind=dp)             :: a,b, x1,x2,x3
        
        include "include/detQ_xi0_sol.f90"
    end

    function detQ_coefs(nlm, xi, omega, rho, lam,mu,eta1,eta2,eta3) result(ri)

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: xi ! Wave vector x-component: k = [xi,0,q] (q solved for here)
        real(kind=dp), intent(in)    :: omega, rho, lam,mu,eta1,eta2,eta3
        real(kind=dp)                :: ri(0:6)
        real(kind=dp)                :: B(6,6), Mv(6), M(3,3), w
        real(kind=dp)                :: a11,a22,a33,a23,a13,a12
        real(kind=dp)                :: b11,b22,b33,b23,b13,b12
        real(kind=dp)                :: c11,c22,c33,c23,c13,c12
        
        call f_ev_ck_Mandel(nlm, Mv, B) ! Structure tensors in Mandel notation: Mv = ev_c2_Mandel, B = ev_c4_Mandel
        M = vec_to_mat(Mv) ! = a2
        w = omega**2 * rho
        include "include/detQ_coefs.f90"
    end


end module wavepropagation 
