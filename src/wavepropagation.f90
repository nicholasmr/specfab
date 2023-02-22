! Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2023

module wavepropagation  

    use tensorproducts
    use mandel
    use elasticities
    use homogenizations
    use moments
    use dynamics
    use rotation

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    integer, private            :: ii,nn ! loop index
    
!    real(kind=dp), private, parameter :: k(3) = [0.0d0, 0.0d0, 1.0d0]
    real(kind=dp), private, parameter :: kk(3,3) = reshape([0.0d0,0.0d0,0.0d0, 0.0d0,0.0d0,0.0d0, 0.0d0,0.0d0,1.0d0], [3,3]) ! k \otimes k, where k = [0,0,1]
       
contains      

    !---------------------------------
    ! ELASTIC
    !---------------------------------

    !---------------------------------
    ! Orthotropic grains
    !---------------------------------

    !--------------
    ! Parameter definitions
    !--------------
    ! alpha:  Reuss--Voigt homogenization weight
    ! lam(6): Monocrystal parameters (lam11,lam22,lam33,lam12,lam13,lam23)
    ! mu(3):  Monocrystal parameters (mu1,mu2,mu3)  
    !--------------
    ! r_i are crystallographic axes in paper, but here they are *defined* as
    !--------------
    ! r_1 := b
    ! r_2 := n
    ! r_3 := v
    !--------------

    function ai_to_Qnorm_orthotropic(a2_r1,a2_r2,a2_r3, a4v_r1,a4v_r2,a4v_r3, a4_r1r2,a4_r1r3,a4_r2r3, lam,mu) result (Q)
       
        implicit none
        
        real(kind=dp), intent(in) :: a2_r1(3,3),    a2_r2(3,3),    a2_r3(3,3)
        real(kind=dp), intent(in) :: a4v_r1(6,6),   a4v_r2(6,6),   a4v_r3(6,6)
        real(kind=dp), intent(in) :: a4_r1r2(3,3,3,3), a4_r1r3(3,3,3,3), a4_r2r3(3,3,3,3)
        real(kind=dp), intent(in) :: lam(6), mu(3) 
        real(kind=dp)             :: Q(3,3)
    
        Q = 0.0d0 ! initialize
        Q = Q + mu(1)/2*(a2_r1 + matmul(kk,a2_r1) + matmul(a2_r1,kk) + doubleinner22(a2_r1,kk)*identity)
        Q = Q + mu(2)/2*(a2_r2 + matmul(kk,a2_r2) + matmul(a2_r2,kk) + doubleinner22(a2_r2,kk)*identity)
        Q = Q + mu(3)/2*(a2_r3 + matmul(kk,a2_r3) + matmul(a2_r3,kk) + doubleinner22(a2_r3,kk)*identity)
        Q = Q + lam(1)*vec_to_mat(matmul(a4v_r1,mat_to_vec(kk)))
        Q = Q + lam(2)*vec_to_mat(matmul(a4v_r2,mat_to_vec(kk)))
        Q = Q + lam(3)*vec_to_mat(matmul(a4v_r3,mat_to_vec(kk)))
        Q = Q + 2*lam(4)*( doubleinner42_firstlast_symmetric(a4_r1r2,kk) )
        Q = Q + 2*lam(5)*( doubleinner42_firstlast_symmetric(a4_r1r3,kk) )
        Q = Q + 2*lam(6)*( doubleinner42_firstlast_symmetric(a4_r2r3,kk) )
    end
    
    !--------------
    ! For CPO defined by distributions r1,r2,r3
    !--------------
    
    function Vi_elastic_orthotropic(nlm_r1,nlm_r2,nlm_r3, alpha,lam,mu,rho, theta_n,phi_n) result (vi)
    
        ! Elastic phase velocities (qP, qS1, qS2) given orientation distributions nlm_r1, nlm_r2, nlm_r3 (of crystallographic axes r1,r2,r3) and propagation directions (theta_n,phi_n) 
        
        ! *** nlm_r3 IS NOT USED BUT ESTIMATED FROM nlm_1 AND nlm_2 ***
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
        real(kind=dp), intent(in)    :: alpha, lam(6), mu(3), rho
        real(kind=dp), intent(in)    :: theta_n(:), phi_n(:) ! arrays of theta and phi values to calculate phase velocities (vj) along
        real(kind=dp)                :: Qnorm(3,3)
        real(kind=dp)                :: vi(3,size(theta_n)) ! qS1, qS2, qP phase velocities (in that order)
        
        complex(kind=dp)             :: nlm4_r1(15), nlm4_rot_r1(15) ! truncated at L=4, and its rotated form
        complex(kind=dp)             :: nlm4_r2(15), nlm4_rot_r2(15) 
        complex(kind=dp)             :: nlm4_r3(15), nlm4_rot_r3(15) 
        
        nlm4_r1 = nlm_r1(:(I_l6-1)) ! l=0,2,4 coefficients
        nlm4_r2 = nlm_r2(:(I_l6-1)) 
        nlm4_r3 = nlm_r3(:(I_l6-1)) 

        do nn = 1,size(theta_n)
            nlm4_rot_r1 = rotate_nlm4(nlm4_r1, -theta_n(nn), -phi_n(nn)) ! negative angles because we are are rotating the specified direction (back) into the vertical orientation
            nlm4_rot_r2 = rotate_nlm4(nlm4_r2, -theta_n(nn), -phi_n(nn)) 
            nlm4_rot_r3 = rotate_nlm4(nlm4_r3, -theta_n(nn), -phi_n(nn)) 
            Qnorm = Qnorm_orthotropic(nlm4_rot_r1,nlm4_rot_r2,nlm4_rot_r3, alpha,lam,mu) ! effective acoustic tensor
            vi(:,nn) = Qnorm_to_vi(Qnorm, rho)
        end do
    end
    
    function Qnorm_orthotropic(nlm_r1,nlm_r2,nlm_r3, alpha,lam,mu) result (Qnorm)

        ! Normalized Voigt--Reuss-averaged acoustic tensor, Qnorm, such that Q = k^2 * Qnorm
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
        real(kind=dp), intent(in)    :: alpha, lam(6), mu(3)
        real(kind=dp)                :: Qnorm(3,3), Qnorm_Reuss(3,3), Qnorm_Voigt(3,3)

!        if (alpha .lt. 1.0d0) Qn_Reuss = Qnorm_orthotropic_Reuss(nlm, lam,mu, Elam,Emu,Egam)
        Qnorm_Reuss = 0.0d0
        if (alpha .gt. 0.0d0) Qnorm_Voigt = Qnorm_orthotropic_Voigt(nlm_r1,nlm_r2,nlm_r3, lam,mu)
        Qnorm = (1-alpha)*Qnorm_Reuss + alpha*Qnorm_Voigt
    end
    
    function Qnorm_orthotropic_Voigt(nlm_r1,nlm_r2,nlm_r3, lam,mu) result (Q)

        ! Acoustic tensor, Q=Q_Voigt, such that div(stress) = matmul(Q,u) where u is a plane wave in the strain field
        ! Assumes k = [0,0,kz]
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
        real(kind=dp), intent(in)    :: lam(6), mu(3) 
        real(kind=dp)                :: Q(3,3)
        
        real(kind=dp) :: a2_r1(3,3), a2v_r1(6), a4v_r1(6,6)
        real(kind=dp) :: a2_r2(3,3), a2v_r2(6), a4v_r2(6,6)
        real(kind=dp) :: a2_r3(3,3), a2v_r3(6), a4v_r3(6,6) 

        real(kind=dp) :: a4_r1r2(3,3,3,3), a4_r1r3(3,3,3,3), a4_r2r3(3,3,3,3)

        call f_ev_ck_Mandel(nlm_r1, a2v_r1, a4v_r1) ! Structure tensors in Mandel notation: a2v = ev_c2_Mandel, a4v = ev_c4_Mandel
        call f_ev_ck_Mandel(nlm_r2, a2v_r2, a4v_r2) 
        call f_ev_ck_Mandel(nlm_r3, a2v_r3, a4v_r3) 
        a2_r1 = vec_to_mat(a2v_r1) 
        a2_r2 = vec_to_mat(a2v_r2) 
        a2_r3 = vec_to_mat(a2v_r3) 

        ! Preferred method for cross structure tensors
        a4_r1r2 = a4_joint(nlm_r1,nlm_r2)
        a4_r1r3 = a4_jointcross(nlm_r1,nlm_r2) ! r3 := r1 x r2
        a4_r2r3 = a4_jointcross(nlm_r2,nlm_r1) ! r3 := r1 x r2

        ! Debug
        ! This should reproduce the discrete method's result for ODFs that are delta funcs (=single grain behavior)
!        a4_r1r2 = outerprod22(a2_r1,a2_r2)
!        a4_r1r3 = outerprod22(a2_r1,a2_r3)
!        a4_r2r3 = outerprod22(a2_r2,a2_r3)

        ! Debug
!        a4_r1r2 = a4_joint(nlm_r1,nlm_r2)
!        a4_r1r3 = a4_joint(nlm_r1,nlm_r3)
!        a4_r2r3 = a4_joint(nlm_r2,nlm_r3)

        Q = ai_to_Qnorm_orthotropic(a2_r1,a2_r2,a2_r3, a4v_r1,a4v_r2,a4v_r3, a4_r1r2,a4_r1r3,a4_r2r3, lam,mu) 
    end
    
    !--------------
    ! For CPO defined by discrete ensemble of axes r1,r2,r3
    !--------------
    
    function Vi_elastic_orthotropic__discrete(r1,r2,r3, alpha,lam,mu,rho, theta_n,phi_n) result (vi)
    
        ! Elastic phase velocities (qP, qS1, qS2) given ensemble of crystallographic axes (r1,r2,r3) and propagation directions (theta_n,phi_n) 
        
        implicit none
        
        real(kind=dp), intent(in) :: r1(:,:), r2(:,:), r3(:,:) ! (3,N)
        real(kind=dp), intent(in) :: alpha, lam(6), mu(3), rho
        real(kind=dp), intent(in) :: theta_n(:), phi_n(:) ! arrays of theta and phi values to calculate phase velocities (vj) along
        real(kind=dp)             :: Qnorm(3,3)
        real(kind=dp)             :: vi(3,size(theta_n)) ! qS1, qS2, qP phase velocities (in that order)

        real(kind=dp), dimension(3,size(r1,dim=2)) :: r1_rot, r2_rot, r3_rot
        
        do nn = 1,size(theta_n)
            do ii = 1,size(r1,dim=2)
                r1_rot(:,ii) = rotate_vector(r1(:,ii), -theta_n(nn), -phi_n(nn)) ! negative angles because we are are rotating the specified direction (back) into the vertical orientation
                r2_rot(:,ii) = rotate_vector(r2(:,ii), -theta_n(nn), -phi_n(nn)) 
                r3_rot(:,ii) = rotate_vector(r3(:,ii), -theta_n(nn), -phi_n(nn))
            end do
            Qnorm = Qnorm_orthotropic_Voigt__discrete(r1_rot,r2_rot,r3_rot, lam,mu) ! effective acoustic tensor
            vi(:,nn) = Qnorm_to_vi(Qnorm, rho)
        end do
    end
    
    function Qnorm_orthotropic_Voigt__discrete(r1,r2,r3, lam,mu) result (Q)

        ! Acoustic tensor, Q=Q_Voigt, such that div(stress) = matmul(Q,u) where u is a plane wave in the strain field
        ! Assumes k = [0,0,kz]
        
        implicit none
        
        real(kind=dp), intent(in) :: r1(:,:), r2(:,:), r3(:,:) ! (3,N)
        real(kind=dp), intent(in) :: lam(6), mu(3)
        real(kind=dp)             :: Q(3,3)
        integer                   :: N

        real(kind=dp) :: a2_r1_(3,3),      a2_r2_(3,3),      a2_r3_(3,3)
        real(kind=dp) :: a2_r1(3,3),       a2_r2(3,3),       a2_r3(3,3)
        real(kind=dp) :: a4_r1(3,3,3,3),   a4_r2(3,3,3,3),   a4_r3(3,3,3,3)
        real(kind=dp) :: a4_r1r2(3,3,3,3), a4_r1r3(3,3,3,3), a4_r2r3(3,3,3,3)

        a2_r1 = 0.0d0
        a2_r2 = 0.0d0
        a2_r3 = 0.0d0
        
        a4_r1 = 0.0d0
        a4_r2 = 0.0d0
        a4_r3 = 0.0d0
        
        a4_r1r2 = 0.0d0
        a4_r1r3 = 0.0d0
        a4_r2r3 = 0.0d0

        N = size(r1, dim=2) 
        do ii=1,N
            a2_r1_ = outerprod(r1(:,ii),r1(:,ii))
            a2_r2_ = outerprod(r2(:,ii),r2(:,ii))
            a2_r3_ = outerprod(r3(:,ii),r3(:,ii))
            
            a2_r1 = a2_r1 + a2_r1_/N
            a2_r2 = a2_r2 + a2_r2_/N
            a2_r3 = a2_r3 + a2_r3_/N
            
            a4_r1 = a4_r1 + outerprod22(a2_r1_, a2_r1_)/N
            a4_r2 = a4_r2 + outerprod22(a2_r2_, a2_r2_)/N
            a4_r3 = a4_r3 + outerprod22(a2_r3_, a2_r3_)/N
            
            a4_r1r2 = a4_r1r2 + outerprod22(a2_r1_, a2_r2_)/N
            a4_r1r3 = a4_r1r3 + outerprod22(a2_r1_, a2_r3_)/N
            a4_r2r3 = a4_r2r3 + outerprod22(a2_r2_, a2_r3_)/N
        end do

        Q = ai_to_Qnorm_orthotropic(a2_r1,a2_r2,a2_r3, a4_to_mat(a4_r1),a4_to_mat(a4_r2),a4_to_mat(a4_r3), a4_r1r2,a4_r1r3,a4_r2r3, lam,mu) 
    end
    
    !---------------------------------
    ! Transversely isotropic grains
    !---------------------------------

    function Vi_elastic_tranisotropic(nlm, alpha, lam,mu,Elam,Emu,Egam, rho, theta_n,phi_n) result (vj)
    
        ! Elastic phase velocities (qP, qS1, qS2) given nlm (i.e. ODF) and propagation directions (theta_n,phi_n) 
        ! Assumes transversely isotropic grains with parameters lam,mu,Elam,Emu,Egam
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: alpha, lam,mu,Elam,Emu,Egam, rho
        real(kind=dp), intent(in)    :: theta_n(:), phi_n(:) ! arrays of theta and phi values to calculate phase velocities (vj) along

        real(kind=dp)                :: vj(3,size(theta_n)) ! qS1, qS2, qP phase velocities
        
        complex(kind=dp)             :: nlm4(15), nlm4_rot(15) ! nlm truncated at L=4, and its rotated form
        real(kind=dp)                :: Qn(3,3), ri(0:6) ! Qn = Qnorm
        complex(kind=dp)             :: a,b,denom,D, x1,x2,x3 ! vars for polynomial solver
        complex(kind=dp)             :: qj(3) ! z-comp of wavevector, k=[0,0,q], for qS1, qS2, qP waves
        real(kind=dp)                :: omega ! Wave angular velocity 
        
        nlm4 = nlm(:(I_l6-1)) ! l=0,2,4 coefficients
        
        omega = 1.0d0 ! Eigenvalues of problem are rho*V^2 where V=omega/k, so the problem does not, in fact, depend on omega. It is included here anyway for readability (and because of the way the analytical solutions were solved for in Mathematica).
        
        do nn = 1,size(theta_n)
            nlm4_rot = rotate_nlm4(nlm4, -theta_n(nn), -phi_n(nn)) ! negative angles because we are are rotating the specified direction (back) into the vertical orientation
            Qn = Qnorm(nlm4_rot, alpha, lam,mu,Elam,Emu,Egam) ! effective acoustic tensor
            ri = detQ_polycoefs(Qn, omega, rho) ! polynomial coefs of Christoffel problem

            ! xi=0 solutions from depressed cubic equation r0 + r2*q^2 + r4*q^4 + r6*q^6 = 0 ==> x^3 + a*x + b = 0
            ! Exported from Mathematica
            a = (3.0d0*ri(6)*ri(2) - ri(4)**2)/(3.0d0*ri(6)**2)
            b = (2.0d0*ri(4)**3 - 9.0d0*ri(6)*ri(4)*ri(2) + 27.0d0*ri(6)**2*ri(0))/(27.0d0*ri(6)**3)
            denom = (6.0d0,0.0d0)**0.6666666666666666*((-9.0d0,0.0d0)*b + Sqrt((12.0d0,0.0d0)*a**(3.0d0,0.0d0) + (81.0d0,0.0d0)*b**(2.0d0,0.0d0)))**0.3333333333333333
            D = ((-9.0d0,0.0d0)*b + Sqrt((12.0d0,0.0d0)*a**(3.0d0,0.0d0) + (81.0d0,0.0d0)*b**(2.0d0,0.0d0)))**0.6666666666666666
            x1 = -ri(4)/(3.0d0*ri(6)) + ((-2.0d0,0.0d0)*(3.0d0,0.0d0)**0.3333333333333333*a + (2.0d0,0.0d0)**0.3333333333333333*D)/denom
            x2 = -ri(4)/(3.0d0*ri(6)) + ((-1.0d0,0.0d0)**0.3333333333333333*((2.0d0,0.0d0)*(3.0d0,0.0d0)**0.3333333333333333*a + (-2.0d0,0.0d0)**0.3333333333333333*D))/denom
            x3 = -ri(4)/(3.0d0*ri(6)) + ((-1.0d0,0.0d0)**1.3333333333333333*((2.0d0,0.0d0)*(-3.0d0,0.0d0)**0.3333333333333333*a + (2.0d0,0.0d0)**0.3333333333333333*D))/denom
            qj(1) = sqrt(x3)
            qj(2) = sqrt(x1)
            qj(3) = sqrt(x2)
            
            ! phase velocities are then by definition
            vj(:,nn) = omega/abs(qj(:))
        end do
    end

    function Qnorm(nlm, alpha, lam,mu,Elam,Emu,Egam) result (Qn)

        ! Normalized Voigt--Reuss-averaged acoustic tensor, Qnorm, such that Q = k^2 * Qnorm
        ! alpha is the linear Voigt--Reuss weight (defined below)
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: alpha, lam,mu,Elam,Emu,Egam ! Reuss--Voigt weight (alpha) and monocrystal parameters
        real(kind=dp)                :: Qn(3,3), Qn_Reuss(3,3), Qn_Voigt(3,3)

        Qn_Reuss = 0.0d0
        Qn_Voigt = 0.0d0

        if (alpha .lt. 1.0d0) then
            Qn_Reuss = Qnorm_Reuss(nlm, lam,mu, Elam,Emu,Egam)
        end if

        if (alpha .gt. 0.0d0) then
            Qn_Voigt = Qnorm_Voigt(nlm, lam,mu, Elam,Emu,Egam)
        end if 
        
        Qn = (1-alpha)*Qn_Reuss + alpha*Qn_Voigt
    end
    
    function Qnorm_Voigt(nlm, lam,mu,Elam,Emu,Egam) result (Q)

        ! Acoustic tensor, Q=Q_Voigt, such that div(stress) = matmul(Q,u) where u is a plane wave in the strain field
        ! Assumes k = [0,0,kz]
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: lam,mu,Elam,Emu,Egam ! Reuss--Voigt weight (alpha) and monocrystal parameters
        real(kind=dp)                :: Q(3,3)

        real(kind=dp) :: k1,k2,k3,k4,k5
        real(kind=dp) :: a2mat(3,3), a2v(6), a4v(6,6)
        real(kind=dp) :: KM_anticomm(3,3)

        call elas_revparams_tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5) ! (k1, ..., k5) are effective elastic coefficients, *not* related to wave vector.
        call f_ev_ck_Mandel(nlm, a2v, a4v) ! Structure tensors in Mandel notation: a2v = ev_c2_Mandel, B = ev_c4_Mandel
        a2mat = vec_to_mat(a2v) ! = a2
        KM_anticomm = matmul(kk,a2mat) + matmul(a2mat,kk) ! = {a2,kk} for k=kvert=[0,0,1]
        Q = k2/2*identity + (k1+k2/2)*kk + k3*KM_anticomm + k4*vec_to_mat(matmul(a4v,mat_to_vec(kk))) + 0.5d0*k5*(a2mat + doubleinner22(a2mat,kk)*identity + KM_anticomm)
    end
    
    function Qnorm_Reuss(nlm, lam,mu, Elam,Emu,Egam) result(Q)

        ! Acoustic tensor, Q=Q_Reuss, such that div(stress) = matmul(Q,u) where u is a plane wave in the strain field
        ! Assumes k = [0,0,kz]
    
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: lam,mu, Elam,Emu,Egam ! monocrystal parameters
        real(kind=dp)                :: Q(3,3) 

        real(kind=dp)    :: k1,k2,k3,k4,k5 ! aux params (not related to wave vector)
        complex(kind=dp) :: n00, n2m(-2:2), n4m(-4:4)
        real(kind=dp)    :: ev_c2(3,3), ev_c4(3,3,3,3) 
        real(kind=dp)    :: P(9,9), Pinv(9,9) 

        call elas_fwdparams_tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5) ! (k1, ..., k5) are effective elastic coefficients, *not* related to wave vector.

        n00 = nlm(1)
        n2m = nlm(I_l2:(I_l4-1))
        n4m = nlm(I_l4:(I_l6-1))
        ev_c2 = f_ev_c2(n00,n2m)
        ev_c4 = f_ev_c4(n00,n2m,n4m)

        ! Set P such that strain = matmul(P,stress)
        include "include/P_Reuss.f90"        
          
        Pinv = matinv(P)
        
        Q(1,:) = 1/2.0d0 * [Pinv(3,3)+Pinv(3,7), Pinv(3,6)+Pinv(3,8), 2*Pinv(3,9)]
        Q(2,:) = 1/2.0d0 * [Pinv(6,3)+Pinv(6,7), Pinv(6,6)+Pinv(6,8), 2*Pinv(6,9)]
        Q(3,:) = 1/2.0d0 * [Pinv(9,3)+Pinv(9,7), Pinv(9,6)+Pinv(9,8), 2*Pinv(9,9)]
    end
    
    !---------------------------------
    ! ELECTROMAGNETIC
    !---------------------------------

    ! N/A

    !---------------------------------
    ! SUPPORTING ROUTINES
    !---------------------------------
    
    function Qnorm_to_vi(Qnorm, rho) result(vi)
    
        ! Solve eigenvalue problem for qP, qS1, and qS2 phase velocities 

        implicit none

        real(kind=dp), intent(in)    :: Qnorm(3,3), rho ! Acoustic tensor, mass density
        complex(kind=dp)             :: qi(3) ! z-comp of wavevector (k=[0,0,q]) for qS1, qS2, qP waves
        real(kind=dp)                :: omega, ri(0:6), vi(3) ! vi = (qS1, qS2, qP) phase velocities
        complex(kind=dp)             :: a,b,denom,D, x1,x2,x3 ! vars for polynomial solver

        ! Eigenvalues of problem are rho*v^2 where v=omega/k, so the problem does not, in fact, depend on omega. 
        ! It is included here anyway for readability (and because of the way the analytical solutions were solved for in Mathematica).
        omega = 1.0d0 

        ri = detQ_polycoefs(Qnorm, omega, rho) ! Polynomial coefficients of Christoffel problem

        ! xi=0 solutions from depressed cubic equation r0 + r2*q^2 + r4*q^4 + r6*q^6 = 0 ==> x^3 + a*x + b = 0
        ! --- Exported from Mathematica ---
        a = (3.0d0*ri(6)*ri(2) - ri(4)**2)/(3.0d0*ri(6)**2)
        b = (2.0d0*ri(4)**3 - 9.0d0*ri(6)*ri(4)*ri(2) + 27.0d0*ri(6)**2*ri(0))/(27.0d0*ri(6)**3)
        denom = (6.0d0,0.0d0)**0.6666666666666666*((-9.0d0,0.0d0)*b + Sqrt((12.0d0,0.0d0)*a**(3.0d0,0.0d0) + (81.0d0,0.0d0)*b**(2.0d0,0.0d0)))**0.3333333333333333
        D = ((-9.0d0,0.0d0)*b + Sqrt((12.0d0,0.0d0)*a**(3.0d0,0.0d0) + (81.0d0,0.0d0)*b**(2.0d0,0.0d0)))**0.6666666666666666
        x1 = -ri(4)/(3.0d0*ri(6)) + ((-2.0d0,0.0d0)*(3.0d0,0.0d0)**0.3333333333333333*a + (2.0d0,0.0d0)**0.3333333333333333*D)/denom
        x2 = -ri(4)/(3.0d0*ri(6)) + ((-1.0d0,0.0d0)**0.3333333333333333*((2.0d0,0.0d0)*(3.0d0,0.0d0)**0.3333333333333333*a + (-2.0d0,0.0d0)**0.3333333333333333*D))/denom
        x3 = -ri(4)/(3.0d0*ri(6)) + ((-1.0d0,0.0d0)**1.3333333333333333*((2.0d0,0.0d0)*(-3.0d0,0.0d0)**0.3333333333333333*a + (2.0d0,0.0d0)**0.3333333333333333*D))/denom
        qi(1) = sqrt(x3)
        qi(2) = sqrt(x1)
        qi(3) = sqrt(x2)
        ! --- End ---
        
        ! Phase velocities are then by definition
        vi(:) = omega/abs(qi(:))
    end

    function detQ_polycoefs(Qn, omega, rho) result(ri)

        implicit none
        
        real(kind=dp), intent(in) :: Qn(3,3), omega, rho ! Qn = Qnorm
        real(kind=dp)             :: ri(0:6) ! polynomial coefs
        real(kind=dp)             :: w, I1,I2,I3
        
        w = omega**2 * rho
        
        ! Polynomial coefs depend on the 3D tensor invariants
        I1 = Qn(1,1) + Qn(2,2) + Qn(3,3)
        I2 = (-1.0d0)*Qn(1,2)**(2.0d0) + (-1.0d0)*Qn(1,3)**(2.0d0) + Qn(1,1)*Qn(2,2) + (-1.0d0)*Qn(2,3)**(2.0d0) + Qn(1,1)*Qn(3,3) + Qn(2,2)*Qn(3,3)
        I3 = (-1.0d0)*Qn(1,3)**(2.0d0)*Qn(2,2) + (2.0d0)*Qn(1,2)*Qn(1,3)*Qn(2,3) + (-1.0d0)*Qn(1,1)*Qn(2,3)**(2.0d0) + (-1.0d0)*Qn(1,2)**(2.0d0)*Qn(3,3) + Qn(1,1)*Qn(2,2)*Qn(3,3)

        ri(0) = -w**3
        ri(1) = 0
        ri(2) = I1*w**2
        ri(3) = 0
        ri(4) = -(I2*w)
        ri(5) = 0
        ri(6) = I3
    end

    function matinv(A) result(Ainv)

      ! Source: https://fortranwiki.org/fortran/show/Matrix+inversion

      real(dp), dimension(:,:), intent(in)     :: A
      real(dp), dimension(size(A,1),size(A,2)) :: Ainv

      real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1))  :: ipiv   ! pivot indices
      integer :: n, info

!      ! External procedures defined in LAPACK
!      external DGETRF
!      external DGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
         stop 'matinv(A) in wavepropagation.f90: Matrix is numerically singular!'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
         stop 'matinv(A) in wavepropagation.f90: Matrix inversion failed!'
      end if
      
    end 

end module wavepropagation 
