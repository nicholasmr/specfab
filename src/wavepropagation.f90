! Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2022

module wavepropagation  

    use tensorproducts
    use mandel
    use elasticities
    use homogenizations
    use moments
    use dynamics

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    integer, private            :: nn ! loop index
    
    real(kind=dp), private, parameter :: kk(3,3) = reshape([0.0d0,0.0d0,0.0d0, 0.0d0,0.0d0,0.0d0, 0.0d0,0.0d0,1.0d0], [3,3]) ! k \otimes k, where k = [0,0,1]
       
contains      

    !---------------------------------
    ! ELASTIC
    !---------------------------------

    function Vi_elastic_tranisotropic(nlm, alpha, lam,mu,Elam,Emu,Egam, rho, theta_n,phi_n) result (vj)
    
        ! Elastic phase velocities (qP, qS1, qS2) given ODF in terms of nlm and propagation directions (theta_n,phi_n) 
    
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: alpha, lam,mu,Elam,Emu,Egam, rho
        real(kind=dp), intent(in)    :: theta_n(:), phi_n(:) ! arrays of theta and phi values to calculate phase velocities (vj) along

        real(kind=dp)                :: vj(3,size(theta_n)) ! qS1, qS2, qP phase velocities
        
        complex(kind=dp)             :: nlm4(15), nlm4_rot(15) ! nlm truncated at L=4, and its rotated form
        real(kind=dp)                :: Qn(3,3), ri(0:6) ! Qn = Qnorm
        complex(kind=dp)             :: a,b,denom,D, x1,x2,x3 ! vars for polynomial solver
        complex(kind=dp)             :: qj(3) ! z-comp of wavector, k=[0,0,q], for qS1, qS2, qP waves
        real(kind=dp)                :: omega ! Wave angular velocity 
        
        nlm4 = nlm(:(I_l6-1)) ! l=0,2,4 coefficients
        
        omega = 1.0d0 ! Eigenvalues of problem are rho*V^2 where V=omega/k, so the problem does not, in fact, depend on omega. It is included here anyway for readability (and because of the way the analytical solutions were solved for in Mathematica).
        
        do nn = 1,size(theta_n)
            nlm4_rot = rotate_nlm4(nlm4, -theta_n(nn), -phi_n(nn)) ! negative angles because we are are rotating the specified direction (back) into the vertical orientation
            Qn = Qnorm(nlm4_rot, alpha, lam,mu,Elam,Emu,Egam) ! effective acoustic tensor
            ri = detQ_polycoefs(Qn, omega, rho) ! polynomial coefs of problem: det(Q) = 0 

            ! xi=0 solutions from depressed cubic r0 + r2*q^2 + r4*q^4 + r6*q^6 = 0 ==> x^3 + a*x + b = 0
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
            
            ! phase velocities are then
            vj(:,nn) = omega/abs(qj(:))
        end do
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

        call elas_revparams__tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5) ! (k1, ..., k5) are effective elastic coefficients, *not* related to wave vector.
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

        call elas_fwdparams__tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5) ! (k1, ..., k5) are effective elastic coefficients, *not* related to wave vector.

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

    !---------------------------------
    ! ELECTROMAGNETIC
    !---------------------------------

    ! N/A

end module wavepropagation 
