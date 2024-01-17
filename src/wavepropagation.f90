! Nicholas M. Rathmann <rathmann@nbi.ku.dk>, 2021-2024

module wavepropagation  

    use header
    use tensorproducts
    use mandel
    use elasticities
    use homogenizations
    use moments
    use rotation

    implicit none 

!    real(kind=dp), private, parameter :: k(3) = [0.0d0, 0.0d0, 1.0d0]
    real(kind=dp), private, parameter :: ksq(3,3) = reshape([0.0d0,0.0d0,0.0d0, 0.0d0,0.0d0,0.0d0, 0.0d0,0.0d0,1.0d0], [3,3]) ! k outer k, where k = [0,0,1]
       
contains      

    !-----------------------------------------------
    ! ELASTIC
    !-----------------------------------------------

    !---------------------------------
    ! Orthotropic grains
    !---------------------------------

    !--------------
    ! Crystallographic axis definitions
    !--------------
    ! qlm_1 := b
    ! qlm_2 := n
    ! qlm_3 := v
    !--------------
    ! Parameter definitions
    !--------------
    ! alpha:      Reuss--Voigt homogenization weight
    ! Lame_grain: Monocrystal Lame parameters (lam11,lam22,lam33, lam23,lam13,lam12, mu1,mu2,mu3) w.r.t. qi=(b,n,v) axes.
    !--------------

    function Vi_elastic_orthotropic(qlm_1,qlm_2,qlm_3, alpha,Lame_grain,rho, theta,phi) result (Vi)

        ! Elastic phase velocities (qP, qS1, qS2) given orientation distributions qlm_1, qlm_2, qlm_3 and propagation directions (theta, phi) 
        ! *** Only Voigt homogenization is supported (alpha is dummy variable for now) ***

        ! NOTE: Subscripts of Lame parameters (lam11,lam22,lam33, lam23,lam13,lam12, mu1,mu2,mu3) refer to the qi'th axis (i.e. distributions qlm_i). 
        ! E.g. lam11 is NOT the Lame parameter along m1', but q1(=b). So unless A-type fabrics are assumed, the entries in Lame_grain should be re-ordered.
        
        implicit none
        
        complex(kind=dp), intent(in) :: qlm_1(:), qlm_2(:), qlm_3(:)
        real(kind=dp), intent(in)    :: alpha, Lame_grain(9), rho ! grain parameters
        real(kind=dp), intent(in)    :: theta(:), phi(:) ! array of theta and phi values to calculate phase velocities (vi) along
        real(kind=dp)                :: Qnorm_Voigt(3,3)
        real(kind=dp)                :: Vi(3,size(theta)) ! qS1, qS2, qP phase velocities (in that order)
        complex(kind=dp)             :: qlm4_rot_1(L4len), qlm4_rot_2(L4len), qlm4_rot_3(L4len) ! truncated at L=4 and rotated
        real(kind=dp)                :: a2_i(3, 3,3), a4_ii(3, 3,3,3,3), a4_jk(3, 3,3,3,3)
        
        do nn = 1,size(theta)
            qlm4_rot_1 = rotate_nlm4(qlm_1(:L4len), -theta(nn), -phi(nn)) ! negative angles because we are are rotating the specified direction (back) into the vertical orientation
            qlm4_rot_2 = rotate_nlm4(qlm_2(:L4len), -theta(nn), -phi(nn)) 
            qlm4_rot_3 = rotate_nlm4(qlm_3(:L4len), -theta(nn), -phi(nn)) 
            
            call ai_orthotropic(qlm4_rot_1, qlm4_rot_2, qlm4_rot_3, a2_i, a4_ii, a4_jk)
            Qnorm_Voigt = ai_to_Qnorm_orthotropic(a2_i, a4_ii, a4_jk, Lame_grain)
            Vi(:,nn) = Qnorm_to_vi(Qnorm_Voigt, rho)
        end do
    end
    
    function Vi_elastic_orthotropic__discrete(mi, alpha,Lame_grain,rho, theta,phi) result (Vi)

        ! Elastic phase velocities (qP, qS1, qS2) given ensemble of crystallographic axes (m1', m2', m3') and propagation directions (theta, phi) 
        ! *** Only Voigt homogenization is supported (alpha is dummy variable for now) ***
                
        implicit none
        
        real(kind=dp), intent(in) :: mi(:,:,:) ! (3,3,N) = (m'_i, xyz comp., grain no.) 
        real(kind=dp), intent(in) :: alpha, Lame_grain(9), rho
        real(kind=dp), intent(in) :: theta(:), phi(:) ! arrays of theta and phi values to calculate phase velocities (vi) along
        real(kind=dp)             :: Qnorm_Voigt(3,3)
        real(kind=dp)             :: Vi(3,size(theta)) ! qS1, qS2, qP phase velocities (in that order)
        real(kind=dp)             :: mi_rot(3,3,size(mi,dim=3))
        real(kind=dp)             :: a2_i(3, 3,3), a4_ii(3, 3,3,3,3), a4_jk(3, 3,3,3,3)
        
        do nn = 1,size(theta)
            do jj = 1,size(mi,dim=3) ! grain no
                do ii = 1,3 ! m'_i
                    mi_rot(ii,:,jj) = rotate_vector(mi(ii,:,jj), -theta(nn), -phi(nn)) ! negative angles because we are are rotating the specified direction (back) into the vertical orientation
                end do
            end do

            call ai_orthotropic_discrete(mi_rot, a2_i, a4_ii, a4_jk)
            Qnorm_Voigt = ai_to_Qnorm_orthotropic(a2_i, a4_ii, a4_jk, Lame_grain)
            Vi(:,nn) = Qnorm_to_vi(Qnorm_Voigt, rho)
        end do
    end

    function ai_to_Qnorm_orthotropic(a2_i, a4_ii, a4_jk, Lame_grain) result (Qnorm)
       
        implicit none
        
        real(kind=dp), intent(in) :: Lame_grain(9), a2_i(3, 3,3), a4_ii(3, 3,3,3,3), a4_jk(3, 3,3,3,3)
        real(kind=dp)             :: lam(6), mu(3), Qnorm(3,3)

        ! Grain Lame parameters
        lam(:) = Lame_grain(1:6) 
        mu(:)  = Lame_grain(7:9)

        Qnorm = 0.0d0 ! initialize
        do ii=1,3
            Qnorm = Qnorm + mu(ii)/2  * (a2_i(ii,:,:) + anticommutator(a2_i(ii,:,:),ksq) + doubleinner22(a2_i(ii,:,:),ksq)*identity)
            Qnorm = Qnorm + lam(ii)   * doubleinner42(a4_ii(ii,:,:,:,:),ksq) ! 11=bb, 22=nn, 33=vv terms
            Qnorm = Qnorm + lam(ii+3) * doubleinner42_firstlast(2*a4_sym2(a4_jk(ii,:,:,:,:)),ksq) ! 23=nv, 13=bv, 12=nb terms
        end do
    end

    !---------------------------------
    ! Transversely isotropic grains
    !---------------------------------

    function Vi_elastic_tranisotropic(nlm, alpha, Lame_grain, rho, theta,phi) result (Vi)

        ! Elastic phase velocities (qP, qS1, qS2) given nlm and propagation directions (theta,phi) 
        ! Assumes transversely isotropic grains with parameters Lame_grain=(lam,mu,Elam,Emu,Egam)
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: alpha, Lame_grain(5), rho
        real(kind=dp), intent(in)    :: theta(:), phi(:) ! arrays of theta and phi values to calculate phase velocities (vi) along

        real(kind=dp)                :: Vi(3,size(theta)) ! qS1, qS2, qP phase velocities
        integer, parameter           :: n4mlen = nlm_lenvec(4)
        complex(kind=dp)             :: nlm4(n4mlen), nlm4_rot(n4mlen) ! nlm truncated at L=4, and its rotated form
        real(kind=dp)                :: Qnorm(3,3)
        
        nlm4 = nlm(:n4mlen) ! l=0,2,4 coefficients

        do nn = 1,size(theta)
            nlm4_rot = rotate_nlm4(nlm4, -theta(nn), -phi(nn)) ! negative angles because we are are rotating the specified direction (back) into the vertical orientation
            Qnorm    = Qnorm_tranisotropic(nlm4_rot, alpha, Lame_grain) ! effective acoustic tensor
            Vi(:,nn) = Qnorm_to_vi(Qnorm, rho)
        end do
    end

    function Qnorm_tranisotropic(nlm, alpha, Lame_grain) result (Qnorm)

        ! Normalized Voigt--Reuss averaged acoustic tensor such that Q = k^2 * Qnorm
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: alpha, Lame_grain(5) ! Reuss--Voigt weight (alpha) and monocrystal parameters
        real(kind=dp)                :: lam,mu,Elam,Emu,Egam ! unpacked Lame parameters
        real(kind=dp)                :: Qnorm(3,3), Qnorm_Reuss(3,3), Qnorm_Voigt(3,3)

        lam  = Lame_grain(1)
        mu   = Lame_grain(2)
        Elam = Lame_grain(3)
        Emu  = Lame_grain(4)
        Egam = Lame_grain(5)

        Qnorm_Reuss = 0.0d0
        Qnorm_Voigt = 0.0d0
        if (alpha .lt. 1.0d0) Qnorm_Reuss = Qnorm_tranisotropic_Reuss(nlm, lam,mu,Elam,Emu,Egam)
        if (alpha .gt. 0.0d0) Qnorm_Voigt = Qnorm_tranisotropic_Voigt(nlm, lam,mu,Elam,Emu,Egam)
        Qnorm = (1-alpha)*Qnorm_Reuss + alpha*Qnorm_Voigt
    end

    function Qnorm_tranisotropic_Voigt(nlm, lam,mu,Elam,Emu,Egam) result (Qnorm)

        ! Acoustic tensor, Qnorm=Qnorm_Voigt, such that div(stress) = k^2*matmul(Qnorm,u) where u is a plane wave in the strain field with k = [0,0,kz]
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: lam,mu,Elam,Emu,Egam ! Reuss--Voigt weight (alpha) and monocrystal parameters
        real(kind=dp)                :: Qnorm(3,3)
        real(kind=dp)                :: k1,k2,k3,k4,k5
        real(kind=dp)                :: a2mat(3,3), a2v(6), a4v(6,6), KM_anticomm(3,3)

        call elas_revparams_tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5) ! (k1, ..., k5) are effective elastic coefficients, *not* related to wave vector.
        call f_ev_ck_Mandel(nlm, a2v, a4v) ! Structure tensors in Mandel notation: a2v = ev_c2_Mandel, B = ev_c4_Mandel
        a2mat = vec_to_mat(a2v) ! = a2
        KM_anticomm = matmul(ksq,a2mat) + matmul(a2mat,ksq) ! = {a2,ksq} for k=kvert=[0,0,1]
        Qnorm = k2/2*identity + (k1+k2/2)*ksq + k3*KM_anticomm + k4*vec_to_mat(matmul(a4v,mat_to_vec(ksq))) + 0.5d0*k5*(a2mat + doubleinner22(a2mat,ksq)*identity + KM_anticomm)
    end

    function Qnorm_tranisotropic_Reuss(nlm, lam,mu,Elam,Emu,Egam) result(Qnorm)

        ! Acoustic tensor, Qnorm=Qnorm_Reuss, such that div(stress) = k^2*matmul(Qnorm,u) where u is a plane wave in the strain field with k = [0,0,kz]

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: lam,mu,Elam,Emu,Egam ! monocrystal parameters
        real(kind=dp)                :: Qnorm(3,3) 
        real(kind=dp)                :: k1,k2,k3,k4,k5 ! aux params (not related to wave vector)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        real(kind=dp)                :: ev_c2(3,3), ev_c4(3,3,3,3), P(9,9), Pinv(9,9) 

        call elas_fwdparams_tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5) ! (k1, ..., k5) are effective elastic coefficients, *not* related to wave vector.
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        ev_c2 = f_ev_c2(n00,n2m)
        ev_c4 = f_ev_c4(n00,n2m,n4m)
        include "include/P_Reuss.f90" ! Set P such that strain = matmul(P,stress). Requires ki, ev_c2, ev_c4
        Pinv = matinv(P)
        Qnorm(1,:) = 1/2.0d0 * [Pinv(3,3)+Pinv(3,7), Pinv(3,6)+Pinv(3,8), 2*Pinv(3,9)]
        Qnorm(2,:) = 1/2.0d0 * [Pinv(6,3)+Pinv(6,7), Pinv(6,6)+Pinv(6,8), 2*Pinv(6,9)]
        Qnorm(3,:) = 1/2.0d0 * [Pinv(9,3)+Pinv(9,7), Pinv(9,6)+Pinv(9,8), 2*Pinv(9,9)]
    end
        
    !-----------------------------------------------
    ! ELECTROMAGNETIC
    !-----------------------------------------------

    ! N/A

    !-----------------------------------------------
    ! SUPPORTING ROUTINES
    !-----------------------------------------------
        
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
