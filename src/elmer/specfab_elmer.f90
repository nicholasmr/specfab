! D. A. Lilien and N. M. Rathmann, 2019-2022

!-------------------
! AUX
!-------------------

function vectorize9(M) result(vec)
    implicit none
    real(kind=dp) :: M(3,3)
    real(kind=dp) :: vec(9)
    vec = reshape(M, [9])
end

function outerprod_to_Mandel(A,B) result (M)

    implicit none
    real(kind=dp), intent(in) :: A(3,3), B(3,3)
    real(kind=dp)             :: M(6,6)
    real(kind=dp), parameter  :: s = sqrt(2.0d0) 

    M(1,:) = [  A(1,1)*B(1,1),   A(1,1)*B(2,2),   A(1,1)*B(3,3), s*A(1,1)*B(2,3), s*A(1,1)*B(1,3), s*A(1,1)*B(1,2)]
    M(2,:) = [  A(2,2)*B(1,1),   A(2,2)*B(2,2),   A(2,2)*B(3,3), s*A(2,2)*B(2,3), s*A(2,2)*B(1,3), s*A(2,2)*B(1,2)]
    M(3,:) = [  A(3,3)*B(1,1),   A(3,3)*B(2,2),   A(3,3)*B(3,3), s*A(3,3)*B(2,3), s*A(3,3)*B(1,3), s*A(3,3)*B(1,2)]
    M(4,:) = [s*A(2,3)*B(1,1), s*A(2,3)*B(2,2), s*A(2,3)*B(3,3), 2*A(2,3)*B(2,3), 2*A(2,3)*B(1,3), 2*A(2,3)*B(1,2)]
    M(5,:) = [s*A(1,3)*B(1,1), s*A(1,3)*B(2,2), s*A(1,3)*B(3,3), 2*A(1,3)*B(2,3), 2*A(1,3)*B(1,3), 2*A(1,3)*B(1,2)]
    M(6,:) = [s*A(1,2)*B(1,1), s*A(1,2)*B(2,2), s*A(1,2)*B(3,3), 2*A(1,2)*B(2,3), 2*A(1,2)*B(1,3), 2*A(1,2)*B(1,2)]
end

include "elmer/include/IBOF.f90"

!-------------------
! VECTORIZED RHEOLOGIES
!-------------------
    
function Cmat_inverse_orthotropic(eps, A,n, m1,m2,m3, Eij) result(C)

    ! Vectorized inverse orthotropic flow law "tau_of_eps__orthotropic()".
    ! Returns 9x9 matrix "C" such that vec(tau) = matmul(C, vec(eps)), where vec(tau) and vec(eps) are 9x1 column vectors.

    implicit none
    real(kind=dp), intent(in)     :: A, m1(3),m2(3),m3(3), Eij(3,3)
    integer, intent(in)           :: n
    real(kind=dp)                 :: eps(3,3), C(9,9)
    real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
    real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
    real(kind=dp)                 :: J1,J2,J3,J4,J5,J6, J23,J31,J12
    real(kind=dp)                 :: viscosity

    call orthotropic_coefs(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)
    call orthotropic_tensors_and_invars(eps, m1,m2,m3, M11,M22,M33,M23,M31,M12, J1,J2,J3,J4,J5,J6)
    call orthotropic_auxinvars(J1,J2,J3, J23,J31,J12)

    viscosity = A**(-1.d0/n) * ( &
        + lam1/gam * J23**2 &
        + lam2/gam * J31**2 & 
        + lam3/gam * J12**2 &
        + 4 * (1/lam4) * J4**2 &
        + 4 * (1/lam5) * J5**2 &
        + 4 * (1/lam6) * J6**2 &
    )**((1-n)/(2.d0*n))

    C = viscosity * ( &
        -3/4.0d0 * lam1/gam * outerprod9(vectorize9(identity - 3*M11), vectorize9(M11)) &
        -3/4.0d0 * lam2/gam * outerprod9(vectorize9(identity - 3*M22), vectorize9(M22)) &
        -3/4.0d0 * lam3/gam * outerprod9(vectorize9(identity - 3*M33), vectorize9(M33)) &
        + 2/lam4            * outerprod9(vectorize9(M23 + transpose(M23)), vectorize9(M23)) &
        + 2/lam5            * outerprod9(vectorize9(M31 + transpose(M31)), vectorize9(M31)) &
        + 2/lam6            * outerprod9(vectorize9(M12 + transpose(M12)), vectorize9(M12)) &
    )
end

function Cmandel_inverse_orthotropic(eps, A,n, m1,m2,m3, Eij) result(Cmandel)

    ! Vectorized inverse orthotropic flow law "tau_of_eps__orthotropic()" in Mandel's form.
    
    ! Returns 6x6 matrix "Cmandel" such that: 
    !       mat_to_vec(tau) = matmul(Cmandel, mat_to_vec(eps))
    ! where 
    !       mat_to_vec(H) = [H_xx, H_yy, H_zz, sqrt(2)*H_yz, sqrt(2)*H_xz, sqrt(2)*H_xy] 
    ! is a 6x1 Mandel vector (see mandel.f90).
    
    ! To get the usual 3x3 form of tau, use: tau = vec_to_mat( matmul(Cmandel, mat_to_vec(eps)) )

    implicit none
    real(kind=dp), intent(in)     :: A, m1(3),m2(3),m3(3), Eij(3,3)
    integer, intent(in)           :: n
    real(kind=dp)                 :: eps(3,3), Cmandel(6,6)
    real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
    real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
    real(kind=dp)                 :: J1,J2,J3,J4,J5,J6, J23,J31,J12
    real(kind=dp)                 :: viscosity

    call orthotropic_coefs(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)
    call orthotropic_tensors_and_invars(eps, m1,m2,m3, M11,M22,M33,M23,M31,M12, J1,J2,J3,J4,J5,J6)
    call orthotropic_auxinvars(J1,J2,J3, J23,J31,J12)

    viscosity = A**(-1.d0/n) * ( &
        + lam1/gam * J23**2 &
        + lam2/gam * J31**2 & 
        + lam3/gam * J12**2 &
        + 4 * (1/lam4) * J4**2 &
        + 4 * (1/lam5) * J5**2 &
        + 4 * (1/lam6) * J6**2 &
    )**((1.0d0-n)/(2.d0*n))

    Cmandel = viscosity * ( &
        -3/4.0d0 * lam1/gam * outerprod_to_Mandel(identity-3*M11,M11) &
        -3/4.0d0 * lam2/gam * outerprod_to_Mandel(identity-3*M22,M22) &
        -3/4.0d0 * lam3/gam * outerprod_to_Mandel(identity-3*M33,M33) &
        + 1/lam4 * outerprod_to_Mandel(M23+transpose(M23), M23+transpose(M23)) &
        + 1/lam5 * outerprod_to_Mandel(M31+transpose(M31), M31+transpose(M31)) &
        + 1/lam6 * outerprod_to_Mandel(M12+transpose(M12), M12+transpose(M12)) &
    )
end

subroutine Cmat_inverse_orthotropic_dimless(eps, n, m1,m2,m3, Eij, MinInVar, viscosity, C6)

    ! Vectorized inverse orthotropic flow law "tau_of_eps__orthotropic()" in Mandel's form.
    
    ! Returns 6x6 matrix "Cmandel" such that: 
    !       mat_to_vec(tau) = matmul(Cmandel, mat_to_vec(eps))
    ! where 
    !       mat_to_vec(H) = [H_xx, H_yy, H_zz, sqrt(2)*H_yz, sqrt(2)*H_xz, sqrt(2)*H_xy] 
    ! is a 6x1 Mandel vector (see mandel.f90).
    
    ! To get the usual 3x3 form of tau, use: tau = vec_to_mat( matmul(Cmandel, mat_to_vec(eps)) )

    implicit none
    real(kind=dp), intent(in)     :: m1(3),m2(3),m3(3), Eij(3,3), MinInVar
    real(kind=dp), intent(out)    :: viscosity, C6(6,6) 
    integer, intent(in)           :: n
    real(kind=dp)                 :: eps(3,3), Cmandel(6,6)
    real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
    real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
    real(kind=dp)                 :: J1,J2,J3,J4,J5,J6, J23,J31,J12

    call orthotropic_coefs(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)
    call orthotropic_tensors_and_invars(eps, m1,m2,m3, M11,M22,M33,M23,M31,M12, J1,J2,J3,J4,J5,J6)
    call orthotropic_auxinvars(J1,J2,J3, J23,J31,J12)

    viscosity = ( &
        + lam1/gam * J23**2 &
        + lam2/gam * J31**2 & 
        + lam3/gam * J12**2 &
        + 4 * (1/lam4) * J4**2 &
        + 4 * (1/lam5) * J5**2 &
        + 4 * (1/lam6) * J6**2 &
    )

    IF (viscosity.LT.MinInVar) viscosity = MinInVar
    viscosity = (2.0_dp * viscosity) ** ((1.0_dp-n)/(2.0_dp*n))

    Cmandel = ( &
        -3/4.0d0 * lam1/gam * outerprod_to_Mandel(identity-3*M11,M11) &
        -3/4.0d0 * lam2/gam * outerprod_to_Mandel(identity-3*M22,M22) &
        -3/4.0d0 * lam3/gam * outerprod_to_Mandel(identity-3*M33,M33) &
        + 1/lam4 * outerprod_to_Mandel(M23+transpose(M23), M23+transpose(M23)) &
        + 1/lam5 * outerprod_to_Mandel(M31+transpose(M31), M31+transpose(M31)) &
        + 1/lam6 * outerprod_to_Mandel(M12+transpose(M12), M12+transpose(M12)) &
    )
    C6 = CMandel([1, 2, 3, 6, 4, 5], [1, 2, 3, 6, 4, 5])
    C6(1:3, 1:3) = 2.0_dp * C6(1:3, 1:3)
    C6(4:6, 1:3) = C6(4:6, 1:3) / sqrt(2.0_dp)
    C6(1:3, 4:6) = C6(1:3, 4:6) * sqrt(2.0_dp)
end

!-------------------
! STRUCTURE TENSORS
!-------------------

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
!    integer :: i,j,k,l

    a4 = 0.0 ! init
    include "elmer/include/ae4_to_a4__body.f90"
end

