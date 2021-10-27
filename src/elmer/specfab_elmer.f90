include "elmer/include/IBOF.f90"

function outerprod9(a,b) result(outerprod)
    ! a_i b_j = rank-2 tensor
    implicit none
    real(kind=dp), intent(in) :: a(9), b(9)
    real(kind=dp) :: outerprod(9,9)
!    outerprod = reshape( [( [( a(ii)*b(jj), jj=1,9)], ii=1,9)] , [9,9])
    outerprod = reshape( [( [( a(ii)*b(jj), ii=1,9)], jj=1,9)] , [9,9])
end

function vectorize9(M) result(vec)
    implicit none
    real(kind=dp) :: M(3,3)
    real(kind=dp) :: vec(9)
    vec = reshape(M, [9])
end
    
function Cmat_inverse_orthotropic(eps, A,n, m1,m2,m3, Eij) result(C)

    ! Vectorized inverse orthotropic flow law "tau_of_eps__orthotropic()".
    ! Returns 9x9 matrix "C" such that vec(tau) = matmul(C, vec(eps)), where vec(tau) and vec(eps) are 9x1 column vectors.

    implicit none
    real(kind=dp), intent(in)     :: A, m1(3),m2(3),m3(3), Eij(3,3)
    integer, intent(in)           :: n
    real(kind=dp)                 :: eps(3,3), C(9,9)
    real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
    real(kind=dp)                 :: eps11,eps22,eps33,eps23,eps31,eps12
    real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
    real(kind=dp)                 :: J1,J2,J3,J4,J5,J6, J23,J31,J12
    real(kind=dp)                 :: viscosity

    call orthotropic_coefs(n, Eij, lam1,lam2,lam3,lam4,lam5,lam6, gam)
    call orthotropic_tensors_and_invars(eps, m1,m2,m3, M11,M22,M33,M23,M31,M12, eps11,eps22,eps33,eps23,eps31,eps12)
    
    J1 = (eps22 - eps33)/2.0d0 
    J2 = (eps33 - eps11)/2.0d0 
    J3 = (eps11 - eps22)/2.0d0 
    J4 = eps23 ! J4 := (eps23 + eps32)/2 (but eps is symmetric)
    J5 = eps31 ! J5 := (eps31 + eps13)/2 (but eps is symmetric)
    J6 = eps12 ! J6 := (eps12 + eps21)/2 (but eps is symmetric)
    
    J23 = -3/2.0d0 * eps11 ! J23 := J2-J3 = (eps22 + eps33 - 2*eps11)/2.0d0 = -3/2.0d0 * eps11
    J31 = -3/2.0d0 * eps22 ! J31 := J3-J1 = (eps11 + eps33 - 2*eps22)/2.0d0 = -3/2.0d0 * eps22
    J12 = -3/2.0d0 * eps33 ! J12 := J1-J2 = (eps11 + eps22 - 2*eps33)/2.0d0 = -3/2.0d0 * eps33

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

