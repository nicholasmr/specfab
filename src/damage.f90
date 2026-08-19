! N. M. Rathmann <rathmann@nbi.ku.dk> 2026-

! Damage dynamics

module damage  

    use header
    use tensorproducts
    use moments ! used by tensorial dynamics routines
    use gaunt

    implicit none 

!    integer, parameter, private :: x = 1, y = 2, z = 3 ! Matrix indices (used for readability)
!    complex(kind=dp), parameter, private :: r = (1,0), i = (0,1) ! real and imag units

    integer, private :: Lcap    ! Truncation "L" of expansion series (internal copy of what was passed to the init routine).
    integer, private :: nlm_len ! Total number of expansion coefficients (i.e. DOFs)

    ! Static (constant) matrices used for spectral dynamics
!    real(kind=dp), parameter :: Ldiag(nlm_lenmax) = [( (-ll*(ll+1),mm=-ll,ll), ll=0,  Lcap__max,2)] ! Diagonal entries of Laplacian diffusion operator.
    
    integer, parameter :: SHI_DAMAGE=1+5
    integer, parameter :: SHI_L4=1+5+9
    
    
contains      

    !---------------------------------
    ! INIT
    !---------------------------------

    subroutine initdamage(Lcap_)

        ! Needs to be called once before using the module routines.

        implicit none    
        integer, intent(in) :: Lcap_ ! Truncation "Lcap"
        
        ! Save internal copy
        Lcap    = Lcap_ 
        nlm_len = (Lcap+1)*(Lcap+2)/2 

        ! Set gaunt coefficients (overlap integrals involving three spherical harmonics)
        call set_gaunts()
    end
    

    function cart2mag(A) result(Alm)
        ! to magnetic basis
        implicit none
        real(kind=dp), intent(in) :: A(3,3)
        complex(kind=dp)          :: Alm(SHI_DAMAGE)
        include "include/cart2mag.f90"
    end

    function mag2cart(Alm) result(A)
        ! to magnetic basis
        implicit none
        complex(kind=dp), intent(in) :: Alm(SHI_DAMAGE)
        real(kind=dp)                :: A(3,3)
        include "include/mag2cart.f90"
    end


    subroutine eig_sym_3x3(A, V, eval)
        ! Eigen-decomposition of symmetric 3x3 matrix
        implicit none
        real(kind=dp), intent(in)  :: A(3,3)
        real(kind=dp), intent(out) :: V(3,3), eval(3) ! eigenvectors (columns), eigenvalues
        integer, parameter :: lwork = 20 ! >= 3*n-1 = 8 for n=3 (matrix dimensions)
        real(kind=dp)      :: work(lwork), Acopy(3,3)
        integer            :: info

        Acopy = A
        V = Acopy
        call dsyev('V','U',3,V,3,eval,work,lwork,info)
        if (info /= 0) then
            stop "DSYEV failed in eig_sym_3x3"
        end if
        
    end
    

    function matpos(A) result(Ap)
        ! Positive (tensile) part
        implicit none
        real(kind=dp), intent(in) :: A(3,3)
        real(kind=dp)             :: Ap(3,3), V(3,3), eval(3), lam
        integer                   :: ii,jj,kk

        call eig_sym_3x3(A, V, eval)
        Ap = 0.0d0
        do ii = 1,3
           lam = max(eval(ii), 0.0d0)
           do jj = 1,3
              do kk = 1,3
                 Ap(jj,kk) = Ap(jj,kk) + lam * V(jj,ii) * V(kk,ii)
              end do
           end do
        end do
    end

    function matneg(A) result(Am)
        ! Negative (compressive) part
        implicit none
        real(kind=dp), intent(in) :: A(3,3)
        real(kind=dp)             :: Am(3,3), V(3,3), eval(3), lam
        integer                   :: ii,jj,kk

        call eig_sym_3x3(A, V, eval)
        Am = 0.0d0
        do ii = 1,3
           lam = min(eval(ii), 0.0d0)
           do jj = 1,3
              do kk = 1,3
                 Am(jj,kk) = Am(jj,kk) + lam * V(jj,ii) * V(kk,ii)
              end do
           end do
        end do
    end

    function damageffect(D) result(G)
        ! Computes damage effect tensor G = (I - D)^(-1) using spectral decomposition
        implicit none
        real(kind=dp), intent(in) :: D(3,3)
        real(kind=dp)             :: G(3,3)
        real(kind=dp) :: V(3,3), eval(3), lam, denom
        integer       :: i,j,k

        call eig_sym_3x3(D, V, eval)

        G = 0.0d0

        do i = 1,3
            lam = eval(i)
            denom = 1.0d0 - lam

!            print *, lam
            if (abs(denom) < 1.0d-3) then
                stop "Matrix (I - D) is singular or nearly singular"
            end if

            do j = 1,3
                do k = 1,3
                   G(j,k) = G(j,k) + (1.0d0/denom) * V(j,i) * V(k,i)
                end do
            end do
        end do
    end

    function netstresstensor(stress, D)
        ! Net stress tensor (sigma tilde)
        implicit none
        real(kind=dp), intent(in) :: stress(3,3), D(3,3) ! stress tensor, damage tensor
        real(kind=dp)             :: G(3,3), netstresstensor(3,3)

        ! Proper
        G = damageffect(D)
        netstresstensor = (matmul(stress, G) + matmul(G, stress))/2
        
!        ! Linearization
!        netstresstensor = stress + (matmul(stress,D)+matmul(D,stress))/2
    end

    function g_SI(stress, D) result(g)
        ! Harmonic couplig coefficients for 
        ! mode I fracture suseptibility (normalized)
        implicit none
        real(kind=dp), intent(in) :: stress(3,3), D(3,3) ! stress tensor, damage tensor
        complex(kind=dp)          :: Alm(SHI_DAMAGE), g(SHI_L4)
        real(kind=dp)             :: stresspos(3,3), netstresspos(3,3), k, norm

        netstresspos = matpos(netstresstensor(stress,D))
        Alm = cart2mag(netstresspos) ! magnetic basis
        g(:) = 0.0d0 ! initialize g before setting components
        include "include/SI-coupling-weights.f90"
        stresspos = matpos(stress)
        norm = 2.0d0/15 * doubleinner22(stresspos,stresspos)
        if (norm > 1e-16) g(:) = k*g(:)/norm
    end
    

    function g_SII(stress, D) result(g)
        ! Harmonic couplig coefficients for 
        ! mode II fracture suseptibility (normalized)
        implicit none
        real(kind=dp), intent(in) :: stress(3,3), D(3,3) ! stress tensor, damage tensor
        complex(kind=dp)          :: Alm(SHI_DAMAGE), g(SHI_L4)
        real(kind=dp)             :: netstress(3,3), k, norm

        netstress = netstresstensor(stress,D)
        Alm = cart2mag(netstress) ! magnetic basis
        g(:) = 0.0d0 ! initialize g before setting components
        include "include/SII-coupling-weights.f90"
        norm = 1.0d0/5 * doubleinner22(stress,stress)
        if (norm > 1e-16) g(:) = k*g(:)/norm
    end
    
    function damagetensor(nlm) result(D)    
        implicit none
        complex(kind=dp), intent(in) :: nlm(nlm_len)
        complex(kind=dp)             :: nlm_temp(nlm_len)
        real(kind=dp)                :: D(3,3)
    
        nlm_temp = nlm 
        nlm_temp(1) = 1
        D = a2(nlm_temp) * REAL(nlm(1)) ! equal to <n^2>, but not normalized
    end
    
    function M_FRSU_I(nlm, stress) result(M)
    
        ! Dynamical operator for mode I fracture susceptibility
    
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(nlm_len)
        real(kind=dp), intent(in)    :: stress(3,3)
        complex(kind=dp)             :: M(nlm_len,nlm_len)
        real(kind=dp)                :: D(3,3), stresspos(3,3), norm, integrity
        complex(kind=dp)             :: g(SHI_L4)
        integer                      :: ii
        
        D = damagetensor(nlm)
        g  = g_SI(stress,D)

        stresspos = matpos(stress)
        norm = 2.0d0/15 * doubleinner22(stresspos,stresspos)

        integrity = 1 - (D(1,1)+D(2,2)+D(3,3)) / 1 ! integrity

        M = 0.0d0
        do ii = 1,SHI_L4
            M(ii,ii) = integrity * norm * g(ii) ! user should multiply operator by S0/critstress
        end do
        
    end
    
end module damage
