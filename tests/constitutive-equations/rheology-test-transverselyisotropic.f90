! N. M. Rathmann <rathmann@nbi.ku.dk>, 2021-2022

program demo

    use specfab    
    
    implicit none

!    integer, parameter :: dp = 8
    real(kind=dp)      :: A,n, Eij(2), m(3), t(3)
    real(kind=dp), dimension(3), parameter :: x1 = [1,0,0], x2 = [0,1,0], x3 = [0,0,1] ! x,y,z dir.

    call initspecfab(4) ! Not going to use specfab, but is good practice to initialize it even when using static routines (as is the case here).

    !-------------------------------------------------------------------
    ! E_ij
    !-------------------------------------------------------------------

    ! Set a synthetic enhancement-factor (fluidity) structure that we later want to recover (test self consistency)
   
    Eij = [0.5, 10.0] ! (Emm, Emt)
    
    !-------------------------------------------------------------------
    ! m_i
    !-------------------------------------------------------------------

    ! Set orthonormal triad of symmetry axes
    m = x1 
    t = x2 

    ! Alternative (45 deg rotated)
    m = [+1.0d0/sqrt(2.0d0), 0.0d0, 1.0d0/sqrt(2.0d0)] 
    t = [-1.0d0/sqrt(2.0d0), 0.0d0, 1.0d0/sqrt(2.0d0)] 
    
    ! Alternative (45 deg rotated)
!    m = [+1.0d0/sqrt(2.0d0), 0.0d0, 1.0d0/sqrt(2.0d0)] 
!    t = [0.0d0, 1.0d0, 0.0d0]
    
    !-------------------------------------------------------------------
    ! Test viscosity structure is correct
    !-------------------------------------------------------------------
    
    A = 3.5d0 ! Rate factor (arbitrary value)
    n = 3.0d0 ! Flow law exponent

    print *, '--------------------------------------------------------'    
    print *, 'Test that enhancement factors are correctly recovered for longitudinal and shear deformation along m and t'
    print *, '--------------------------------------------------------'
        
    print *, 'eps_{mm}/eps_{mm}^{Glen} = ', eps_ratio(tau_vv(m),   A,n, m,Eij, m,m), ' --- should be E_{mm} = ', Eij(1)
    print *, 'eps_{mt}/eps_{mt}^{Glen} = ', eps_ratio(tau_vw(m,t), A,n, m,Eij, m,t), ' --- should be E_{mt} = ', Eij(2)
   
    print *, "... where m, t = "
    print *, m
    print *, t
   
    !-------------------------------------------------------------------
    ! Test inverse rheology --> check that eps_in = eps(tau(eps_in))
    !-------------------------------------------------------------------

    print *, '--------------------------------------------------------'    
    print *, 'Test self consistency of forward and inverse rheologies (verify that tau0 = tau(eps(tau0)) for different stress tensors)'
    print *, '--------------------------------------------------------'
    print *, '...compression and shear stress tensors w.r.t m and t:'
    call tau_of_eps_of_tau(tau_vv(x1), A,n, m,Eij) 
    call tau_of_eps_of_tau(tau_vv(x2), A,n, m,Eij) 
    call tau_of_eps_of_tau(tau_vv(x3), A,n, m,Eij) 
    call tau_of_eps_of_tau(tau_vw(x1,x2), A,n, m,Eij) 
    call tau_of_eps_of_tau(tau_vw(x1,x3), A,n, m,Eij) 
    call tau_of_eps_of_tau(tau_vw(x2,x3), A,n, m,Eij) 
    print *, '...an arbitrary, non-trivial stress tensor:'
    call tau_of_eps_of_tau(0.1*tau_vw(x1,x2) + 0.33*tau_vw(x2,x3) + 0.51*tau_vw(x1,x2) + 0.43*tau_vv(x3), A,n, m,Eij) 
   
contains

function eps_ratio(tau, A,n, m,Eij, v,w) 

    implicit none

    integer, parameter        :: dp = 8
    real(kind=dp), intent(in) :: tau(3,3), v(3), w(3)
    real(kind=dp), intent(in) :: A,n, m(3), Eij(2)
    real(kind=dp)             :: eps_ratio

    eps_ratio = doubleinner22(rheo_fwd_tranisotropic(tau, A,n, m,Eij), outerprod(v,w)) / &
                doubleinner22(rheo_fwd_isotropic(    tau, A,n),            outerprod(v,w))
end

subroutine tau_of_eps_of_tau(tau_in, A,n, m,Eij) 

    implicit none

    integer, parameter        :: dp = 8
    real(kind=dp), intent(in) :: tau_in(3,3), A,n, m(3), Eij(2)
    real(kind=dp)             :: eps(3,3), tau(3,3)

    eps = rheo_fwd_tranisotropic(tau_in,  A,n, m,Eij)
    tau = rheo_rev_tranisotropic(eps,     A,n, m,Eij)
    
    print *, 'tau0             = ', tau_in
    print *, 'tau(eps(tau0))   = ', tau
    print *, 'norm2(abs(diff)) = ', norm2(abs(tau-tau_in))
    print *, '----------------'
        
end
    
end program

