! N. M. Rathmann <rathmann@nbi.ku.dk>, 2021

program demo

    use specfab    
    
    implicit none

    integer, parameter :: dp = 8
    real(kind=dp)      :: A, Eij(3,3), m1(3),m2(3),m3(3)
    integer            :: n
    real(kind=dp), dimension(3), parameter :: x1 = [1,0,0], x2 = [0,1,0], x3 = [0,0,1] ! x,y,z dir.

    call initspecfab(4) ! Not going to use specfab, but is good practice to initialize it even when using static routines (as is the case here).

    !-------------------------------------------------------------------
    ! E_ij
    !-------------------------------------------------------------------

    ! Set a synthetic enhancement-factor (fluidity) structure that we later want to recover (test self consistency)
   
    Eij(1,1) = 5  ! m1--m1 longitudinal enhancement 
    Eij(2,2) = 10 ! m2--m2 longitudinal enhancement 
    Eij(3,3) = 20 ! m3--m3 longitudinal enhancement 

    Eij(2,3) = 50 ! m2--m3 shear enhancement 
    Eij(3,1) = 60 ! m3--m1 shear enhancement 
    Eij(1,2) = 70 ! m1--m2 shear enhancement 

    ! Not strictly needed in this demo, but Eij should be symmetric
    Eij(3,2) = Eij(2,3) 
    Eij(1,3) = Eij(3,1)
    Eij(2,1) = Eij(1,2)
    
    !-------------------------------------------------------------------
    ! m_i
    !-------------------------------------------------------------------

    ! Set orthonormal triad of symmetry axes
    m1 = x1 
    m2 = x2 
    m3 = x3 

    ! Alternative (45 deg rotated)
!    m1 = [+1.0d0/sqrt(2.0d0), 0.0d0, 1.0d0/sqrt(2.0d0)] 
!    m2 = [-1.0d0/sqrt(2.0d0), 0.0d0, 1.0d0/sqrt(2.0d0)] 
!    m3 = [0.0d0, 1.0d0, 0.0d0]
    
    ! Alternative (45 deg rotated)
!    m1 = [+1.0d0/sqrt(2.0d0), 0.0d0, 1.0d0/sqrt(2.0d0)] 
!    m2 = [0.0d0, 1.0d0, 0.0d0]
!    m3 = [-1.0d0/sqrt(2.0d0), 0.0d0, 1.0d0/sqrt(2.0d0)] 
    
    !-------------------------------------------------------------------
    ! Test orthtropic viscosity structure is correct
    !-------------------------------------------------------------------
    
    A = 3.5d0 ! Rate factor (arbitrary value)
    n = 3     ! Flow law exponent

    print *, '--------------------------------------------------------'    
    print *, 'Test that enhancement factors are correctly recovered for longitudinal and shear deformation along m_1, m_2, m_3'
    print *, '--------------------------------------------------------'
        
    print *, '(i,j)=(1,1) :: eps_{m_i m_j}/eps_{m_i m_j}^{Glen} = ', eps_ratio(tau_vv(m1),    A,n, m1,m2,m3, Eij, m1,m1), ' --- should be E_{ij} = ', Eij(1,1)
    print *, '(i,j)=(2,2) :: eps_{m_i m_j}/eps_{m_i m_j}^{Glen} = ', eps_ratio(tau_vv(m2),    A,n, m1,m2,m3, Eij, m2,m2), ' --- should be E_{ij} = ', Eij(2,2)
    print *, '(i,j)=(3,3) :: eps_{m_i m_j}/eps_{m_i m_j}^{Glen} = ', eps_ratio(tau_vv(m3),    A,n, m1,m2,m3, Eij, m3,m3), ' --- should be E_{ij} = ', Eij(3,3)
    print *, '(i,j)=(2,3) :: eps_{m_i m_j}/eps_{m_i m_j}^{Glen} = ', eps_ratio(tau_vw(m2,m3), A,n, m1,m2,m3, Eij, m2,m3), ' --- should be E_{ij} = ', Eij(2,3)
    print *, '(i,j)=(3,1) :: eps_{m_i m_j}/eps_{m_i m_j}^{Glen} = ', eps_ratio(tau_vw(m1,m3), A,n, m1,m2,m3, Eij, m3,m1), ' --- should be E_{ij} = ', Eij(3,1)
    print *, '(i,j)=(1,2) :: eps_{m_i m_j}/eps_{m_i m_j}^{Glen} = ', eps_ratio(tau_vw(m1,m2), A,n, m1,m2,m3, Eij, m1,m2), ' --- should be E_{ij} = ', Eij(1,2)
   
    print *, "... where m_1, m_2, m_3 = "
    print *, m1
    print *, m2
    print *, m3
   
    !-------------------------------------------------------------------
    ! Test inverse rheology --> check that eps_in = eps(tau(eps_in))
    !-------------------------------------------------------------------

    print *, '--------------------------------------------------------'    
    print *, 'Test self consistency of forward and inverse rheologies (verify that tau0 = tau(eps(tau0)) for different stress tensors)'
    print *, '--------------------------------------------------------'
    print *, '...compression and shear stress tensors w.r.t m_1, m_2, and m_3:'
    call tau_of_eps_of_tau(tau_vv(x1), A, n, m1,m2,m3, Eij) 
    call tau_of_eps_of_tau(tau_vv(x2), A, n, m1,m2,m3, Eij) 
    call tau_of_eps_of_tau(tau_vv(x3), A, n, m1,m2,m3, Eij) 
    call tau_of_eps_of_tau(tau_vw(x1,x2), A, n, m1,m2,m3, Eij) 
    call tau_of_eps_of_tau(tau_vw(x1,x3), A, n, m1,m2,m3, Eij) 
    call tau_of_eps_of_tau(tau_vw(x2,x3), A, n, m1,m2,m3, Eij) 
    print *, '...an arbitrary, non-trivial stress tensor:'
    call tau_of_eps_of_tau(0.1*tau_vw(x1,x2) + 0.33*tau_vw(x2,x3) + 0.51*tau_vw(x1,x2) + 0.43*tau_vv(x3), A, n, m1,m2,m3, Eij) 
   
contains

function eps_ratio(tau,A,n, m1,m2,m3, Eij, v,w) 

    implicit none

    integer, parameter        :: dp = 8
    real(kind=dp), intent(in) :: tau(3,3), v(3), w(3)
    real(kind=dp), intent(in) :: A, Eij(3,3), m1(3),m2(3),m3(3)
    integer, intent(in)       :: n
    real(kind=dp)             :: eps_ratio

    eps_ratio = doubleinner22(eps_of_tau__orthotropic(tau,A,n, m1,m2,m3, Eij), outerprod(v,w)) / &
                doubleinner22(eps_of_tau__isotropic(  tau,A,n),                outerprod(v,w))
end

subroutine tau_of_eps_of_tau(tau_in, A, n, m1,m2,m3, Eij) 

    implicit none

    integer, parameter        :: dp = 8
    real(kind=dp), intent(in) :: tau_in(3,3), A, m1(3),m2(3),m3(3), Eij(3,3)
    integer, intent(in)       :: n
    real(kind=dp)             :: eps(3,3), tau(3,3)

    eps = eps_of_tau__orthotropic(tau_in, A, n, m1,m2,m3, Eij)
    tau = tau_of_eps__orthotropic(eps,    A, n, m1,m2,m3, Eij)
    
    print *, 'tau0             = ', tau_in
    print *, 'tau(eps(tau0))   = ', tau
    print *, 'norm2(abs(diff)) = ', norm2(abs(tau-tau_in))
    print *, '----------------'
        
end
    
end program

