! N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

program demo

    use specfab    
    
    implicit none

    integer, parameter :: dp = 8
    real(kind=dp)      :: mu,lam, eta1,eta2,eta3, m(3)
    real(kind=dp)      :: Elam,Emu,Egam, eta1_check,eta2_check,eta3_check ! for verifying conversion between aniso params etai <--> Ei

    ! For testing Reuss homogenization
    complex(kind=dp)   :: nlm(15)
    complex(kind=dp)   :: nlmiso(15)     
    real(kind=dp)      :: a2tensor(3,3), a4tensor(3,3,3,3), strain_input(3,3), stress_true(3,3), stress_homo(3,3)
    
    real(kind=dp), dimension(3), parameter :: x1 = [1,0,0], x2 = [0,1,0], x3 = [0,0,1] ! x,y,z dir.
    
    call initspecfab(4) ! Not going to use specfab, but is good practice to initialize it even when using static routines (as is the case here).
   
    ! Monocrystal parameters
    mu   = mu_Bennet68
    lam  = lam_Bennet68
    eta1 = eta1_Bennet68
    eta2 = eta2_Bennet68
    eta3 = eta3_Bennet68
    call elastic_tranisotropic__etai_to_Ei(lam,mu,eta1,eta2,eta3, Elam,Emu,Egam)
    call elastic_tranisotropic__Ei_to_etai(lam,mu,Elam,Emu,Egam, eta1_check,eta2_check,eta3_check)
    print *, 'Anisotropic params (Bennet 1968): eta1, eta2, eta3 = ', eta1, eta2, eta3
    print *, 'Anisotropic params (from etai):   Elam, Emu, Egam  = ', Elam, Emu, Egam
    print *, 'Anisotropic params (from Ei):     eta1, eta2, eta3 = ', eta1_check, eta2_check, eta3_check

    ! Symmetry axis 
    m = x1
    m = [+1.0d0/sqrt(2.0d0), 0.0d0, 1.0d0/sqrt(2.0d0)] ! Alternative (45 deg rotated)
   
   
    !-------------------------------------------------------------------
    ! Test inverse rheology --> check that eps_in = eps(tau(eps_in))
    !-------------------------------------------------------------------

    print *, '--------------------------------------------------------'    
    print *, 'Test self consistency of forward and inverse elasticities (verify that tau0 = tau(eps(tau0)) for different stress tensors)'
    print *, '--------------------------------------------------------'
    print *, '...compression and shear stress tensors w.r.t m and t:'
    call tau_of_eps_of_tau(tau_vv(x1), lam,mu,Elam,Emu,Egam,m)  
    call tau_of_eps_of_tau(tau_vv(x2), lam,mu,Elam,Emu,Egam,m) 
    call tau_of_eps_of_tau(tau_vv(x3), lam,mu,Elam,Emu,Egam,m) 
    call tau_of_eps_of_tau(tau_vw(x1,x2), lam,mu,Elam,Emu,Egam,m) 
    call tau_of_eps_of_tau(tau_vw(x1,x3), lam,mu,Elam,Emu,Egam,m) 
    call tau_of_eps_of_tau(tau_vw(x2,x3), lam,mu,Elam,Emu,Egam,m) 
    print *, '...an arbitrary, non-trivial stress tensor:'
    call tau_of_eps_of_tau(0.1*tau_vw(x1,x2) + 0.33*tau_vw(x2,x3) + 0.51*tau_vw(x1,x2) + 0.43*tau_vv(x3) + 2*outerprod(x1,x1), lam,mu,Elam,Emu,Egam,m) ! tau does not have to be traceless, here we also add a (last) term breaking that
   
   
    !-------------------------------------------------------------------
    ! Test Reuss homogenization
    !-------------------------------------------------------------------
    
    print *, '--------------------------------------------------------'    
    print *, 'Test Reuss homogenization by considering a perfect single maximum'
    print *, '--------------------------------------------------------'
    
    ! ODF = delta(r-z) 
    a4tensor = 0.0d0 ! init

!    m = x1 
!    a4tensor(1,1,1,1) = 1 ! m = x1

    m = x2 
    a4tensor(2,2,2,2) = 1 ! m = x2

!    m = x3 
!    a4tensor(3,3,3,3) = 1 ! m = x3

    nlm = a4_to_nlm(a4tensor)
!    print *, 'nlm=', nlm

    strain_input = tau_vv(x3)
!    strain_input = tau_vw(x1,x3)
    print *, 'strain_input = ', strain_input
    stress_true = stress_of_strain__tranisotropic(strain_input, lam,mu,Elam,Emu,Egam,m)
    stress_homo = ev_sigprime_Reuss(strain_input, nlm, lam,mu,Elam,Emu,Egam)
   
    print *, 'stress_true = ', stress_true
    print *, 'stress_homo = ', stress_homo
    print *, 'norm2(abs(diff)) [%] = ', 100*norm2(abs(stress_homo-stress_true))/norm2(abs(stress_true))
    print *, '----------------'
    
    !-------------------------------------------------------------------
    ! Test acoustic tensor, Q
    !-------------------------------------------------------------------
    
    print *, ' '
    print *, '--------------------------------------------------------'    
    print *, 'Test acoustic tensor, Q, agrees between Reuss and Voigt approximations in limits where they should match (isotropy, of perfect single max)'
    print *, '--------------------------------------------------------'

    ! init
    a2tensor = 0.0d0 
    a4tensor = 0.0d0 
    nlm(:) = 0.0d0
    
    !---- Fabric = isotropic -----

    print *, '---------------------'    
    print *, '*** ISOTROPIC FABRIC'
    print *, '---------------------'
    
    nlmiso = 0.0d0
    nlmiso(1) = 1/sqrt(4*3.14159265359)
 
    print *, '::: lam,mu,Ei=1 :::'
    call test_Q(nlmiso, 1.0d0,1.0d0,1.0d0,1.0d0,1.0d0) 

    print *, '::: Ei=1 :::'
    call test_Q(nlmiso, lam,mu,1.0d0,1.0d0,1.0d0) 
    
    print *, '::: full aniso params :::'
    call test_Q(nlmiso, lam,mu,Elam,Emu,Egam) 
    
    !---- Fabric = single maximum -----

    print *, '---------------------'   
    print *, '*** ANISOTROPIC FABRIC'
    print *, '---------------------'
       
!    a4tensor(1,1,1,1) = 1 ! m = x1
!    a4tensor(2,2,2,2) = 1 ! m = x2
!!    a4tensor(3,3,3,3) = 1 ! m = x3
!    nlm = a4_to_nlm(a4tensor)

    a2tensor(2,2) = 0.5 
    a2tensor(3,3) = 0.5
!    a2tensor(3,3) = 1
    nlm(1:6) = a2_to_nlm(a2tensor)

    print *, '::: lam,mu,Ei=1 :::'
    call test_Q(nlm, 1.0d0,1.0d0,1.0d0,1.0d0,1.0d0) 

    print *, '::: Ei=1 :::'
    call test_Q(nlm, lam,mu,1.0d0,1.0d0,1.0d0) 

    print *, '::: full aniso params :::'
    call test_Q(nlm, lam,mu,Elam,Emu,Egam) 

contains

subroutine tau_of_eps_of_tau(tau_in, lam,mu,Elam,Emu,Egam,m) 

    implicit none

    integer, parameter        :: dp = 8
    real(kind=dp), intent(in) :: tau_in(3,3), lam,mu,Elam,Emu,Egam,m(3)
    real(kind=dp)             :: eps(3,3), tau(3,3)

    eps = strain_of_stress__tranisotropic(tau_in,  lam,mu,Elam,Emu,Egam,m)
    tau = stress_of_strain__tranisotropic(eps,     lam,mu,Elam,Emu,Egam,m)
    
    print *, 'tau0             = ', tau_in
    print *, 'tau(eps(tau0))   = ', tau
    print *, 'norm2(abs(diff)) = ', norm2(abs(tau-tau_in))
    print *, '----------------'
        
end

subroutine test_Q(nlm, lam,mu,Elam,Emu,Egam) 

    implicit none

    integer, parameter           :: dp = 8
    complex(kind=dp), intent(in) :: nlm(15)
    real(kind=dp), intent(in)    :: lam,mu,Elam,Emu,Egam
    real(kind=dp)                :: Qr(3,3), Qv(3,3)

    Qr = Qnorm_Reuss(nlm, lam,mu, Elam,Emu,Egam)
    Qv = Qnorm_Voigt(nlm, lam,mu, Elam,Emu,Egam)
    
    print *, 'alpha = 1 (Voigt): Q=', Qv
    print *, 'alpha = 0 (Reuss): Q=', Qr
    print *, 'norm2(abs(diff)) = ', norm2(abs(Qv-Qr)), ' or in pct = ', 100*norm2(abs(Qr-Qv))/norm2(abs(Qv))
    print *, '----------------'
end
    
end program

