! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

! The structure and declarations in this file might seem odd, but allow a simple wrapping of specfab.f90 into a Python module using F2PY

module specfabpy
    
    use specfab    

    implicit none

contains

    !------------------
    ! Setup
    !------------------

    subroutine init(Lcap_, nlmlen) 
        implicit none
        integer, intent(in) :: Lcap_
        integer, intent(out) :: nlmlen
        call initspecfab(Lcap_)
        nlmlen = nlm_len
    end
    
    function get_lm(nlmlen)
        implicit none
        integer, intent(in) :: nlmlen
        integer :: get_lm(2,nlmlen)
        get_lm(:,:) = lm(:,1:nlmlen)
    end
    
    !------------------
    ! Time evolution matrices
    !------------------
    
    function get_dndt_ij_lrot(nlmlen, eps,omg, tau,Aprime,Ecc,Eca, beta)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        integer, intent(in) :: nlmlen
        real(kind=dp), intent(in) :: eps(3,3), omg(3,3), tau(3,3) ! strain-rate (eps) and spin (omg)
        real(kind=dp), intent(in) :: Aprime, Ecc, Eca, beta
        complex(kind=dp) :: get_dndt_ij_lrot(nlmlen,nlmlen)
        get_dndt_ij_lrot = dndt_ij_LROT(eps,omg, tau,Aprime,Ecc,Eca, beta)
    end

    function get_dndt_ij_ddrx(nlmlen, tau)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        integer, intent(in) :: nlmlen
        real(kind=dp), intent(in) ::  tau(3,3)
        complex(kind=dp) :: get_dndt_ij_ddrx(nlmlen,nlmlen)
        get_dndt_ij_ddrx = dndt_ij_DDRX(tau)
    end
    
    function get_dndt_ij_regl(nlmlen, eps, tau, Aprime, beta)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        integer, intent(in) :: nlmlen
        real(kind=dp), intent(in) ::  eps(3,3), tau(3,3), Aprime, beta
        complex(kind=dp) :: get_dndt_ij_regl(nlmlen,nlmlen)
        get_dndt_ij_regl = dndt_ij_REGL(eps, tau,Aprime, beta) 
    end
    
    !------------------
    ! Eigen frame related
    !------------------
    
    subroutine get_eigenframe(nlm, e1,e2,e3, eigvals)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in)  :: nlm(:) ! f2py does not need to know the dimension
        real(kind=dp), dimension(3), intent(out) :: e1,e2,e3, eigvals
        call eigenframe(nlm, e1,e2,e3, eigvals)
    end
    
    subroutine get_pqframe(nlm, p23,p12,p13, q23,q12,q13)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), dimension(3), intent(out) :: p23,p12,p13, q23,q12,q13
        call pqframe(nlm, p23,p12,p13, q23,q12,q13)
    end
    
    subroutine get_mtframe_twodimensional(nlm, m,t, am,at, Emm,Emt, Exx,Exz, fabtype,  Ecc,Eca,alpha,nprime)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(out) :: m(2),t(2), am,at, Emm,Emt, Exx,Exz, fabtype ! (m,t) = two-dimension vectors (x,z coords)
        real(kind=dp), intent(in) :: Ecc, Eca, alpha
        integer, intent(in) :: nprime
        call mtframe_twodimensional(nlm, m,t, am,at, Emm,Emt, Exx,Exz, fabtype,  Ecc,Eca,alpha,nprime)
    end
    
    subroutine get_mtframe_threedimensional(nlm, m,t, am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype,  Ecc,Eca,alpha,nprime)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(out) :: m(3),t(3), am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype ! (m,t) = two-dimension vectors (x,z coords)
        real(kind=dp), intent(in) :: Ecc, Eca, alpha
        integer, intent(in) :: nprime
        call mtframe_threedimensional(nlm, m,t, am,at1,at2, Emm,Emt, Exx,Exy,Exz, fabtype,  Ecc,Eca,alpha,nprime)
    end
    
    !------------------
    ! Enhancement factors
    !------------------
    
    function get_Eeiej(nlm, Ecc,Eca,alpha,nprime)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in) :: Ecc, Eca, alpha
        integer, intent(in) :: nprime
        real(kind=dp) :: get_Eeiej(3,3)
        get_Eeiej = Eeiej(nlm, Ecc,Eca,alpha,nprime)
    end
    
    function get_Exixj(nlm, Ecc,Eca,alpha,nprime)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in) :: Ecc, Eca, alpha
        integer, intent(in) :: nprime
        real(kind=dp) :: get_Exixj(3,3)
        get_Exixj = Exixj(nlm, Ecc,Eca,alpha,nprime)
    end

    function get_Epijqij(nlm, Ecc,Eca,alpha,nprime)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in) :: Ecc, Eca, alpha
        integer, intent(in) :: nprime    
        real(kind=dp) :: get_Epijqij(3)
        get_Epijqij = Epijqij(nlm, Ecc,Eca,alpha,nprime)
    end
    
    function get_Evw(ev_c2,ev_c4,ev_c6,ev_c8, vw,tau, Ecc,Eca,alpha,nprime)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        real(kind=dp), intent(in) :: Ecc, Eca, alpha, vw(3,3), tau(3,3)
        integer, intent(in) :: nprime
        real(kind=dp) :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
        real(kind=dp) :: get_Evw
        get_Evw = (1-alpha)*Evw_Sac(vw, tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca, nprime) + alpha*Evw_Tay(vw, tau, ev_c2,ev_c4, Ecc,Eca, nprime)
    end
    
    !------------------
    ! Rheology
    !------------------
    
    function get_tau_of_eps(eps,m,A,Emm,Emt,n)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        real(kind=dp), intent(in) :: eps(3,3),m(3),A,Emm,Emt
        integer, intent(in) :: n
        real(kind=dp) :: get_tau_of_eps(3,3)
        get_tau_of_eps = tau_of_eps(eps,m,A,Emm,Emt,n)
    end
    
    function get_eps_of_tau(tau,m,A,Emm,Emt,n)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        real(kind=dp), intent(in) :: tau(3,3),m(3),A,Emm,Emt
        integer, intent(in) :: n
        real(kind=dp) :: get_eps_of_tau(3,3)
        get_eps_of_tau = eps_of_tau(tau,m,A,Emm,Emt,n)
    end
    
    !------------------
    ! AUX
    !------------------
    
    function get_nu(eps)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        real(kind=dp), intent(in) :: eps(3,3)
        real(kind=dp) :: get_nu
        get_nu = f_nu(eps)
    end
    
    subroutine get_ck_moments(nlm, ev_c2,ev_c4,ev_c6,ev_c8)
        implicit none
        integer, parameter :: dp = 8 ! Default precision
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(out) :: ev_c2(3,3),ev_c4(3,3,3,3),ev_c6(3,3,3,3, 3,3),ev_c8(3,3,3,3, 3,3,3,3)
        call ck_moments(nlm, ev_c2,ev_c4,ev_c6,ev_c8)
    end
        
end module specfabpy 

