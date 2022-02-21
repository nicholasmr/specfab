! N. M. Rathmann <rathmann@nbi.ku.dk>, 2019-2022

module enhancementfactors  

    use tensorproducts
    use moments
    use homogenizations

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real,    parameter, private :: Pi = 3.141592653589793
    
    integer, private       :: Lcap
    real(kind=dp), private :: ev_c2_iso(3,3), ev_c4_iso(3,3,3,3), ev_c6_iso(3,3,3,3, 3,3), ev_c8_iso(3,3,3,3, 3,3,3,3) ! <c^k> for isotropic n(theta,phi)
    
contains      

    !---------------------------------
    ! INIT
    !---------------------------------
       
    subroutine initenhancementfactors(Lcap_)

        ! Needs to be called once before using the module routines.

        implicit none    
        integer, intent(in) :: Lcap_ ! Truncation "Lcap"
        complex(kind=dp)    :: nlm_iso(1+5+9+13+17) ! nlm for L=8 is needed
                
        Lcap = Lcap_ ! Save internal copy

        ! Calculate structure tensors for an isotropic fabric
        nlm_iso(:) = 0.0d0
        nlm_iso(1) = 1/Sqrt(4*Pi) ! Normalized ODF 
        call f_ev_ck(nlm_iso, 'f', ev_c2_iso,ev_c4_iso,ev_c6_iso,ev_c8_iso) ! Sets <c^i> := a^(i) for i=2,4,6,8
    end

    !---------------------------------
    ! ENHANCEMENT-FACTORS
    !---------------------------------

    function Evw(vw, tau, nlm, Ecc, Eca, alpha, nprime)

        ! Generalized enhancement factor, E_{vw}
        ! Assumes a transversely isotropic grain rheology (Ecc,Eca,alpha,nprime)

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: Ecc, Eca, alpha, vw(3,3), tau(3,3)
        integer, intent(in)          :: nprime
        real(kind=dp)                :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
        real(kind=dp)                :: Evw

        if (nprime .eq. 1) then 
            ! Linear grain rheology (n'=1) relies only on <c^2> and <c^4>
            call f_ev_ck(nlm, 'r', ev_c2,ev_c4,ev_c6,ev_c8) ! Calculate structure tensors of orders 2,4 (6,8 are assumed zero for faster evaluation)
        else if (nprime .eq. 3) then
            ! Nonlinear grain rheology (n'=3) relies on <c^k> for k=2,4,6,8
            call f_ev_ck(nlm, 'f', ev_c2,ev_c4,ev_c6,ev_c8) ! Calculate structure tensors of orders 2,4,6,8
        end if

        Evw = (1-alpha)*Evw_Sac(vw, tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca, nprime) &
                + alpha*Evw_Tay(vw, tau, ev_c2,ev_c4,             Ecc,Eca, nprime)
    end

    function Eeiej(nlm, e1,e2,e3, Ecc,Eca,alpha,nprime)

        ! Enhancement factors in directions (ei,ej), not necessarily a2 eigen directions
        ! (3x3 symmetric matrix of enhancement factors)

        implicit none

        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), dimension(3)  :: e1,e2,e3
        real(kind=dp), intent(in)    :: Ecc, Eca, alpha
        integer, intent(in)          :: nprime
        real(kind=dp)                :: Eeiej(3,3)
        
        ! Longitudinal
        Eeiej(1,1) = Evw(outerprod(e1,e1), tau_vv(e1),     nlm, Ecc,Eca,alpha,nprime) 
        Eeiej(2,2) = Evw(outerprod(e2,e2), tau_vv(e2),     nlm, Ecc,Eca,alpha,nprime)
        Eeiej(3,3) = Evw(outerprod(e3,e3), tau_vv(e3),     nlm, Ecc,Eca,alpha,nprime)    
        
        ! Shear
        Eeiej(1,2) = Evw(outerprod(e1,e2), tau_vw(e1,e2),  nlm, Ecc,Eca,alpha,nprime) 
        Eeiej(1,3) = Evw(outerprod(e1,e3), tau_vw(e1,e3),  nlm, Ecc,Eca,alpha,nprime) 
        Eeiej(2,3) = Evw(outerprod(e2,e3), tau_vw(e2,e3),  nlm, Ecc,Eca,alpha,nprime)
        
        ! Symmetric matrix
        Eeiej(2,1) = Eeiej(1,2)
        Eeiej(3,1) = Eeiej(1,3) 
        Eeiej(3,2) = Eeiej(2,3)   
    end

    function Evw_Sac(vw, tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca, nprime)

        implicit none
        
        real(kind=dp), intent(in) :: Ecc, Eca, vw(3,3), tau(3,3)
        integer, intent(in)       :: nprime
        real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
        real(kind=dp)             :: Evw_Sac
        
        Evw_Sac = doubleinner22(ev_epsprime_Sac(tau, ev_c2,    ev_c4,    ev_c6,    ev_c8,     Ecc,Eca,nprime), vw) / &
                  doubleinner22(ev_epsprime_Sac(tau, ev_c2_iso,ev_c4_iso,ev_c6_iso,ev_c8_iso, Ecc,Eca,nprime), vw)
    end

    function Evw_Tay(vw, tau, ev_c2,ev_c4, Ecc,Eca, nprime)

        implicit none
        
        real(kind=dp), intent(in) :: Ecc, Eca, vw(3,3), tau(3,3)
        integer, intent(in)       :: nprime
        real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3)
        real(kind=dp)             :: Evw_Tay
        
        Evw_Tay = doubleinner22(ev_epsprime_Tay(tau, ev_c2,    ev_c4,     Ecc,Eca,nprime), vw) / &
                  doubleinner22(ev_epsprime_Tay(tau, ev_c2_iso,ev_c4_iso, Ecc,Eca,nprime), vw)
    end

    !---------------------------------
    ! SYNTHETIC STRESS STATES
    !---------------------------------

    function tau_vv(v) 
        ! v--v compression/extension
        implicit none
        real(kind=dp), intent(in) :: v(3)
        integer, parameter :: identity(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
        real(kind=dp) :: tau_vv(3,3)
        tau_vv = identity/3.0d0 - outerprod(v,v)
    end

    function tau_vw(v,w)
        ! v--w shear
        implicit none
        real(kind=dp), intent(in) :: v(3), w(3)
        real(kind=dp) :: tau_vw(3,3)
        tau_vw = outerprod(v,w) + outerprod(w,v)
    end

end module enhancementfactors 
