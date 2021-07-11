! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2019-2021

module specfabpy_const

    implicit none
    integer, parameter :: dp = 8 ! Default precision
    real, parameter    :: Pi = 3.1415927
    integer, parameter :: x = 1, y = 2, z = 3 ! Matrix indices
    
end module specfabpy_const


module specfabpy
    
    use specfabpy_const
    use specfab, &
        a2__sf => a2, &
        a4__sf => a4, &
        a2_to_nlm__sf => a2_to_nlm, &
        a4_to_nlm__sf => a4_to_nlm, &
        tau_of_eps__sf => tau_of_eps, &
        eps_of_tau__sf => eps_of_tau, &
        frame__sf => frame, &
        Eeiej__sf => Eeiej, &
        Evw__sf => Evw, &
        Eca_opt_lin__sf    => Eca_opt_lin, &
        Ecc_opt_lin__sf    => Ecc_opt_lin, &
        alpha_opt_lin__sf  => alpha_opt_lin, &
        Eca_opt_nlin__sf   => Eca_opt_nlin, &
        Ecc_opt_nlin__sf   => Ecc_opt_nlin, &
        alpha_opt_nlin__sf => alpha_opt_nlin
        
    implicit none
    
    ! Optimal n'=1 (lin) grain parameters 
    real(kind=dp), parameter :: Eca_opt_lin   = Eca_opt_lin__sf
    real(kind=dp), parameter :: Ecc_opt_lin   = Ecc_opt_lin__sf
    real(kind=dp), parameter :: alpha_opt_lin = alpha_opt_lin__sf
    ! Optimal n'=3 (nlin) grain parameters 
    real(kind=dp), parameter :: Eca_opt_nlin   = Eca_opt_nlin__sf
    real(kind=dp), parameter :: Ecc_opt_nlin   = Ecc_opt_nlin__sf
    real(kind=dp), parameter :: alpha_opt_nlin = alpha_opt_nlin__sf

contains

    !---------------------------------
    ! SETUP
    !---------------------------------

    subroutine init(L, nlmlen) 
        implicit none
        integer, intent(in)  :: L
        integer, intent(out) :: nlmlen
        
        call initspecfab(L)
        nlmlen = nlm_len
    end
    
    function get_lm(nlmlen)
        implicit none
        integer, intent(in) :: nlmlen
        integer             :: get_lm(2,nlmlen)
        
        get_lm(:,:) = lm(:,1:nlmlen)
    end
        
    !---------------------------------
    ! FABRIC DYNAMICS
    !---------------------------------
    
    function dndt_LATROT(nlm, eps,omg) 
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: eps(3,3), omg(3,3)
        complex(kind=dp)             :: dndt_LATROT(size(nlm),size(nlm))
        
        dndt_LATROT = dndt_ij_LATROT(eps,omg, 0*eps,0d0,1d0,1d0, 1d0) ! only the Taylor model is provided for the python interface
    end
    
    function dndt_DDRX(nlm, tau)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: tau(3,3)
        complex(kind=dp)             :: dndt_DDRX(size(nlm),size(nlm))
        
        dndt_DDRX = dndt_ij_DDRX(nlm, tau)
    end

    function dndt_DDRX_src(nlmlen, tau)
        use specfabpy_const
        implicit none
        integer, intent(in)       :: nlmlen
        real(kind=dp), intent(in) :: tau(3,3)
        complex(kind=dp)          :: dndt_DDRX_src(nlmlen,nlmlen)
        
        dndt_DDRX_src = dndt_ij_DDRX_src(tau)
    end
    
    function dndt_REG(nlm)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        complex(kind=dp)             :: dndt_REG(size(nlm),size(nlm))
        
        dndt_REG = dndt_ij_REG() 
    end    
    
    function nu(nu0, eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: nu0, eps(3,3)
        real(kind=dp)             :: nu
        
        nu = f_nu_eps(nu0, eps)
    end

    !---------------------------------
    ! FABRIC FRAME
    !---------------------------------
    
    subroutine frame(nlm, ftype, e1,e2,e3, eigvals)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:) 
        character*1, intent(in)      :: ftype 
        real(kind=dp), intent(out)   :: e1(3),e2(3),e3(3), eigvals(3)
        
        call frame__sf(nlm, ftype, e1,e2,e3, eigvals)
    end
    
    ! JOSEF (Rathmann and Lilien, 2021) fabric interface   
!    include "frame_josef.f90"
    
    !------------------
    ! Enhancement factors
    !------------------
    
    function enhfac_vw(a2,a4,a6,a8, vw,tau, Ecc,Eca,alpha,nprime) result(Evw)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: vw(3,3), tau(3,3), Ecc,Eca,alpha
        integer, intent(in)       :: nprime
        real(kind=dp)             :: a2(3,3), a4(3,3,3,3), a6(3,3,3,3, 3,3), a8(3,3,3,3, 3,3,3,3)
        real(kind=dp)             :: Evw
        
        Evw = (1-alpha)*Evw_Sac(vw, tau, a2,a4,a6,a8, Ecc,Eca,nprime) &
                + alpha*Evw_Tay(vw, tau, a2,a4,       Ecc,Eca,nprime)
    end
    
    function enhfac_eiej(nlm, e1,e2,e3, Ecc,Eca,alpha,nprime) result(Eeiej)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: e1(3),e2(3),e3(3)
        real(kind=dp), intent(in)    :: Ecc, Eca, alpha
        integer, intent(in)          :: nprime
        real(kind=dp)                :: Eeiej(3,3)
        
        Eeiej = Eeiej__sf(nlm, e1,e2,e3, Ecc,Eca,alpha,nprime)
    end

    !---------------------------------
    ! STRUCTURE TENSORS
    !---------------------------------

    function a2(nlm) 
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a2(3,3)
      
        a2 = a2__sf(nlm)
    end    
    
    function a4(nlm) 
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a4(3,3,3,3)

        a4 = a4__sf(nlm)
    end    
    
    subroutine ai(nlm, a2,a4,a6,a8)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(out)   :: a2(3,3), a4(3,3,3,3), a6(3,3,3,3, 3,3), a8(3,3,3,3, 3,3,3,3)

        call f_ev_ck(nlm, a2,a4,a6,a8) ! recall notation: ai := a^(i) := ev_ci := <c^i> 
    end
    
    function a2_to_nlm(a2) 
        use specfabpy_const
        implicit none
        integer, parameter        :: nlmlen = 1+5
        real(kind=dp), intent(in) :: a2(3,3)
        real(kind=dp)             :: a2_to_nlm(nlmlen)

        a2_to_nlm = a2_to_nlm__sf(a2)
    end
        
    function a4_to_nlm(a2,a4)
        use specfabpy_const
        implicit none
        integer, parameter        :: nlmlen = 1+5+9
        real(kind=dp), intent(in) :: a2(3,3), a4(3,3,3,3)
        real(kind=dp)             :: a4_to_nlm(nlmlen)
        
        a4_to_nlm = a4_to_nlm__sf(a2,a4)
    end
        
    !---------------------------------
    ! TRANSVERSELY ISOTROPIC RHEOLOGY 
    !---------------------------------
    
    function tau_of_eps(eps, m,A,Emm,Emt,n)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), m(3),A,Emm,Emt 
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: tau_of_eps(3,3)
        
        tau_of_eps = tau_of_eps__sf(eps, m,A,Emm,Emt,n)
    end
    
    function eps_of_tau(tau, m,A,Emm,Emt,n)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), m(3),A,Emm,Emt
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: eps_of_tau(3,3)
        
        eps_of_tau = eps_of_tau__sf(tau, m,A,Emm,Emt,n)
    end
        
    !---------------------------------
    ! AUX
    !---------------------------------
    
    function DDRX_decayrate(nlm, tau)
    
        use specfabpy_const
        implicit none
        
        integer, parameter           :: lmDDRX_len = 1+5+9 ! Scope of harmonic interactions is local in wave space (this is NOT a free parameter)
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: tau(3,3) ! DDRX rate contant, dev. stress tensor
        complex(kind=dp)             :: DDRX_decayrate(lmDDRX_len)
        !...
        complex(kind=dp) :: g(lmDDRX_len)
        real(kind=dp)    :: Davg
        complex(kind=dp) :: qt(-2:2)
        real(kind=dp)    :: k

        ! Quadric expansion coefficients of tau
        qt = quad_rr(tau)
        
        ! Harmonic interaction weights 
        include "include/DDRX__body.f90"
        g = k*g/doubleinner22(tau,tau) ! = D
        
        Davg = doubleinner22(matmul(tau,tau), a2__sf(nlm)) - doubleinner22(tau,doubleinner42(a4__sf(nlm),tau)) ! (tau.tau):a2 - tau:a4:tau 
        Davg = Davg/doubleinner22(tau,tau) ! = <D>
        
        ! Davg is the "isotropic reference to be subtracted"; that is, we subtract a constant value from "D" (on S^2), which amounts to adjusting the value of the isotropic speactral expansion coefficient (l,m=0,0). 
        ! The factor of Sqrt(4*Pi) ensures that the isotropic contribution is indeed n00*Y00 = Sqrt(4*Pi)*Davg*Y00 = Davg
        g(1) = g(1) - Sqrt(4*Pi)*Davg 
        
        ! Calculate D - <D>
        DDRX_decayrate = g 
    end
    
    !---------------------------------
    ! EXTERNAL, NON-CORE FEATURES
    !---------------------------------

    ! Elmer ice flow model 
    include "elmer/specfabpy_elmer.f90"

    ! JOSEF ice flow model 
    include "josef/specfabpy_josef.f90"
        
end module specfabpy 

