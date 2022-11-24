! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2022

! This file is a wrapper for the Fortran routines to be provided in specfabpy

module specfabpy_const

    implicit none
    integer, parameter :: dp = 8 ! Default precision
    real, parameter    :: Pi = 3.141592653589793
    integer, parameter :: x = 1, y = 2, z = 3 ! Matrix indices
    
end module specfabpy_const

module specfabpy
    
    use specfabpy_const
    
    use specfab, &
    
        ! Fabric processes
        M_LROT__sf => M_LROT, &
        M_DDRX__sf => M_DDRX, M_DDRX_src__sf => M_DDRX_src, & 
        M_CDRX__sf => M_CDRX, &
        M_REG__sf  => M_REG, &

        ! Structure tensors 
        a2__sf => a2, a4__sf => a4, & 
        a2_to_nlm__sf => a2_to_nlm, a4_to_nlm__sf => a4_to_nlm, & ! Canonical representations
        ae2_to_a2__sf => ae2_to_a2, ae4_to_a4__sf => ae4_to_a4, & ! Elmer representations
        a4_to_mat__sf => a4_to_mat, & ! a4 in Mandel notation
        a4_eigentensors__sf => a4_eigentensors, & ! 3x3 eigentensors of a4
        a4_IBOF__sf => a4_IBOF, & ! From Elmer        
        
        ! Rheologies
        tau_of_eps__isotropic__sf     => tau_of_eps__isotropic,     eps_of_tau__isotropic__sf     => eps_of_tau__isotropic, &
        tau_of_eps__tranisotropic__sf => tau_of_eps__tranisotropic, eps_of_tau__tranisotropic__sf => eps_of_tau__tranisotropic, &
        tau_of_eps__orthotropic__sf   => tau_of_eps__orthotropic,   eps_of_tau__orthotropic__sf   => eps_of_tau__orthotropic, &
        tau_of_eps__orthotropic_Martin__sf => tau_of_eps__orthotropic_Martin, eps_of_tau__orthotropic_Martin__sf => eps_of_tau__orthotropic_Martin, &
        tau_of_eps__orthotropic_Pettit__sf => tau_of_eps__orthotropic_Pettit, eps_of_tau__orthotropic_Pettit__sf => eps_of_tau__orthotropic_Pettit, &
        eps_of_tau__linearTaylorSachs__sf => eps_of_tau__linearTaylorSachs, &
        
        ! Elasticities 
        stress_of_strain__tranisotropic__sf => stress_of_strain__tranisotropic, &
        strain_of_stress__tranisotropic__sf => strain_of_stress__tranisotropic, &
        Vi_elastic_tranisotropic__sf => Vi_elastic_tranisotropic, &
        Qnorm__sf => Qnorm, &
        
        ! Fluid enhancement factors
        frame__sf => frame, &
        Eeiej__sf => Eeiej, &
        Evw__sf   => Evw, &

        ! Optimal transversely isotropic (monocrystal) fluid enhancement factors for ice grains 
        Eca_opt_lin__sf  => Eca_opt_lin,  Ecc_opt_lin__sf  => Ecc_opt_lin,  alpha_opt_lin__sf  => alpha_opt_lin, &
        Eca_opt_nlin__sf => Eca_opt_nlin, Ecc_opt_nlin__sf => Ecc_opt_nlin, alpha_opt_nlin__sf => alpha_opt_nlin, &
        
        ! nlm and rnlm 
        lm__sf => lm, &
        nlm_len__sf => nlm_len, rnlm_len__sf => rnlm_len, &
        nlm_to_rnlm__sf => nlm_to_rnlm, rnlm_to_nlm__sf => rnlm_to_nlm, &
        reduce_M__sf => reduce_M, &
        rotate_nlm__sf => rotate_nlm, &
        Sl__sf => Sl, & ! Power spectrum
        
        ! Numerics
        Ldiag__sf => Ldiag, &
        apply_bounds__sf => apply_bounds, & 

        ! Fabric processes
        Gamma0__sf => Gamma0

        
    implicit none
    
    ! Optimal n'=1 (lin) grain parameters 
    real(kind=dp), parameter :: Eca_opt_lin   = Eca_opt_lin__sf
    real(kind=dp), parameter :: Ecc_opt_lin   = Ecc_opt_lin__sf
    real(kind=dp), parameter :: alpha_opt_lin = alpha_opt_lin__sf
    
    ! Optimal n'=3 (nlin) grain parameters 
    real(kind=dp), parameter :: Eca_opt_nlin   = Eca_opt_nlin__sf
    real(kind=dp), parameter :: Ecc_opt_nlin   = Ecc_opt_nlin__sf
    real(kind=dp), parameter :: alpha_opt_nlin = alpha_opt_nlin__sf
    
    integer :: Lcap

contains

    !---------------------------------
    ! SETUP
    !---------------------------------

    subroutine init(L, lm, nlm_len) 
        implicit none
        integer, intent(in)  :: L
        integer, intent(out) :: lm(2,(L+1)*(L+2)/2)
        integer, intent(out) :: nlm_len
        
        Lcap = L
        call initspecfab(L)
        nlm_len = nlm_len__sf
        lm = lm__sf(:,1:nlm_len)
    end

    !---------------------------------
    ! FABRIC DYNAMICS
    !---------------------------------
    
    function M_LROT(nlm, eps,omg) 
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: eps(3,3), omg(3,3)
        complex(kind=dp)             :: M_LROT(size(nlm),size(nlm))
        
        M_LROT = M_LROT__sf(eps,omg, 0*eps,0d0,1d0,1d0, 1d0) ! only the Taylor model is provided for the python interface
    end
    
    function M_DDRX(nlm, tau)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: tau(3,3)
        complex(kind=dp)             :: M_DDRX(size(nlm),size(nlm))
        
        M_DDRX = M_DDRX__sf(nlm, tau)
    end

    function M_DDRX_src(nlmlen, tau)
        use specfabpy_const
        implicit none
        integer, intent(in)       :: nlmlen
        real(kind=dp), intent(in) :: tau(3,3)
        complex(kind=dp)          :: M_DDRX_src(nlmlen,nlmlen)
        
        M_DDRX_src = M_DDRX_src__sf(tau)
    end
    
    function M_CDRX(nlm)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        complex(kind=dp)             :: M_CDRX(size(nlm),size(nlm))
        
        M_CDRX = M_CDRX__sf()
    end
    
    function M_REG(nlm, eps)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: eps(3,3)
        real(kind=dp)                :: M_REG(size(nlm),size(nlm))
        
        M_REG = M_REG__sf(eps) 
    end    
    
    function Lmat(nlm)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: Lmat(size(nlm),size(nlm))
        integer                      :: ii
        
        Lmat = 0.0
        do ii = 1, size(nlm)
            Lmat(ii,ii) = Ldiag__sf(ii) ! Laplacian (diagonal) matrix
        end do
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
    
    subroutine a4_eigentensors(nlm, Q1,Q2,Q3,Q4,Q5,Q6, eigvals)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:) 
        real(kind=dp), intent(out)   :: Q1(3,3),Q2(3,3),Q3(3,3),Q4(3,3),Q5(3,3),Q6(3,3)
        real(kind=dp), intent(out)   :: eigvals(6)
        
        call a4_eigentensors__sf(nlm, Q1,Q2,Q3,Q4,Q5,Q6, eigvals)
    end
        
    !------------------
    ! ENHANCEMENT FACTORS
    !------------------
    
    function Evw(nlm, vw,tau, Ecc,Eca,alpha,nprime) 
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in) :: vw(3,3), tau(3,3), Ecc,Eca,alpha
        integer, intent(in)       :: nprime
        real(kind=dp)             :: Evw
        
        Evw = (1-alpha)*Evw_Sac(vw, tau, nlm, Ecc,Eca,nprime) &
                + alpha*Evw_Tay(vw, tau, nlm, Ecc,Eca,nprime)
    end
    
    function Eeiej(nlm, e1,e2,e3, Ecc,Eca,alpha,nprime) 
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

        call f_ev_ck(nlm, 'f', a2,a4,a6,a8) ! recall notation: ai := a^(i) := ev_ci := <c^i> 
    end
    
    function a2_to_nlm(a2) result(nlm)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: a2(3,3)
        integer, parameter        :: nlmlen = 1+5
        complex(kind=dp)          :: nlm(nlmlen)

        nlm = a2_to_nlm__sf(a2)
    end
        
    function a4_to_nlm(a4) result(nlm)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: a4(3,3,3,3)
        integer, parameter        :: nlmlen = 1+5+9
        complex(kind=dp)          :: nlm(nlmlen)
        
        nlm = a4_to_nlm__sf(a4)
    end
    
   function a4_to_mat(A) result (M)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: A(3,3,3,3)
        real(kind=dp)             :: M(6,6)
        
        M = a4_to_mat__sf(A)
    end
        
    !---------------------------------
    ! GRAIN-AVERAGED RHEOLOGY (SACHS, TAYLOR)
    !---------------------------------
       
    function eps_of_tau__linearTaylorSachs(tau, nlm, Aprime,Ecc,Eca,alpha) result(eps_of_tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in)    :: tau(3,3), Aprime, Ecc,Eca,alpha
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: eps_of_tau(3,3)
    
        eps_of_tau = eps_of_tau__linearTaylorSachs__sf(tau, nlm, Aprime,Ecc,Eca,alpha)
    end
        
    !---------------------------------
    ! TRANSVERSELY ISOTROPIC RHEOLOGY 
    !---------------------------------

    function eps_of_tau__tranisotropic(tau, A,n, m,Emm,Emt) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, m(3),Emm,Emt
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: eps(3,3)
        
        eps = eps_of_tau__tranisotropic__sf(tau, A,n, m,Emm,Emt)
    end
    
    function tau_of_eps__tranisotropic(eps, A,n, m,Emm,Emt) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, m(3),Emm,Emt 
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: tau(3,3)
        
        tau = tau_of_eps__tranisotropic__sf(eps, A,n, m,Emm,Emt)
    end
    
    !---------------------------------
    ! ORTHOTROPIC RHEOLOGY 
    !---------------------------------
    
    function eps_of_tau__orthotropic(tau, A,n, m1,m2,m3, Eij) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, m1(3),m2(3),m3(3), Eij(3,3)
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: eps(3,3)
        
        eps = eps_of_tau__orthotropic__sf(tau, A,n, m1,m2,m3, Eij)
    end
    
    function tau_of_eps__orthotropic(eps, A,n, m1,m2,m3, Eij) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, m1(3),m2(3),m3(3), Eij(3,3)
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: tau(3,3)
        
        tau = tau_of_eps__orthotropic__sf(eps, A,n, m1,m2,m3, Eij)
    end
    
    ! Pettit's hypothesis that the fluidity is orientation independent.
    function eps_of_tau__orthotropic_Pettit(tau, A,n, m1,m2,m3, Eij) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, m1(3),m2(3),m3(3), Eij(3,3)
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: eps(3,3)
        
        eps = eps_of_tau__orthotropic_Pettit__sf(tau, A,n, m1,m2,m3, Eij)
    end
    
    function tau_of_eps__orthotropic_Pettit(eps, A,n, m1,m2,m3, Eij) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, m1(3),m2(3),m3(3), Eij(3,3)
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: tau(3,3)
        
        tau = tau_of_eps__orthotropic_Pettit__sf(eps, A,n, m1,m2,m3, Eij)
    end
    
    ! Martin's hypothesis that the viscosity is orientation independent.
    function eps_of_tau__orthotropic_Martin(tau, A,n, m1,m2,m3, Eij) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, m1(3),m2(3),m3(3), Eij(3,3)
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: eps(3,3)
        
        eps = eps_of_tau__orthotropic_Martin__sf(tau, A,n, m1,m2,m3, Eij)
    end
    
    function tau_of_eps__orthotropic_Martin(eps, A,n, m1,m2,m3, Eij) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, m1(3),m2(3),m3(3), Eij(3,3)
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: tau(3,3)
        
        tau = tau_of_eps__orthotropic_Martin__sf(eps, A,n, m1,m2,m3, Eij)
    end
    
    !---------------------------------
    ! ISOTROPIC (GLEN) RHEOLOGY 
    !---------------------------------

    function eps_of_tau__isotropic(tau, A,n) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: eps(3,3)
        
        eps = eps_of_tau__isotropic__sf(tau, A,n)
    end
    
    function tau_of_eps__isotropic(eps, A,n) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A
        integer, intent(in)       :: n ! Glen exponent
        real(kind=dp)             :: tau(3,3)
        
        tau = tau_of_eps__isotropic__sf(eps, A,n)
    end
    
    !---------------------------------
    ! TRANSVERSELY ISOTROPIC ELASTICITY 
    !---------------------------------
    
    function stress_of_strain__tranisotropic(strain, lam,mu, Elam,Emu,Egam,m) result(stress)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: strain(3,3), lam,mu, Elam,Emu,Egam,m(3)
        real(kind=dp)             :: stress(3,3)
        
        stress = stress_of_strain__tranisotropic__sf(strain, lam,mu, Elam,Emu,Egam,m)
    end
    
    function strain_of_stress__tranisotropic(stress, lam,mu, Elam,Emu,Egam,m) result(strain)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: stress(3,3), lam,mu, Elam,Emu,Egam,m(3)
        real(kind=dp)             :: strain(3,3)
        
        strain = strain_of_stress__tranisotropic__sf(stress, lam,mu, Elam,Emu,Egam,m)
    end
    
    subroutine Cij_to_Lame(C11,C33,C55,C12,C13, lam,mu,Elam,Emu,Egam)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in)  :: C11,C33,C55,C12,C13
        real(kind=dp), intent(out) :: lam,mu, Elam,Emu,Egam ! Lame parameters and their directional enhancements
        
        call Cij_to_Lame__traniso(C11,C33,C55,C12,C13, lam,mu,Elam,Emu,Egam)
    end
    
    !---------------------------------
    ! ELASTIC PHASE VELOCITIES
    !---------------------------------
    ! For a composite material consisting of transversely isotropic grains
    
    function Vi_elastic_tranisotropic(nlm, alpha, lam,mu,Elam,Emu,Egam, rho, theta_n,phi_n) result(Vi)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: alpha, lam,mu,Elam,Emu,Egam, rho
        real(kind=dp), intent(in)    :: theta_n(:), phi_n(:) ! arrays of theta and phi values to calculate phase velocities (vj) along
        real(kind=dp)                :: Vi(3,size(theta_n)) ! qS1, qS2, qP phase velocities

        Vi = Vi_elastic_tranisotropic__sf(nlm, alpha, lam,mu,Elam,Emu,Egam, rho, theta_n,phi_n) 
    end
  
    function Qnorm(nlm, alpha, lam,mu,Elam,Emu,Egam) result(Qn)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: alpha, lam,mu,Elam,Emu,Egam
        real(kind=dp)                :: Qn(3,3)

        Qn = Qnorm__sf(nlm, alpha, lam,mu,Elam,Emu,Egam)
    end
    
    !---------------------------------
    ! REDUCED FORM
    !---------------------------------

    subroutine get_rnlm_len(rnlm_len) 
        implicit none
        integer, intent(out) :: rnlm_len
        
        rnlm_len = rnlm_len__sf
    end

    function nlm_to_rnlm(nlm, rnlm_len) result(rnlm) 
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        integer, intent(in)          :: rnlm_len
        complex(kind=dp)             :: rnlm(rnlm_len)
        
        rnlm = nlm_to_rnlm__sf(nlm)
    end
    
    function rnlm_to_nlm(rnlm, nlm_len) result (nlm)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: rnlm(:)
        integer, intent(in)          :: nlm_len
        complex(kind=dp)             :: nlm(nlm_len)
        
        nlm = rnlm_to_nlm__sf(rnlm)
    end
    
    subroutine reduce_M(M,rnlm_len, Mrr,Mri,Mir,Mii)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in)                             :: M(:,:)
        integer, intent(in)                                      :: rnlm_len
        real(kind=dp), intent(out), dimension(rnlm_len,rnlm_len) :: Mrr,Mri,Mir,Mii
        
        call reduce_M__sf(M, Mrr,Mri,Mir,Mii)
    end
    
    !---------------------------------
    ! AUX
    !---------------------------------
    
    function rotate_nlm(nlm, theta,phi) result (nlm_rot)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:) 
        real(kind=dp), intent(in)    :: theta, phi 
        complex(kind=dp)             :: nlm_rot(size(nlm))
        
        nlm_rot = rotate_nlm__sf(nlm, theta,phi)
    end
  
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
        include "include/ddrx-coupling-weights.f90"
        g = k*g/doubleinner22(tau,tau) ! = D
        
        Davg = doubleinner22(matmul(tau,tau), a2__sf(nlm)) - doubleinner22(tau,doubleinner42(a4__sf(nlm),tau)) ! (tau.tau):a2 - tau:a4:tau 
        Davg = Davg/doubleinner22(tau,tau) ! = <D>
        
        ! Davg is the "isotropic reference to be subtracted"; that is, we subtract a constant value from "D" (on S^2), which amounts to adjusting the value of the isotropic speactral expansion coefficient (l,m=0,0). 
        ! The factor of Sqrt(4*Pi) ensures that the isotropic contribution is indeed n00*Y00 = Sqrt(4*Pi)*Davg*Y00 = Davg
        g(1) = g(1) - Sqrt(4*Pi)*Davg 
        
        ! Calculate D - <D>
        DDRX_decayrate = g 
    end
    
    function Gamma0(eps, T, A, Q)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3) ! strain-rate tensor
        real(kind=dp), intent(in) :: T, A, Q
        real(kind=dp)             :: Gamma0
        
        Gamma0 = Gamma0__sf(eps, T, A, Q)
    end
    
    function apply_bounds(nlm) result (nlm_bounded)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        complex(kind=dp)             :: nlm_bounded(size(nlm))
    
        nlm_bounded = apply_bounds__sf(nlm)
    end
    
    function Sl(nlm, l)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        integer, intent(in)          :: l
        real(kind=dp)                :: Sl
       
       Sl = Sl__sf(nlm, l) 
    end
    
    !---------------------------------
    ! EXTERNAL, NON-CORE FEATURES
    !---------------------------------

    ! Elmer ice flow model 
    include "elmer/specfabpy_elmer.f90"

    ! JOSEF ice flow model 
    include "josef/specfabpy_josef.f90"
        
end module specfabpy 

