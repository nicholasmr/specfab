! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2024

! This file is a wrapper for the Fortran routines to be provided in specfabpy

module specfabpy_const

    implicit none
    integer, parameter :: dp = 8 ! Default precision
    real(kind=dp), parameter :: Pi = 3.141592653589793d0
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
        a2__sf => a2, a4__sf => a4, a6__sf => a6, & 
        a2_to_nlm__sf => a2_to_nlm, a4_to_nlm__sf => a4_to_nlm, a6_to_nlm__sf => a6_to_nlm, & ! canonical representations
        ae2_to_a2__sf => ae2_to_a2, ae4_to_a4__sf => ae4_to_a4, & ! Elmer representations
        a4_to_mat__sf => a4_to_mat, & ! a4 in Mandel notation
        a4_eigentensors__sf => a4_eigentensors, & ! 3x3 eigentensors of a4
        a4_IBOF__sf => a4_IBOF, & ! from Elmer
        a2_orth__sf => a2_orth, a4_orth__sf => a4_orth, &

        ! Rheologies
        rheo_rev_isotropic__sf     => rheo_rev_isotropic, &
        rheo_fwd_isotropic__sf     => rheo_fwd_isotropic, &
        rheo_rev_tranisotropic__sf => rheo_rev_tranisotropic, &
        rheo_fwd_tranisotropic__sf => rheo_fwd_tranisotropic, &
        rheo_rev_orthotropic__sf   => rheo_rev_orthotropic, &
        rheo_fwd_orthotropic__sf   => rheo_fwd_orthotropic, &
        rheo_rev_orthotropic_Martin__sf => rheo_rev_orthotropic_Martin, &
        rheo_fwd_orthotropic_Martin__sf => rheo_fwd_orthotropic_Martin, &
        rheo_rev_orthotropic_Pettit__sf => rheo_rev_orthotropic_Pettit, &
        rheo_fwd_orthotropic_Pettit__sf => rheo_fwd_orthotropic_Pettit, &
        rheo_fwd_tranisotropic_sachshomo__sf  => rheo_fwd_tranisotropic_sachshomo, &
        rheo_fwd_tranisotropic_taylorhomo__sf => rheo_fwd_tranisotropic_taylorhomo, &
        
        ! Elasticity
        elas_rev_tranisotropic__sf => elas_rev_tranisotropic, &
        elas_fwd_tranisotropic__sf => elas_fwd_tranisotropic, &
        Vi_elastic_tranisotropic__sf => Vi_elastic_tranisotropic, &
        Vi_elastic_orthotropic__sf => Vi_elastic_orthotropic, &
        Vi_elastic_orthotropic__discrete__sf => Vi_elastic_orthotropic__discrete, &
        Qnorm_tranisotropic__sf => Qnorm_tranisotropic, &
        Cij_to_Lame_tranisotropic__sf => Cij_to_Lame_tranisotropic, &
        Cij_to_Lame_orthotropic__sf   => Cij_to_Lame_orthotropic, &
        
        ! Electromagnetic
        Vi_electromagnetic_tranisotropic__sf => Vi_electromagnetic_tranisotropic, &
        
        ! Fluid enhancement factors
        Eij_tranisotropic__sf => Eij_tranisotropic, &
!        Eij_tranisotropic_HSP__sf => Eij_tranisotropic_HSP, &
        Eij_orthotropic__sf   => Eij_orthotropic, &
        Evw_tranisotropic__sf => Evw_tranisotropic, &
        Evw_tranisotropic_HSP__sf => Evw_tranisotropic_HSP, &
        Evw_orthotropic__sf   => Evw_orthotropic, &
        Evw_orthotropic_discrete__sf => Evw_orthotropic_discrete, &
        frame__sf => frame, &
        E_EIE__sf => E_EIE, E_CAFFE__sf => E_CAFFE, E_ESTAR__sf => E_ESTAR, &
        
        ! nlm and rnlm representations
        lm__sf => lm, &
        nlm_len__sf => nlm_len, rnlm_len__sf => rnlm_len, &
        nlm_to_rnlm__sf => nlm_to_rnlm, rnlm_to_nlm__sf => rnlm_to_nlm, &
        reduce_M__sf => reduce_M, &
        rotate_nlm__sf => rotate_nlm, &
        rotate_vector__sf => rotate_vector, &
        Sl__sf => Sl, & ! power spectrum
        
        ! Numerics
        Ldiag__sf => Ldiag, &
        apply_bounds__sf => apply_bounds, & 

        ! Fabric processes
        Gamma0__sf => Gamma0, &

        ! Deformation modes (code in include/specfabpy_deformationmodes.f90)
        F_to_strain__sf => F_to_strain, ugrad_to_D_and_W__sf => ugrad_to_D_and_W, &
        pureshear_b__sf => pureshear_b, pureshear_strainii_to_t__sf => pureshear_strainii_to_t, &
        pureshear_F__sf => pureshear_F, pureshear_ugrad__sf => pureshear_ugrad, &
        simpleshear_gamma__sf => simpleshear_gamma, simpleshear_gamma_to_t__sf => simpleshear_gamma_to_t, &
        simpleshear_F__sf => simpleshear_F, simpleshear_ugrad__sf => simpleshear_ugrad, &

        ! AUX
        vec_to_mat_voigt__sf => vec_to_mat_voigt, &
        nhat40_empcorr_ice__sf => nhat40_empcorr_ice, &
        nlm_ideal__sf => nlm_ideal, &        
        nlm_isvalid__sf => nlm_isvalid, &
        Lame_olivine_A2X__sf => Lame_olivine_A2X, &
        L2len__sf => L2len, L4len__sf => L4len, L6len__sf => L6len, L8len__sf => L8len
        
    implicit none
    
    integer :: Lcap

    integer, parameter :: L2len = L2len__sf
    integer, parameter :: L4len = L4len__sf
    integer, parameter :: L6len = L6len__sf
    integer, parameter :: L8len = L8len__sf
    
    integer, parameter :: I20 = ILm(2)-1
    integer, parameter :: I40 = ILm(4)-1
    integer, parameter :: I60 = ILm(6)-1
    integer, parameter :: I80 = ILm(8)-1    

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

    function nlm_len()
        implicit none
        integer :: nlm_len
        nlm_len = nlm_len__sf
    end
    
    !---------------------------------
    ! FABRIC DYNAMICS
    !---------------------------------
    
    function M_LROT(nlm, eps, omg, iota, zeta) 
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: eps(3,3), omg(3,3), iota, zeta
        complex(kind=dp)             :: M_LROT(size(nlm),size(nlm))
        
        M_LROT = M_LROT__sf(eps, omg, iota, zeta)
    end
    
    function M_DDRX(nlm, tau)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: tau(3,3)
        complex(kind=dp)             :: M_DDRX(size(nlm),size(nlm))
        
        M_DDRX = M_DDRX__sf(nlm, tau)
    end

    function M_DDRX_src(nlm, tau)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in) :: tau(3,3)
        complex(kind=dp)          :: M_DDRX_src(size(nlm),size(nlm))
        
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
!        integer                      :: ii
        
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
    
    function Evw_tranisotropic(nlm, v,w,tau, Eij_grain,alpha,n_grain) result(Evw)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in) :: v(3),w(3), tau(3,3), Eij_grain(2),alpha
        integer, intent(in)       :: n_grain
        real(kind=dp)             :: Evw
        
        Evw = Evw_tranisotropic__sf(v,w, tau, nlm, Eij_grain,alpha,n_grain)
    end
    
    function Eij_tranisotropic(nlm, e1,e2,e3, Eij_grain,alpha,n_grain) result(Eij)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: e1(3),e2(3),e3(3)
        real(kind=dp), intent(in)    :: Eij_grain(2), alpha
        integer, intent(in)          :: n_grain
        real(kind=dp)                :: Eij(6)
        
        Eij = Eij_tranisotropic__sf(nlm, e1,e2,e3, Eij_grain,alpha,n_grain)
    end

    function Evw_tranisotropic_HSP(nlm, v,w,tau, Eij_grain,beta,n_grain) result(Evw)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in) :: v(3),w(3), tau(3,3), Eij_grain(2),beta
        integer, intent(in)       :: n_grain
        real(kind=dp)             :: Evw
        
        Evw = Evw_tranisotropic_HSP__sf(v,w, tau, nlm, Eij_grain,beta,n_grain)
    end

    function Evw_orthotropic(nlm_r1, nlm_r2, nlm_r3, v,w,tau, Eij_grain, alpha, n_grain) result(Evw)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
        real(kind=dp), intent(in)    :: Eij_grain(6), alpha, v(3),w(3), tau(3,3)
        integer, intent(in)          :: n_grain
        real(kind=dp)                :: Evw
        
        Evw = Evw_orthotropic__sf(v,w, tau, nlm_r1,nlm_r2,nlm_r3, Eij_grain,alpha,n_grain)
    end
    
    function Evw_orthotropic_discrete(mi, v,w,tau, Eij_grain, alpha, n_grain) result(Evw)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: mi(:,:,:) ! (3,3,N) = (m'_i, xyz comp., grain no.) 
        real(kind=dp), intent(in) :: Eij_grain(6), alpha, v(:,:),w(:,:), tau(:,:,:)
        integer, intent(in)       :: n_grain
        real(kind=dp)             :: Evw(size(v)/3)
        
        Evw = Evw_orthotropic_discrete__sf(v,w, tau, mi, Eij_grain,alpha,n_grain)
    end
    
    function Eij_orthotropic(nlm_r1, nlm_r2, nlm_r3, e1,e2,e3, Eij_grain,alpha,n_grain) result(Eij)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
        real(kind=dp), intent(in)    :: e1(3),e2(3),e3(3)
        real(kind=dp), intent(in)    :: Eij_grain(6), alpha
        integer, intent(in)          :: n_grain
        real(kind=dp)                :: Eij(6)
        
        Eij = Eij_orthotropic__sf(nlm_r1,nlm_r2,nlm_r3, e1,e2,e3, Eij_grain,alpha,n_grain)
    end
    
    ! Equivalent isotropic enhancement factors
    
    function E_EIE(eps, nglen, q, nlm, Eij_grain,alpha,n_grain) result(E)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in)    :: eps(3,3), nglen, q, Eij_grain(2), alpha
        integer, intent(in)          :: n_grain
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: E
        
        E = E_EIE__sf(eps, nglen, q, nlm, Eij_grain,alpha,n_grain)
    end
    
    function E_CAFFE(eps, nlm, Emin, Emax)  result(E)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in)    :: eps(3,3), Emin, Emax
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: E
        
        E = E_CAFFE__sf(eps, nlm, Emin, Emax)
    end
    
!    function Eiso_ESTAR(eps, ...)  result(E)
!        use specfabpy_const
!        implicit none
!        real(kind=dp), intent(in)    :: eps(3,3), Emin, Emax
!        complex(kind=dp), intent(in) :: nlm(:)
!        real(kind=dp)                :: E
!        
!        E = Eiso_ESTAR__sf(eps, ...)
!    end
    
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
    
    function a6(nlm) 
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a6(3,3,3,3,3,3)

        a6 = a6__sf(nlm)
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

    function a6_to_nlm(a6) result(nlm)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: a6(3,3,3,3,3,3)
        integer, parameter        :: nlmlen = 1+5+9+13
        complex(kind=dp)          :: nlm(nlmlen)
        
        nlm = a6_to_nlm__sf(a6)
    end
    
    function a4_to_mat(A) result (M)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: A(3,3,3,3)
        real(kind=dp)             :: M(6,6)
        
        M = a4_to_mat__sf(A)
    end
        
    function a2_orth(blm,nlm) result(a2)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: a2(3,3)
        
        a2 = a2_orth__sf(blm,nlm)
    end
  
    function a4_orth(blm,nlm) result(a4)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: a4(3,3,3,3)
        
        a4 = a4_orth__sf(blm,nlm)
    end
        
    !---------------------------------
    ! ISOTROPIC RHEOLOGY 
    !---------------------------------

    function rheo_fwd_isotropic(tau, A,n) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, n
        real(kind=dp)             :: eps(3,3)
        
        eps = rheo_fwd_isotropic__sf(tau, A,n)
    end
    
    function rheo_rev_isotropic(eps, A,n) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, n
        real(kind=dp)             :: tau(3,3)
        
        tau = rheo_rev_isotropic__sf(eps, A,n)
    end
        
    !---------------------------------
    ! TRANSVERSELY ISOTROPIC RHEOLOGY 
    !---------------------------------

    function rheo_fwd_tranisotropic(tau, A,n, m,Eij) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A,n, m(3),Eij(2)
        real(kind=dp)             :: eps(3,3)
        
        eps = rheo_fwd_tranisotropic__sf(tau, A,n, m,Eij)
    end
    
    function rheo_rev_tranisotropic(eps, A,n, m,Eij) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A,n, m(3),Eij(2)
        real(kind=dp)             :: tau(3,3)
        
        tau = rheo_rev_tranisotropic__sf(eps, A,n, m,Eij)
    end
    
!    function rheo_fwd_tranisotropic_sachshomo(tau, nlm, Eij_grain,n_grain) result(eps)
!        use specfabpy_const
!        implicit none
!        complex(kind=dp), intent(in) :: nlm(:)
!        real(kind=dp), intent(in)    :: Eij_grain(2), tau(3,3)
!        integer, intent(in)          :: n_grain
!        real(kind=dp)                :: eps(3,3)
!        
!        eps = rheo_fwd_tranisotropic_sachshomo__sf(tau, nlm, Eij_grain,n_grain)
!    end
    
    !---------------------------------
    ! ORTHOTROPIC RHEOLOGY 
    !---------------------------------
    
    function rheo_fwd_orthotropic(tau, A,n, m1,m2,m3, Eij) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A,n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: eps(3,3)
        
        eps = rheo_fwd_orthotropic__sf(tau, A,n, m1,m2,m3, Eij)
    end
    
    function rheo_rev_orthotropic(eps, A,n, m1,m2,m3, Eij) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A,n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: tau(3,3)
        
        tau = rheo_rev_orthotropic__sf(eps, A,n, m1,m2,m3, Eij)
    end
    
    ! Pettit's hypothesis that the fluidity is orientation independent.
    function rheo_fwd_orthotropic_Pettit(tau, A,n, m1,m2,m3, Eij) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A,n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: eps(3,3)
        
        eps = rheo_fwd_orthotropic_Pettit__sf(tau, A,n, m1,m2,m3, Eij)
    end
    
    function rheo_rev_orthotropic_Pettit(eps, A,n, m1,m2,m3, Eij) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A,n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: tau(3,3)
        
        tau = rheo_rev_orthotropic_Pettit__sf(eps, A,n, m1,m2,m3, Eij)
    end
    
    ! Martin's hypothesis that the viscosity is orientation independent.
    function rheo_fwd_orthotropic_Martin(tau, A,n, m1,m2,m3, Eij) result(eps)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A,n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: eps(3,3)
        
        eps = rheo_fwd_orthotropic_Martin__sf(tau, A,n, m1,m2,m3, Eij)
    end
    
    function rheo_rev_orthotropic_Martin(eps, A,n, m1,m2,m3, Eij) result(tau)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A,n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: tau(3,3)
        
        tau = rheo_rev_orthotropic_Martin__sf(eps, A,n, m1,m2,m3, Eij)
    end
    
    !---------------------------------
    ! TRANSVERSELY ISOTROPIC ELASTICITY 
    !---------------------------------
    
    function elas_rev_tranisotropic(strain, lam,mu, Elam,Emu,Egam,m) result(stress)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: strain(3,3), lam,mu, Elam,Emu,Egam,m(3)
        real(kind=dp)             :: stress(3,3)
        
        stress = elas_rev_tranisotropic__sf(strain, lam,mu, Elam,Emu,Egam,m)
    end
    
    function elas_fwd_tranisotropic(stress, lam,mu, Elam,Emu,Egam,m) result(strain)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: stress(3,3), lam,mu, Elam,Emu,Egam,m(3)
        real(kind=dp)             :: strain(3,3)
        
        strain = elas_fwd_tranisotropic__sf(stress, lam,mu, Elam,Emu,Egam,m)
    end
    
    subroutine Cij_to_Lame_tranisotropic(Cij, Lame)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in)  :: Cij(5)
        real(kind=dp), intent(out) :: Lame(5)
        
        call Cij_to_Lame_tranisotropic__sf(Cij, Lame)
    end
    
    !---------------------------------
    ! ORTHOTROPIC ELASTICITY 
    !---------------------------------
    
    subroutine Cij_to_Lame_orthotropic(Cij, Lame)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in)  :: Cij(9)
        real(kind=dp), intent(out) :: Lame(9) 
        
        call Cij_to_Lame_orthotropic__sf(Cij, Lame)
    end
    
    !---------------------------------
    ! ELASTIC PHASE VELOCITIES
    !---------------------------------
    
    ! For a composite material consisting of orthotropic grains
    
    function Vi_elastic_orthotropic(nlm_1,nlm_2,nlm_3, alpha,lame_grain,rho, theta,phi) result(Vi)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm_1(:), nlm_2(:), nlm_3(:)
        real(kind=dp), intent(in)    :: alpha, lame_grain(9), rho
        real(kind=dp), intent(in)    :: theta(:), phi(:) ! arrays of theta and phi values to calculate phase velocities along
        real(kind=dp)                :: Vi(3,size(theta)) ! qS1, qS2, qP phase velocities

        Vi = Vi_elastic_orthotropic__sf(nlm_1,nlm_2,nlm_3, alpha,lame_grain,rho, theta,phi) 
    end
    
    function Vi_elastic_orthotropic__discrete(mi, alpha,lame_grain,rho, theta,phi) result(Vi)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: mi(:,:,:)
        real(kind=dp), intent(in) :: alpha, lame_grain(9), rho
        real(kind=dp), intent(in) :: theta(:), phi(:) ! arrays of theta and phi values to calculate phase velocities along
        real(kind=dp)             :: Vi(3,size(theta)) ! qS1, qS2, qP phase velocities

        Vi = Vi_elastic_orthotropic__discrete__sf(mi, alpha,lame_grain,rho, theta,phi) 
    end
    
    ! For a composite material consisting of transversely isotropic grains
    
    function Vi_elastic_tranisotropic(nlm, alpha, lame_grain, rho, theta_n,phi_n) result(Vi)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: alpha, lame_grain(5), rho
        real(kind=dp), intent(in)    :: theta_n(:), phi_n(:) ! arrays of theta and phi values to calculate phase velocities along
        real(kind=dp)                :: Vi(3,size(theta_n)) ! qS1, qS2, qP phase velocities

        Vi = Vi_elastic_tranisotropic__sf(nlm, alpha, lame_grain, rho, theta_n,phi_n) 
    end
  
    function Qnorm_tranisotropic(nlm, alpha, lame_grain) result(Qnorm)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: alpha, lame_grain(5)
        real(kind=dp)                :: Qnorm(3,3)

        Qnorm = Qnorm_tranisotropic__sf(nlm, alpha, lame_grain)
    end
    
    !---------------------------------
    ! ELECTROMAGNETIC PHASE VELOCITIES
    !---------------------------------
    
    function Vi_electromagnetic_tranisotropic(nlm, eps_m, eps_t, mu, theta_n,phi_n) result(Vi)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: eps_m, eps_t, mu
        real(kind=dp), intent(in)    :: theta_n(:), phi_n(:) ! arrays of theta and phi values to calculate phase velocities along
        real(kind=dp)                :: Vi(2,size(theta_n)) ! qS1, qS2 phase velocities

        Vi = Vi_electromagnetic_tranisotropic__sf(nlm, eps_m, eps_t, mu, theta_n,phi_n)
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
    
    function rotate_vector(v, theta, phi) result(w)
        use specfabpy_const    
        implicit none
        real(kind=dp), intent(in) :: v(3), theta, phi
        real(kind=dp) :: w(3)
        
        w = rotate_vector__sf(v, theta, phi)
    end
    
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
        complex(kind=dp) :: g(lmDDRX_len), qt(-2:2)
        real(kind=dp)    :: evD, k

        qt = quad_rr(tau) ! Quadric expansion coefficients
        include "include/ddrx-coupling-weights.f90" ! Harmonic expansion coefficients of D (requires qt, sets g and k)
        g = k*g * 5/doubleinner22(tau,tau) ! Normalize by tau:tau/5

        evD = ev_D(nlm, tau) ! <D>
        
        ! <D> is the isotropic reference level to be subtracted from D (on S^2), amounting to adjusting the isotropic harmonic mode (l,m=0,0).
        g(1) = g(1) - Sqrt(4*Pi)*evD ! Prefactor Sqrt(4*Pi) ensures the isotropic contribution is indeed n00*Y00 = Sqrt(4*Pi)*<D>*Y00 = <D>
        
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
    
    function vec_to_mat_voigt(v) result (M)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: v(6) 
        real(kind=dp)             :: M(3,3) 
        
        M = vec_to_mat_voigt__sf(v)
    end
    
    function nhat40_empcorr_ice(nhat20) result(nhat40)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: nhat20(:)
        real(kind=dp)             :: nhat40(size(nhat20))

        nhat40 = nhat40_empcorr_ice__sf(nhat20)
    end
    
    function nlm_ideal(m, colat, L) result(nlm)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in) :: m(3), colat 
        integer, intent(in)       :: L
        complex(kind=dp)          :: nlm((L+1)*(L+2)/2)
        nlm(:) = nlm_ideal__sf(m, colat, L)
    end
    
    function nlm_isvalid(nhat20, nhat40) result(isvalid)
        use specfabpy_const
        implicit none
        complex(kind=dp), intent(in) :: nhat20(:), nhat40(:) 
        logical                      :: isvalid(size(nhat20))
        isvalid = nlm_isvalid__sf(nhat20, nhat40)
    end 
    
    subroutine Lame_olivine_A2X(LameA, Xtype, LameX)
        use specfabpy_const
        implicit none
        real(kind=dp), intent(in)  :: LameA(9) ! (lam11,lam22,lam33, lam23,lam13,lam12, mu1,mu2,mu3) w.r.t. m1',m2',m3'
        real(kind=dp), intent(out) :: LameX(9)
        character*1, intent(in) :: Xtype
        call Lame_olivine_A2X__sf(LameA, Xtype, LameX)
    end
    
!    ! DELETE ONCE DEBUG FINISHED:
!    subroutine orthstructs(nlm_1,nlm_2,nlm_3, mi,  a2v_nlm, a2v_mi, a4v_nlm, a4v_mi)
!        use specfabpy_const
!        implicit none
!        complex(kind=dp), intent(in) :: nlm_1(:), nlm_2(:), nlm_3(:)
!        real(kind=dp), intent(in)    :: mi(:,:,:) ! (3,3,N) = (m'_i, xyz comp., grain no.) 
!        real(kind=dp)                :: a4v_jk_sym2(3,6,6), a4v_jk_sym4(3,6,6) ! see aiv_orthotropic() for defs
!        real(kind=dp), intent(out)   :: a2v_nlm(3,6), a2v_mi(3,6), a4v_nlm(3,6,6), a4v_mi(3,6,6)
!            
!        call aiv_orthotropic(nlm_1,nlm_2,nlm_3, a2v_nlm,a4v_nlm, a4v_jk_sym2,a4v_jk_sym4)
!        call aiv_orthotropic_discrete(mi, a2v_mi,a4v_mi, a4v_jk_sym2,a4v_jk_sym4)
!    end
    
    !---------------------------------
    ! EXTERNAL, NON-CORE FEATURES
    !---------------------------------

    ! Elmer ice flow model 
    include "elmer/specfabpy_elmer.f90"

    ! JOSEF ice flow model 
    include "josef/specfabpy_josef.f90"
        
    ! Deformation modes module
    include "include/specfabpy_deformationmodes.f90"
        
end module specfabpy 

