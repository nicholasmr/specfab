! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2023

! Bulk polycrystalline homogenization schemes for transversely isotropic and ortotropic grains, both viscously and elastically.

module homogenizations  

    use tensorproducts
    use mandel
    use moments
    use rheologies
    use elasticities

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
!    integer, parameter :: identity(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    real(kind=dp), parameter :: identity(3,3) = reshape([1.0d0,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,0.0d0,1.0d0], [3,3])
    
    ! Aux
    real, parameter, private :: Pi = 3.141592653589793
    integer, private       :: Lcap
    complex(kind=dp)       :: nlm_iso(1+5+9+13+17) ! nlm for L=8 is needed
    real(kind=dp), private :: ev_c2_iso(3,3), ev_c4_iso(3,3,3,3), ev_c6_iso(3,3,3,3, 3,3), ev_c8_iso(3,3,3,3, 3,3,3,3) ! <c^k> for isotropic n(theta,phi)
    real(kind=dp), private :: ev_c4_iso_Mandel(6,6) 
        
    real(kind=dp), private, parameter  :: identity_vec6(6) =  [1,1,1, 0,0,0] ! = mat_to_vec(identity)
    real(kind=dp), private, parameter  :: identity6(6,6) = reshape([1, 0, 0, 0, 0, 0, &
                                                                    0, 1, 0, 0, 0, 0, & 
                                                                    0, 0, 1, 0, 0, 0, & 
                                                                    0, 0, 0, 1, 0, 0, & 
                                                                    0, 0, 0, 0, 1, 0, & 
                                                                    0, 0, 0, 0, 0, 1], [6,6])
        
contains      

    !---------------------------------
    ! INIT
    !---------------------------------
       
    subroutine inithomogenizations(Lcap_)

        ! Needs to be called once before using the module routines.

        implicit none    
        integer, intent(in) :: Lcap_ ! Truncation "Lcap"
                
        Lcap = Lcap_ ! Save internal copy

        ! Calculate structure tensors for an isotropic fabric
        nlm_iso(:) = 0.0d0
        nlm_iso(1) = 1/Sqrt(4*Pi) ! Normalized ODF 
        call f_ev_ck(nlm_iso, 'f', ev_c2_iso,ev_c4_iso,ev_c6_iso,ev_c8_iso) ! a^(k) (k=2,4,6,8) for an isotropic fabric
        ev_c4_iso_Mandel = a4_to_mat(ev_c4_iso)
    end

    !---------------------------------
    ! VISCOUS -- Transversely isotropic grains
    !---------------------------------

    function rheo_fwd_tranisotropic_sachshomo(tau, nlm, Eij_grain,n_grain) result(eps)
    
        ! Grain-averaged forward rheology subject to Sachs homogenization (constant stress over all grains).
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: Eij_grain(2), tau(3,3)
        integer, intent(in)          :: n_grain
        
        real(kind=dp)             :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
        real(kind=dp)             :: ev_etac0 = 0, ev_etac2(3,3), ev_etac4(3,3,3,3) ! <eta*c^k> terms
        real(kind=dp)             :: eps(3,3), coefA,coefB,coefC, tausq(3,3), I2
        real(kind=dp)             :: a2v(6), a4v(6,6), a4tau(3,3)

        call rheo_params_tranisotropic(Eij_grain,3,DFLOAT(n_grain),1, coefA,coefB,coefC)
        I2 = doubleinner22(tau,tau)

        if (n_grain .eq. 1) then
            call f_ev_ck_Mandel(nlm, a2v, a4v) ! Structure tensors in Mandel notation: a2v = ev_c2_Mandel, a4v = ev_c4_Mandel
            ev_etac0 = 1.0d0
            ev_etac2 = vec_to_mat(a2v) ! = a2 = ev_c2
            a4tau    = vec_to_mat(matmul(a4v,mat_to_vec(tau))) ! = a4 : tau = doubleinner42(ev_etac4,tau)

        else if (n_grain .eq. 3) then
            if (size(nlm) < (8 +1)*(8 +2)/2) then ! L < 8 ?
                stop "specfab error: Sachs homogenization with n'=3 requires L >= 8."
            end if
            call f_ev_ck(nlm, 'f', ev_c2,ev_c4,ev_c6,ev_c8) ! Calculate structure tensors of orders 2,4,6,8
            tausq = matmul(tau,tau)
            ev_etac0 = I2*  1.0 + coefB*doubleinner22(doubleinner42(ev_c4,tau),tau) + 2*coefC*doubleinner22(ev_c2,tausq)
            ev_etac2 = I2*ev_c2 + coefB*doubleinner42(doubleinner62(ev_c6,tau),tau) + 2*coefC*doubleinner42(ev_c4,tausq)
            ev_etac4 = I2*ev_c4 + coefB*doubleinner62(doubleinner82(ev_c8,tau),tau) + 2*coefC*doubleinner62(ev_c6,tausq)
            a4tau = doubleinner42(ev_etac4,tau)
            
        ! Same as n'=3 but disregarding the orientation-dependent terms in the nonlinear grain fluidity.        
        else if (n_grain .eq. -3) then 
            call f_ev_ck_Mandel(nlm, a2v, a4v) ! Structure tensors in Mandel notation: a2v = ev_c2_Mandel, a4v = ev_c4_Mandel
            ev_etac0 = I2 * 1.0d0
            ev_etac2 = I2 * vec_to_mat(a2v) ! = I2*a2 = I2*ev_c2
            a4tau    = I2 * vec_to_mat(matmul(a4v,mat_to_vec(tau))) ! = I2 * a4 : tau = I2 * doubleinner42(ev_etac4,tau)
        else
            print *, "n'=", n_grain
            stop "specfab error: unsupported n'"
        end if

        eps = ev_etac0*tau - coefA*doubleinner22(ev_etac2,tau)*identity + coefB*a4tau + coefC*(matmul(tau,ev_etac2)+matmul(ev_etac2,tau))
    end
    
    function rheo_fwd_tranisotropic_taylorhomo(tau, nlm, Eij_grain,n_grain) result(eps)
        
        ! Grain-averaged forward rheology subject to Taylor homogenization (constant strainrate over all grains).
        ! NOTE: Taylor model supports only n'=1. n_grain is required anyway for future compatibility with n'>1.
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: Eij_grain(2), tau(3,3)
        integer, intent(in)          :: n_grain ! Dummy variable: <eps'(tau)> with Taylor hypothesis is implemented only for n' = 1.
        
        real(kind=dp), parameter  :: s = sqrt(2.0)
        real(kind=dp)             :: a2mat(3,3), a2v(6), a4v(6,6)
        real(kind=dp)             :: P(6,6), L(6,6), tau_vec(6,1),  P_reg(6,6),tau_vec_reg(6,1)
        real(kind=dp)             :: eps(3,3), coefA,coefB,coefC
        integer                   :: info
    !    integer :: ipiv(9), work
                                                               
        call rheo_params_tranisotropic(Eij_grain,3,DFLOAT(n_grain),-1, coefA,coefB,coefC)
        
        call f_ev_ck_Mandel(nlm, a2v, a4v) ! Structure tensors in Mandel notation: a2v = ev_c2_Mandel, B = ev_c4_Mandel
        a2mat = vec_to_mat(a2v) ! = a2
        L = anticommutator_Mandel(a2mat)
        
        ! Matrix "P" of the vectorized bulk Taylor rheology: vec(tau) = matmul(P, vec(eps))
        P = identity6 - coefA*outerprod6(identity_vec6,a2v) + coefB*a4v + coefC*L 
        
        ! Solve inverse problem; we are seeking eps given tau.
        tau_vec(:,1) = mat_to_vec(tau)
        call dposv('L', 6, 1, P, 6, tau_vec, 6, info) ! tau_vec is now "eps_vec" solution. For some reason, "U" does not work.

        ! Ill posed? => Regularization needed.
        if (info /= 0) then
            P_reg       = matmul(TRANSPOSE(P),P) + 1e-6*identity6
            tau_vec_reg = matmul(TRANSPOSE(P),tau_vec)
            call dposv('L', 6, 1, P_reg, 6, tau_vec_reg, 6, info)
            tau_vec = tau_vec_reg
            if (info /= 0) then
                stop 'specfab error: Taylor viscosity-matrix inversion failed! Please check the ODF is correct (reducing the fabric integration time-step, and/or increasing regularization, for transient problems can often help).'
            end if
        end if
        
        ! Revert solution to 3x3
        eps = vec_to_mat(tau_vec)
    end

    function rheo_fwd_tranisotropic_sachshomo__isotropic(tau, Eij_grain,n_grain) result(eps)
       
        ! Same as rheo_fwd_tranisotropic_sachshomo() but uses hardcoded isotropic ODF
       
        implicit none
        
        real(kind=dp), intent(in) :: Eij_grain(2), tau(3,3)
        integer, intent(in)       :: n_grain
        real(kind=dp)             :: ev_etac0 = 0.0d0, ev_etac2(3,3), ev_etac4(3,3,3,3) ! <eta*c^k> terms
        real(kind=dp)             :: eps(3,3), coefA,coefB,coefC, tausq(3,3), I2, a4tau(3,3)

        call rheo_params_tranisotropic(Eij_grain,3,DFLOAT(n_grain),1, coefA,coefB,coefC)
        I2 = doubleinner22(tau,tau)
                    
        if (n_grain .eq.  1) then 
            eps = (1 + 2.0d0/15*coefB + 2.0d0/3*coefC)*tau ! n'=1
        else if (n_grain .eq. -3) then 
            eps = I2*eps ! Same as n'=3 but disregarding the orientation-dependent terms in the nonlinear grain fluidity.
        else if (n_grain .eq.  3) then 
            tausq = matmul(tau,tau)
            ev_etac0 = I2*1.0       + coefB*doubleinner22(doubleinner42(ev_c4_iso,tau),tau) + 2*coefC*doubleinner22(ev_c2_iso,tausq)
            ev_etac2 = I2*ev_c2_iso + coefB*doubleinner42(doubleinner62(ev_c6_iso,tau),tau) + 2*coefC*doubleinner42(ev_c4_iso,tausq)
            ev_etac4 = I2*ev_c4_iso + coefB*doubleinner62(doubleinner82(ev_c8_iso,tau),tau) + 2*coefC*doubleinner62(ev_c6_iso,tausq)
            a4tau = doubleinner42(ev_etac4,tau)
            eps = ev_etac0*tau - coefA*doubleinner22(ev_etac2,tau)*identity + coefB*a4tau + coefC*(matmul(tau,ev_etac2)+matmul(ev_etac2,tau))
        else
            print *, "n'=", n_grain
            stop "specfab error: unsupported n'"
        end if
    end
    
    function rheo_fwd_tranisotropic_taylorhomo__isotropic(tau, Eij_grain,n_grain) result(eps)

        ! Same as rheo_fwd_tranisotropic_taylorhomo() but uses hardcoded isotropic ODF

        implicit none
        
        real(kind=dp), intent(in) :: Eij_grain(2), tau(3,3)
        integer, intent(in)       :: n_grain ! Dummy variable, only n' = 1 implemented.
        real(kind=dp)             :: eps(3,3), coefA,coefB,coefC

        call rheo_params_tranisotropic(Eij_grain,3,DFLOAT(n_grain),-1, coefA,coefB,coefC)
        eps = tau/(1 + 2.0d0/15*coefB + 2.0d0/3*coefC) ! Analytical inverse when isotropic.
    end
    
    function anticommutator_Mandel(a2mat) result(L)

        ! Mandel form of anti-commutator: {a2,tau} := a2 . tau + tau . a2 

        ! I.e. the 6x6 tensor L = MandelForm[I \otimes a2 + a2 \otimes I]
        ! so that {a2,tau} = MandelVecToMat[L . MandelVector[tau]]
    
        implicit none
        
        real(kind=dp), intent(in) :: a2mat(3,3)
        real(kind=dp)             :: L(6,6)
        real(kind=dp), parameter  :: s = sqrt(2.0d0)
                                                
        ! L in Mandel notation
        L(1,:) = [2*a2mat(1,1), 0.0d0, 0.0d0, 0.0d0, s*a2mat(1,3), s*a2mat(1,2)]
        L(2,:) = [0.0d0, 2*a2mat(2,2), 0.0d0, s*a2mat(2,3), 0.0d0, s*a2mat(1,2)]
        L(3,:) = [0.0d0, 0.0d0, 2*a2mat(3,3), s*a2mat(2,3), s*a2mat(1,3), 0.0d0]
        L(4,:) = [0.0d0, s*a2mat(2,3), s*a2mat(2,3), a2mat(2,2)+a2mat(3,3), a2mat(1,2), a2mat(1,3)]
        L(5,:) = [s*a2mat(1,3), 0.0d0, s*a2mat(1,3), a2mat(1,2), a2mat(1,1)+a2mat(3,3), a2mat(2,3)]
        L(6,:) = [s*a2mat(1,2), s*a2mat(1,2), 0.0d0, a2mat(1,3), a2mat(2,3), a2mat(1,1)+a2mat(2,2)]
    end    

    !---------------------------------
    ! VISCOUS -- Orthotropic grains
    !---------------------------------
    
    ! *** EXPERIMENTAL / STILL BEING DEVELOPED ***

    function rheo_fwd_orthotropic_sachshomo(tau, nlm_r1, nlm_r2, nlm_r3, Eij_grain, n_grain) result(eps)
    
        ! Grain-averaged forward orthotropic rheology subject to Sachs homogenization (constant stress over all grains).
        
        implicit none

        real(kind=dp), intent(in)    :: tau(3,3), Eij_grain(6)
        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
        integer, intent(in)          :: n_grain
        
        real(kind=dp) :: eps(3,3), lami(6), ci(6), gam
!        real(kind=dp) :: a2v_r1(6),a2v_r2(6),a2v_r3(6), a4v_r1(6,6),a4v_r2(6,6),a4v_r3(6,6), a4v_r1r2(6,6),a4v_r1r3(6,6),a4v_r2r3(6,6)
        real(kind=dp) :: a2v_r1(6),a2v_r2(6),a2v_r3(6), a4v_r1(6,6),a4v_r2(6,6),a4v_r3(6,6), a4_r1r2(3,3,3,3),a4_r1r3(3,3,3,3),a4_r2r3(3,3,3,3)  
        real(kind=dp) :: tauv(6), Q(6,6)
                
!        real(kind=dp) :: a2v_r1(6), a4v_r1(6,6), a2v_r2(6), a4v_r2(6,6), a2v_r3(6), a4v_r3(6,6) 
!        real(kind=dp) :: a4_r1r2(3,3,3,3), a4_r1r3(3,3,3,3), a4_r2r3(3,3,3,3)

!        real(kind=dp) :: ev_r1(3,3), ev_r2(3,3), ev_r3(3,3)   
!        real(kind=dp) :: ev_r1r2_sym(3,3), ev_r1r3_sym(3,3), ev_r2r3_sym(3,3)

        call rheo_params_orthotropic(Eij_grain, DFLOAT(n_grain), lami, gam)
        ci(1:3) = 4.0d0/3 * lami(1:3)
        ci(4:6) =       2 * lami(4:6)

        tauv = mat_to_vec(tau)
            
        ! Linear grain fluidity    
        if (n_grain .eq. 1) then

!            call f_ev_ck_Mandel(nlm_r1, a2v_r1, a4v_r1) ! Structure tensors in Mandel (vectorized) notation
!            call f_ev_ck_Mandel(nlm_r2, a2v_r2, a4v_r2) 
!            ev_r1 = vec_to_mat(matmul(a4v_r1,tauv)) ! <r1^4> : tau
!            ev_r2 = vec_to_mat(matmul(a4v_r2,tauv)) ! <r2^4> : tau
!            a4_r1r2 = a4_joint(nlm_r1,nlm_r2)       ! <r1^2 r2^2>
!            
!            ! Is nlm_r3 given? Then use it, else derive it from nlm_r1 and nlm_r2            
!            if (real(nlm_r3(1)) > 1e-10) then
!                call f_ev_ck_Mandel(nlm_r3, a2v_r3, a4v_r3) 
!                ev_r3 = vec_to_mat(matmul(a4v_r3,tauv)) ! <r3^4> : tau 
!                a4_r1r3 = a4_joint(nlm_r1,nlm_r3)       ! <r1^2 r3^2>
!                a4_r2r3 = a4_joint(nlm_r2,nlm_r3)       ! <r2^2 r3^2>
!            else
!                ev_r3 = doubleinner42(a4_orth(nlm_r1, nlm_r2), tau) ! <r3^4> : tau 
!                a4_r1r3 = a4_jointcross(nlm_r1,nlm_r2) ! <r1^2 r3^2> where r3 = r1 x r2 
!                a4_r2r3 = a4_jointcross(nlm_r2,nlm_r1) ! <r2^2 r3^2> where r3 = r1 x r2 
!            end if

!            ! Symmetric double inner products
!            ev_r1r2_sym = doubleinner24(tau, a4_r1r2) + doubleinner42(a4_r1r2, tau) ! tau : <r1^2 r2^2> + <r2^2 r1^2> : tau
!            ev_r1r3_sym = doubleinner24(tau, a4_r1r3) + doubleinner42(a4_r1r3, tau) ! tau : <r1^2 r3^2> + <r3^2 r1^2> : tau
!            ev_r2r3_sym = doubleinner24(tau, a4_r2r3) + doubleinner42(a4_r2r3, tau) ! tau : <r2^2 r3^2> + <r3^2 r2^2> : tau

!            ! Finally, the strain-rate is:
!            eps = + lam1 * (ev_r2 + ev_r3 - ev_r2r3_sym)/4.0d0 &
!                  + lam2 * (ev_r1 + ev_r3 - ev_r1r3_sym)/4.0d0 &
!                  + lam3 * (ev_r1 + ev_r2 - ev_r1r2_sym)/4.0d0 &
!                  + lam4 * vec_to_mat(matmul(a4_to_mat(a4_symmetrize_orthotropic(a4_r2r3)),tauv)) &
!                  + lam5 * vec_to_mat(matmul(a4_to_mat(a4_symmetrize_orthotropic(a4_r1r3)),tauv)) &
!                  + lam6 * vec_to_mat(matmul(a4_to_mat(a4_symmetrize_orthotropic(a4_r1r2)),tauv))
!!                  + lam4 * (4*doubleinner42_firstlast_symmetric(a4_r2r3, tau))/4.0d0 &
!!                  + lam5 * (4*doubleinner42_firstlast_symmetric(a4_r1r3, tau))/4.0d0 &
!!                  + lam6 * (4*doubleinner42_firstlast_symmetric(a4_r1r2, tau))/4.0d0 

            call structensors_for_homogenized_orthotropic_rheology(nlm_r1,nlm_r2,nlm_r3, a2v_r1,a2v_r2,a2v_r3, a4v_r1,a4v_r2,a4v_r3, a4_r2r3,a4_r1r3,a4_r1r2)
            
            ! Matrix "Q" of the vectorized bulk Sachs rheology: vec(eps) = matmul(Q, vec(tau))
            Q = + ci(1) * (a4v_r2 + a4v_r3 - 2*a4_to_mat(a4_sym(a4_r2r3)))/4.0d0 &
                + ci(2) * (a4v_r1 + a4v_r3 - 2*a4_to_mat(a4_sym(a4_r1r3)))/4.0d0 &
                + ci(3) * (a4v_r1 + a4v_r2 - 2*a4_to_mat(a4_sym(a4_r1r2)))/4.0d0 &
                + ci(4) * a4_to_mat(a4_altsym(a4_r2r3)) &
                + ci(5) * a4_to_mat(a4_altsym(a4_r1r3)) &
                + ci(6) * a4_to_mat(a4_altsym(a4_r1r2))
                
            eps = vec_to_mat(matmul(Q,tauv))
        else
            eps(:,:) = 0.0d0 ! n_grain not supported, silently return 0
        end if
    end

    function rheo_fwd_orthotropic_taylorhomo(tau, nlm_r1, nlm_r2, nlm_r3, Eij_grain, n_grain) result(eps)
        
        ! Grain-averaged forward orthotropic rheology subject to Taylor homogenization (constant strainrate over all grains).
        
        implicit none
        
        real(kind=dp), intent(in)    :: tau(3,3), Eij_grain(6)
        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
        integer, intent(in)          :: n_grain

        real(kind=dp) :: eps(3,3), lami(6), ci(6), gam
        real(kind=dp) :: a2v_r1(6),a2v_r2(6),a2v_r3(6), a4v_r1(6,6),a4v_r2(6,6),a4v_r3(6,6), a4_r1r2(3,3,3,3),a4_r1r3(3,3,3,3),a4_r2r3(3,3,3,3)  
 
        real(kind=dp) :: P(6,6), tau_vec(6,1), P_reg(6,6), tau_vec_reg(6,1)
        integer       :: info
    !    integer :: ipiv(9), work
                                                               
        call rheo_params_orthotropic(Eij_grain, DFLOAT(n_grain), lami, gam)
        ci(1:3) = 4.0d0/3 * lami(1:3)/gam
        ci(4:6) =       2 * 1/lami(4:6)

        ! Linear grain viscosity    
        if (n_grain .eq. 1) then

            call structensors_for_homogenized_orthotropic_rheology(nlm_r1,nlm_r2,nlm_r3, a2v_r1,a2v_r2,a2v_r3, a4v_r1,a4v_r2,a4v_r3, a4_r2r3,a4_r1r3,a4_r1r2)
            
            ! Matrix "P" of the vectorized bulk Taylor rheology: tau = matmul(P, eps)
            P = + ci(1) * (-3.0d0/2) * (outerprod6(identity_vec6,a2v_r1) - 3*a4v_r1)/2 &
                + ci(2) * (-3.0d0/2) * (outerprod6(identity_vec6,a2v_r2) - 3*a4v_r2)/2 &
                + ci(3) * (-3.0d0/2) * (outerprod6(identity_vec6,a2v_r3) - 3*a4v_r3)/2 &
                + ci(4) * a4_to_mat(a4_altsym(a4_r2r3)) &
                + ci(5) * a4_to_mat(a4_altsym(a4_r1r3)) &
                + ci(6) * a4_to_mat(a4_altsym(a4_r1r2))
            
            ! Solve inverse problem; we are seeking eps given tau in tau = matmul(P, eps).
            tau_vec(:,1) = mat_to_vec(tau)
            call dposv('L', 6, 1, P, 6, tau_vec, 6, info) ! tau_vec is now "eps_vec" solution. For some reason, "U" does not work.

            ! Ill posed? => Regularization needed.
            if (info /= 0) then
                P_reg       = matmul(TRANSPOSE(P),P) + 1e-6*identity6
                tau_vec_reg = matmul(TRANSPOSE(P),tau_vec)
                call dposv('L', 6, 1, P_reg, 6, tau_vec_reg, 6, info)
                tau_vec = tau_vec_reg
                if (info /= 0) then
                    stop 'specfab error: Taylor viscosity-matrix inversion failed! Please check the CPO is correct (reducing the fabric integration time-step, and/or increasing regularization, for transient problems can often help).'
                end if
            end if
            
            ! Revert solution to 3x3
            eps = vec_to_mat(tau_vec)
        
        else
            eps(:,:) = 0.0d0 ! n_grain not supported, silently return 0
        end if
    end

    subroutine structensors_for_homogenized_orthotropic_rheology(nlm_r1,nlm_r2,nlm_r3, a2v_r1,a2v_r2,a2v_r3, a4v_r1,a4v_r2,a4v_r3, a4_r2r3,a4_r1r3,a4_r1r2)
    
        ! Structure tensors in Mandel (vectorized) notation used by homogenized orthtropic rheology
    
        implicit none
        
        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
!        real(kind=dp), intent(out)   :: a2v_r1(6),a2v_r2(6),a2v_r3(6), a4v_r1(6,6),a4v_r2(6,6),a4v_r3(6,6), a4v_r2r3(6,6),a4v_r1r3(6,6),a4v_r1r2(6,6) 

        real(kind=dp), intent(out) :: a2v_r1(6),a2v_r2(6),a2v_r3(6), a4v_r1(6,6),a4v_r2(6,6),a4v_r3(6,6)
        real(kind=dp), dimension(3,3,3,3), intent(out) :: a4_r2r3,a4_r1r3,a4_r1r2 
    
        call f_ev_ck_Mandel(nlm_r1, a2v_r1, a4v_r1)   ! <r1^2> and <r1^4>
        call f_ev_ck_Mandel(nlm_r2, a2v_r2, a4v_r2)   ! <r2^2> and <r2^4>
        a4_r1r2 = a4_joint(nlm_r1,nlm_r2) ! <r1^2 r2^2>
        
        ! Is nlm_r3 given? Then use it, else derive it from nlm_r1 and nlm_r2            
        if (real(nlm_r3(1)) > 1e-10) then
            call f_ev_ck_Mandel(nlm_r3, a2v_r3, a4v_r3)   ! <r3^2> and <r3^4>
            a4_r1r3 = a4_joint(nlm_r1,nlm_r3) ! <r1^2 r3^2>
            a4_r2r3 = a4_joint(nlm_r2,nlm_r3) ! <r2^2 r3^2>
        else
            a2v_r3 = mat_to_vec(a2_orth(nlm_r1, nlm_r2)) ! <r3^2> 
            a4v_r3 = a4_to_mat( a4_orth(nlm_r1, nlm_r2)) ! <r3^4> 
            a4_r1r3 = a4_jointcross(nlm_r1,nlm_r2) ! <r1^2 r3^2> 
            a4_r2r3 = a4_jointcross(nlm_r2,nlm_r1) ! <r2^2 r3^2> 
        end if
        
    end

    function a4_sym(a4) result(a4sym)
    
        implicit none
        
        real(kind=dp), intent(in) :: a4(3,3,3,3)
        real(kind=dp) :: a4sym(3,3,3,3)
        integer :: ii,jj,kk,ll ! loop indices
        
        do ii=1,3
            do jj=1,3
                do kk=1,3
                    do ll=1,3
                        a4sym(ii,jj,kk,ll) = (a4(ii,jj,kk,ll) + a4(kk,ll,ii,jj))/2
                    end do
                end do
            end do
        end do
    end

    function a4_altsym(a4) result(a4sym)
    
        implicit none
        
        real(kind=dp), intent(in) :: a4(3,3,3,3)
        real(kind=dp) :: a4sym(3,3,3,3)
        integer :: ii,jj,kk,ll ! loop indices
        
        do ii=1,3
            do jj=1,3
                do kk=1,3
                    do ll=1,3
                        a4sym(ii,jj,kk,ll) = (a4(ii,kk,jj,ll) + a4(kk,jj,ii,ll) + a4(ii,ll,kk,jj) + a4(ll,jj,kk,ii))/4
                    end do
                end do
            end do
        end do    
    end
    
    !---------------------------------
    ! ELASTIC
    !---------------------------------

    !function elas_fwd_tranisotropic_voigthomo(strain, nlm, lam,mu, Elam,Emu,Egam) result(strain)
    
        ! N/A
    
    !end

    function elas_rev_tranisotropic_reusshomo(strain, nlm, lam,mu, Elam,Emu,Egam) result(stress)

        ! Reuss-averaged elastic constitutive equation assuming transversely isotropic grains: sig(<eps>) (given nlm)
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: strain(3,3)
        real(kind=dp), intent(in)    :: lam,mu, Elam,Emu,Egam ! monocrystal parameters
        real(kind=dp)                :: stress(3,3) ! sigma

        real(kind=dp) :: k1,k2,k3,k4,k5 ! aux params
        real(kind=dp) :: a2mat(3,3), a2v(6), a4v(6,6)
        real(kind=dp) :: P(6,6), L(6,6), strain_vec(6,1) !,  P_reg(6,6),strain_vec_reg(6,1)
        integer       :: info

        call elas_fwdparams_tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5)
        
        call f_ev_ck_Mandel(nlm, a2v, a4v) ! Structure tensors in Mandel notation: a2v = ev_c2_Mandel, B = ev_c4_Mandel
        a2mat = vec_to_mat(a2v) ! = a2
        L = anticommutator_Mandel(a2mat)
        
        ! Matrix "P" of the vectorized bulk Reuss rheology: strain = matmul(P, stress)
        P = k1*outerprod6(identity_vec6,identity_vec6) + k2*identity6 &
            + k3*(outerprod6(identity_vec6,a2v) + outerprod6(a2v,identity_vec6)) + k4*a4v + k5*L 

        ! Solve *inverse* problem; we are seeking stress given strain in strain = matmul(P, stress).
        strain_vec(:,1) = mat_to_vec(strain)
        call dposv('L', 6, 1, P, 6, strain_vec, 6, info) ! strain_vec is now "stress_vec" solution. For some reason, "U" does not work.

        ! Ill posed? => Regularization needed.
        if (info /= 0) then
!            P_reg          = matmul(TRANSPOSE(P),P) + 1e-6*identity6
!            strain_vec_reg = matmul(TRANSPOSE(P),strain_vec)
!            call dposv('L', 6, 1, P_reg, 6, strain_vec_reg, 6, info)
!            strain_vec = strain_vec_reg
!            if (info /= 0) then
                stop 'specfab error: Reuss elasticity-matrix inversion failed! Please check the ODF (nlm) is correct.'
!            end if
        end if
        
        ! Revert solution to 3x3
        stress = vec_to_mat(strain_vec)
    end

end module homogenizations 
