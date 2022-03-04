! N. M. Rathmann <rathmann@nbi.ku.dk>, 2019-2022

module homogenizations  

    use tensorproducts
    use mandel
    use moments
    use rheologies

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    integer, parameter :: identity(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    
    ! Optimal n'=1 (lin) grain parameters 
    ! These are the linear mixed Taylor--Sachs best-fit parameters from Rathmann and Lilien (2021)
    real(kind=dp), parameter :: Eca_opt_lin   = 1d3
    real(kind=dp), parameter :: Ecc_opt_lin   = 1d0
    real(kind=dp), parameter :: alpha_opt_lin = 0.0125
    
    ! Optimal n'=3 (nlin) grain parameters 
    ! These are the nonlinear Sachs-only best-fit parameters (Rathmann et al., 2021) 
    real(kind=dp), parameter :: Eca_opt_nlin   = 1d4
    real(kind=dp), parameter :: Ecc_opt_nlin   = 1d0
    real(kind=dp), parameter :: alpha_opt_nlin = 0
        
contains      

    !---------------------------------
    ! FLUID
    !---------------------------------

    function ev_epsprime_Sac(tau, nlm, Ecc,Eca,nprime) result(eps)

        ! Sachs grain-averaged rheology
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: Ecc, Eca, tau(3,3)
        integer, intent(in)          :: nprime
        
        real(kind=dp), parameter  :: d = 3.0d0
        real(kind=dp)             :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
        real(kind=dp)             :: ev_etac0 = 0, ev_etac2(3,3), ev_etac4(3,3,3,3) ! <eta*c^k> terms
        real(kind=dp)             :: eps(3,3), coefA,coefB,coefC, tausq(3,3), I2, B3(3,3), B(6,6), Mv(6)

        call tranisotropic_coefs(Ecc,Eca,d,nprime,1.0d0, coefA,coefB,coefC)

        I2 = doubleinner22(tau,tau)

        ! Linear grain fluidity    
        if (nprime .eq. 1) then
            call f_ev_ck_Mandel(nlm, Mv, B) ! Structure tensors in Mandel notation: Mv = ev_c2_Mandel, B = ev_c4_Mandel
            ev_etac0 = 1.0d0
            ev_etac2 = vec_to_mat(Mv) ! = a2 = ev_c2
            B3 = vec_to_mat(matmul(B,mat_to_vec(tau))) ! = doubleinner42(ev_etac4,tau)

        ! Nonlinear orientation-dependent grain fluidity (Johnson's rheology)        
        else if (nprime .eq. 3) then
            tausq = matmul(tau,tau)
            if (size(nlm) < (8 +1)*(8 +2)/2) then ! L < 8 ?
                stop "specfab error: Sachs homogenization with n'=3 requires L >= 8."
            end if
            call f_ev_ck(nlm, 'f', ev_c2,ev_c4,ev_c6,ev_c8) ! Calculate structure tensors of orders 2,4,6,8
            ev_etac0 = I2*  1.0 + coefB*doubleinner22(doubleinner42(ev_c4,tau),tau) + 2*coefC*doubleinner22(ev_c2,tausq)
            ev_etac2 = I2*ev_c2 + coefB*doubleinner42(doubleinner62(ev_c6,tau),tau) + 2*coefC*doubleinner42(ev_c4,tausq)
            ev_etac4 = I2*ev_c4 + coefB*doubleinner62(doubleinner82(ev_c8,tau),tau) + 2*coefC*doubleinner62(ev_c6,tausq)
            B3 = doubleinner42(ev_etac4,tau)
            
        ! Same as n'=3 but *without* the orientation-dependent terms in the nonlinear grain fluidity.        
        else if (nprime .eq. -3) then 
            call f_ev_ck_Mandel(nlm, Mv, B) ! Structure tensors in Mandel notation: Mv = ev_c2_Mandel, B = ev_c4_Mandel
            ev_etac0 = I2 * 1.0d0
            ev_etac2 = I2 * vec_to_mat(Mv) ! = I2*a2 = I2*ev_c2
            B3       = I2 * vec_to_mat(matmul(B,mat_to_vec(tau))) ! = I2 * doubleinner42(ev_etac4,tau)
            
        end if

        eps = ev_etac0*tau - coefA*doubleinner22(ev_etac2,tau)*identity + coefB*B3 + coefC*(matmul(tau,ev_etac2)+matmul(ev_etac2,tau))
    end

    function ev_epsprime_Tay(tau, nlm, Ecc,Eca,nprime) result(eps)
        
        ! Taylor grain-averaged rheology
        ! NOTE: Taylor model supports only n'=1. nprime is required anyway for future compatibility with n'>1.
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: Ecc, Eca, tau(3,3)
        integer, intent(in)          :: nprime ! Dummy variable: <eps'(tau)> with Taylor hypothesis is implemented only for n' = 1.
        
        real(kind=dp), parameter  :: d = 3.0d0, s = sqrt(2.0)
        real(kind=dp)             :: Mv(6), M(3,3), B(6,6)
        real(kind=dp)             :: P(6,6), L(6,6), tau_vec(6,1),  P_reg(6,6),tau_vec_reg(6,1)
        real(kind=dp)             :: eps(3,3), coefA,coefB,coefC
        integer                   :: info
    !    integer :: ipiv(9), work

        real(kind=dp), parameter  :: identity_vec6(6) =  [1,1,1, 0,0,0] ! = mat_to_vec(identity)
        real(kind=dp), parameter  :: identity6(6,6) = reshape([1, 0, 0, 0, 0, 0, &
                                                               0, 1, 0, 0, 0, 0, & 
                                                               0, 0, 1, 0, 0, 0, & 
                                                               0, 0, 0, 1, 0, 0, & 
                                                               0, 0, 0, 0, 1, 0, & 
                                                               0, 0, 0, 0, 0, 1], [6,6])
                                                               
        call tranisotropic_coefs(Ecc,Eca,d,nprime,-1.0d0, coefA,coefB,coefC)
        
        call f_ev_ck_Mandel(nlm, Mv, B) ! Structure tensors in Mandel notation: Mv = ev_c2_Mandel, B = ev_c4_Mandel
        M = vec_to_mat(Mv) ! = a2
                                                
        ! L in Mandel notation
        L(1,:) = [2*M(1,1), 0.0d0, 0.0d0, 0.0d0, s*M(1,3), s*M(1,2)]
        L(2,:) = [0.0d0, 2*M(2,2), 0.0d0, s*M(2,3), 0.0d0, s*M(1,2)]
        L(3,:) = [0.0d0, 0.0d0, 2*M(3,3), s*M(2,3), s*M(1,3), 0.0d0]
        L(4,:) = [0.0d0, s*M(2,3), s*M(2,3), M(2,2)+M(3,3), M(1,2), M(1,3)]
        L(5,:) = [s*M(1,3), 0.0d0, s*M(1,3), M(1,2), M(1,1)+M(3,3), M(2,3)]
        L(6,:) = [s*M(1,2), s*M(1,2), 0.0d0, M(1,3), M(2,3), M(1,1)+M(2,2)]
        
        P = identity6 - coefA*outerprod6(identity_vec6,Mv) + coefB*B + coefC*L ! matrix of vectorized bulk rheology: tau = matmul(P, eps)
        
        ! Inverse
        tau_vec(:,1) = mat_to_vec(tau)
        call dposv('L', 6, 1, P, 6, tau_vec, 6, info) ! tau_vec is now "eps_vec" solution. For some reason, "U" does not work.

        if (info /= 0) then
            P_reg       = matmul(TRANSPOSE(P),P) + 1e-6*identity6
            tau_vec_reg = matmul(TRANSPOSE(P),tau_vec)
            call dposv('L', 6, 1, P_reg, 6, tau_vec_reg, 6, info)
            tau_vec = tau_vec_reg
            if (info /= 0) then
                stop 'specfab error: Taylor viscosity-matrix inversion failed! Please check the ODF is correct (reducing the fabric integration time-step, and/or increasing regularization, for transient problems can often help).'
            end if
        end if
        
        eps = vec_to_mat(tau_vec)
    end

    function eps_of_tau__linearTaylorSachs(tau, nlm, Aprime,Ecc,Eca,alpha) result(eps)

        ! Mixed Taylor--Sachs grain-averaged rheology:
        !       eps = (1-alpha)*<eps'(tau)> + alpha*eps(<tau'>)
        ! ...where eps'(tau') is the transversely isotropic rheology for n'=1 (linear viscous)

        implicit none

        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: tau(3,3), Aprime, Ecc, Eca, alpha
        real(kind=dp)                :: eps(3,3)
        
        eps = (1-alpha)*Aprime*ev_epsprime_Sac(tau, nlm, Ecc,Eca,1)  &
                + alpha*Aprime*ev_epsprime_Tay(tau, nlm, Ecc,Eca,1) 
    end
    
    !---------------------------------
    ! SOLID
    !---------------------------------

    ! N/A

end module homogenizations 
