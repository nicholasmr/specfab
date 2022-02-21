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

    function ev_epsprime_Sac(tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca,nprime) 

        ! Sachs grain-averaged rheology
        
        implicit none
        
        real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
        real(kind=dp), intent(in) :: Ecc, Eca, tau(3,3)
        integer, intent(in)       :: nprime
        real(kind=dp), parameter  :: d = 3.0d0
        real(kind=dp)             :: ev_etac0 = 0, ev_etac2(3,3), ev_etac4(3,3,3,3)
        real(kind=dp)             :: ev_epsprime_Sac(3,3), coefA,coefB,coefC, tausq(3,3), I2

        call tranisotropic_coefs(Ecc,Eca,d,nprime,1.0d0, coefA,coefB,coefC)

        I2 = doubleinner22(tau,tau)
        tausq = matmul(tau,tau)

        ! Linear grain fluidity    
        if (nprime .eq. 1) then
            ev_etac0 = 1.0
            ev_etac2 = ev_c2
            ev_etac4 = ev_c4

        ! Nonlinear orientation-dependent grain fluidity (Johnson's rheology)        
        else if (nprime .eq. 3) then
            ev_etac0 = I2*  1.0 + coefB*doubleinner22(doubleinner42(ev_c4,tau),tau) + 2*coefC*doubleinner22(ev_c2,tausq)
            ev_etac2 = I2*ev_c2 + coefB*doubleinner42(doubleinner62(ev_c6,tau),tau) + 2*coefC*doubleinner42(ev_c4,tausq)
            ev_etac4 = I2*ev_c4 + coefB*doubleinner62(doubleinner82(ev_c8,tau),tau) + 2*coefC*doubleinner62(ev_c6,tausq)

        ! Same as n'=3 but *without* the orientation-dependent terms in the nonlinear grain fluidity.        
        else if (nprime .eq. -3) then 
            ev_etac0 = I2*  1.0
            ev_etac2 = I2*ev_c2
            ev_etac4 = I2*ev_c4
            
        end if

        ev_epsprime_Sac = ev_etac0*tau - coefA*doubleinner22(ev_etac2,tau)*identity + coefB*doubleinner42(ev_etac4,tau) + coefC*(matmul(tau,ev_etac2)+matmul(ev_etac2,tau))
    end

    function ev_epsprime_Tay(tau, ev_c2,ev_c4, Ecc,Eca,nprime) 
        
        ! Taylor grain-averaged rheology
        ! NOTE: Taylor model supports only n'=1. nprime is required anyway for future compatibility with n'>1.
        
        implicit none
        
        real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3)
        real(kind=dp), intent(in) :: Ecc, Eca, tau(3,3)
        integer, intent(in)       :: nprime ! Dummy variable: <eps'(tau)> with Taylor hypothesis is implemented only for n' = 1.
        real(kind=dp), parameter  :: d = 3.0d0
        real(kind=dp)             :: P(6,6), tau_vec(6,1),  P_reg(6,6),tau_vec_reg(6,1)
        real(kind=dp)             :: ev_epsprime_Tay(3,3), coefA,coefB,coefC
        integer                   :: info
    !    integer :: ipiv(9), work

        real(kind=dp), parameter  :: identity6(6,6) = reshape([1, 0, 0, 0, 0, 0, &
                                                               0, 1, 0, 0, 0, 0, & 
                                                               0, 0, 1, 0, 0, 0, & 
                                                               0, 0, 0, 1, 0, 0, & 
                                                               0, 0, 0, 0, 1, 0, & 
                                                               0, 0, 0, 0, 0, 1], [6,6])

        call tranisotropic_coefs(Ecc,Eca,d,nprime,-1.0d0, coefA,coefB,coefC)

        ! Set P (6x6 matrix)
        P(1,1)=(1.0) + (-1.0)*ev_c2(1,1)*coefA + ev_c4(1,1,1,1)*coefB + (2.0)*ev_c2(1,1)*coefC
        P(1,2)=(-1.0)*ev_c2(2,2)*coefA + ev_c4(1,1,2,2)*coefB
        P(1,3)=(-1.0)*ev_c2(3,3)*coefA + ev_c4(1,1,3,3)*coefB
        P(1,4)=(-1.0)*Sqrt((2.0))*ev_c2(2,3)*coefA + Sqrt((2.0))*ev_c4(1,1,2,3)*coefB
        P(1,5)=(-1.0)*Sqrt((2.0))*ev_c2(1,3)*coefA + Sqrt((2.0))*ev_c4(1,1,3,1)*coefB + Sqrt((2.0))*ev_c2(1,3)*coefC
        P(1,6)=(-1.0)*Sqrt((2.0))*ev_c2(1,2)*coefA + Sqrt((2.0))*ev_c4(1,1,1,2)*coefB + Sqrt((2.0))*ev_c2(1,2)*coefC
        P(2,1)=(-1.0)*ev_c2(1,1)*coefA + ev_c4(2,2,1,1)*coefB
        P(2,2)=(1.0) + (-1.0)*ev_c2(2,2)*coefA + ev_c4(2,2,2,2)*coefB + (2.0)*ev_c2(2,2)*coefC
        P(2,3)=(-1.0)*ev_c2(3,3)*coefA + ev_c4(2,2,3,3)*coefB
        P(2,4)=(-1.0)*Sqrt((2.0))*ev_c2(2,3)*coefA + Sqrt((2.0))*ev_c4(2,2,2,3)*coefB + Sqrt((2.0))*ev_c2(2,3)*coefC
        P(2,5)=(-1.0)*Sqrt((2.0))*ev_c2(1,3)*coefA + Sqrt((2.0))*ev_c4(2,2,3,1)*coefB
        P(2,6)=(-1.0)*Sqrt((2.0))*ev_c2(1,2)*coefA + Sqrt((2.0))*ev_c4(2,2,1,2)*coefB + Sqrt((2.0))*ev_c2(1,2)*coefC
        P(3,1)=(-1.0)*ev_c2(1,1)*coefA + ev_c4(3,3,1,1)*coefB
        P(3,2)=(-1.0)*ev_c2(2,2)*coefA + ev_c4(3,3,2,2)*coefB
        P(3,3)=(1.0) + (-1.0)*ev_c2(3,3)*coefA + ev_c4(3,3,3,3)*coefB + (2.0)*ev_c2(3,3)*coefC
        P(3,4)=(-1.0)*Sqrt((2.0))*ev_c2(2,3)*coefA + Sqrt((2.0))*ev_c4(3,3,2,3)*coefB + Sqrt((2.0))*ev_c2(2,3)*coefC
        P(3,5)=(-1.0)*Sqrt((2.0))*ev_c2(1,3)*coefA + Sqrt((2.0))*ev_c4(3,3,3,1)*coefB + Sqrt((2.0))*ev_c2(1,3)*coefC
        P(3,6)=(-1.0)*Sqrt((2.0))*ev_c2(1,2)*coefA + Sqrt((2.0))*ev_c4(3,3,1,2)*coefB
        P(4,1)=Sqrt((2.0))*ev_c4(2,3,1,1)*coefB
        P(4,2)=Sqrt((2.0))*ev_c4(2,3,2,2)*coefB + Sqrt((2.0))*ev_c2(2,3)*coefC
        P(4,3)=Sqrt((2.0))*ev_c4(2,3,3,3)*coefB + Sqrt((2.0))*ev_c2(2,3)*coefC
        P(4,4)=(1.0) + (2.0)*ev_c4(2,3,2,3)*coefB + (ev_c2(2,2) + ev_c2(3,3))*coefC
        P(4,5)=(2.0)*ev_c4(2,3,3,1)*coefB + ev_c2(1,2)*coefC
        P(4,6)=(2.0)*ev_c4(2,3,1,2)*coefB + ev_c2(1,3)*coefC
        P(5,1)=Sqrt((2.0))*ev_c4(3,1,1,1)*coefB + Sqrt((2.0))*ev_c2(1,3)*coefC
        P(5,2)=Sqrt((2.0))*ev_c4(3,1,2,2)*coefB
        P(5,3)=Sqrt((2.0))*ev_c4(3,1,3,3)*coefB + Sqrt((2.0))*ev_c2(1,3)*coefC
        P(5,4)=(2.0)*ev_c4(3,1,2,3)*coefB + ev_c2(1,2)*coefC
        P(5,5)=(1.0) + (2.0)*ev_c4(3,1,3,1)*coefB + (ev_c2(1,1) + ev_c2(3,3))*coefC
        P(5,6)=(2.0)*ev_c4(3,1,1,2)*coefB + ev_c2(2,3)*coefC
        P(6,1)=Sqrt((2.0))*ev_c4(1,2,1,1)*coefB + Sqrt((2.0))*ev_c2(1,2)*coefC
        P(6,2)=Sqrt((2.0))*ev_c4(1,2,2,2)*coefB + Sqrt((2.0))*ev_c2(1,2)*coefC
        P(6,3)=Sqrt((2.0))*ev_c4(1,2,3,3)*coefB
        P(6,4)=(2.0)*ev_c4(1,2,2,3)*coefB + ev_c2(1,3)*coefC
        P(6,5)=(2.0)*ev_c4(1,2,3,1)*coefB + ev_c2(2,3)*coefC
        P(6,6)=(1.0) + (2.0)*ev_c4(1,2,1,2)*coefB + (ev_c2(1,1) + ev_c2(2,2))*coefC

        tau_vec(:,1) = mat_to_vec(tau)
        call dposv('L', 6, 1, P, 6, tau_vec, 6, info) ! tau_vec is now "eps_vec" solution. For some reason, "U" does not work..

        if (info /= 0) then
            P_reg       = matmul(TRANSPOSE(P),P) + 1e-6*identity6
            tau_vec_reg = matmul(TRANSPOSE(P),tau_vec)
            call dposv('L', 6, 1, P_reg, 6, tau_vec_reg, 6, info)
            tau_vec = tau_vec_reg
            if (info /= 0) then
                stop 'specfab error: Taylor viscosity-matrix inversion failed! Please check the ODF is correct (reducing the fabric integration time-step, and/or increasing regularization, for transient problems can often help).'
            end if
        end if
        
        ev_epsprime_Tay = vec_to_mat(tau_vec)
    end

    function eps_of_tau__linearTaylorSachs(tau, nlm, Aprime,Ecc,Eca,alpha) result(eps)

        ! Mixed Taylor--Sachs grain-averaged rheology:
        !       eps = (1-alpha)*<eps'(tau)> + alpha*eps(<tau'>)
        ! ...where eps'(tau') is the transversely isotropic rheology for n'=1 (linear viscous)

        implicit none

        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: tau(3,3), Aprime, Ecc, Eca, alpha
        real(kind=dp)                :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
        real(kind=dp)                :: eps(3,3)
        
        ! Linear grain rheology (n'=1) relies only on <c^2> and <c^4>
        call f_ev_ck(nlm, 'r', ev_c2,ev_c4,ev_c6,ev_c8) ! Calculate structure tensors of orders 2,4 (6,8 are assumed zero for faster evaluation)

        eps = (1-alpha)*Aprime*ev_epsprime_Sac(tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca,1)  &
                + alpha*Aprime*ev_epsprime_Tay(tau, ev_c2,ev_c4,             Ecc,Eca,1) 
    end
    
    !---------------------------------
    ! SOLID
    !---------------------------------

    ! N/A

end module homogenizations 
