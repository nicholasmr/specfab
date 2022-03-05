subroutine dndt_to_daidt(ddt_nlm, n00, da2dt, da4dt)

    complex(kind=dp), intent(in) :: ddt_nlm(nlm_len), n00
    real(kind=dp), intent(out)   :: da2dt(3,3), da4dt(3,3,3,3)
    complex(kind=dp)             :: ddt_n2m(-2:2), ddt_n4m(-4:4)
    
    ddt_n2m = ddt_nlm(I_l2:(I_l4-1))
    ddt_n4m = ddt_nlm(I_l4:(I_l6-1))
    
    da2dt = f_ev_c2(n00, ddt_n2m)          - f_ev_c2(n00, 0*ddt_n2m)            ! <c^> entries are *linear* combinations of n_l^m, allowing d/dt <c^2> to be calculated by substituting "d/dt n_l^m" for "n_l^m" in "<c^2>(n_l^m)" by removing the constant part(s).
    da4dt = f_ev_c4(n00, ddt_n2m, ddt_n4m) - f_ev_c4(n00, 0*ddt_n2m, 0*ddt_n4m) ! ...same reasoning as above
end

subroutine daidt_LATROT(eps, omg, a2, a4, da2dt, da4dt)
    
    ! Lattice rotation contribution to d/dt a^(2) and d/dt a^(4) 
    ! Assumes Taylor hypothesis (constant strain-rate over polycrystal).
    
    implicit none

    real(kind=dp), intent(in)  :: eps(3,3), omg(3,3), a2(3,3), a4(3,3,3,3)
    real(kind=dp), intent(out) :: da2dt(3,3), da4dt(3,3,3,3)
    complex(kind=dp)           :: nlm(nlm_len), ddt_nlm(nlm_len)
    
    nlm = 0.0
    nlm(1:(I_l6-1)) = a4_to_nlm(a2, a4) ! tensorial --> spectral, truncated at L=4
    ddt_nlm = matmul(dndt_ij_LATROT(eps,omg, 0*eps,0d0,0d0,0d0, 1d0), nlm) ! spectral evolution  --- d/dt nlm_i = M_ij nlm_j 
    call dndt_to_daidt(ddt_nlm, nlm(1), da2dt, da4dt) ! spectral --> tensorial 
end

subroutine daidt_DDRX(tau, Gamma0, a2, a4, da2dt, da4dt)
    
    ! DDRX contribution to d/dt a^(2) and d/dt a^(4) 
    
    implicit none

    real(kind=dp), intent(in)  :: tau(3,3), Gamma0, a2(3,3), a4(3,3,3,3)
    real(kind=dp), intent(out) :: da2dt(3,3), da4dt(3,3,3,3)
    complex(kind=dp)           :: nlm(nlm_len), ddt_nlm(nlm_len)
    
    nlm = 0.0
    nlm(1:(I_l6-1)) = a4_to_nlm(a2, a4) ! tensorial --> spectral, truncated at L=4
    ddt_nlm = Gamma0 * matmul(dndt_ij_DDRX(nlm, tau), nlm) ! spectral evolution  --- d/dt nlm_i = M_ij nlm_j 
    call dndt_to_daidt(ddt_nlm, nlm(1) , da2dt, da4dt) ! spectral --> tensorial 
end

subroutine daidt_REG(eps, a2, a4, da2dt, da4dt)
    
    ! Regularization contribution to d/dt a^(2) and d/dt a^(4) 
    
    implicit none

    real(kind=dp), intent(in)  :: eps(3,3), a2(3,3), a4(3,3,3,3)
    real(kind=dp), intent(out) :: da2dt(3,3), da4dt(3,3,3,3)
    complex(kind=dp)           :: nlm(nlm_len), ddt_nlm(nlm_len)
    
    nlm = 0.0
    nlm(1:(I_l6-1)) = a4_to_nlm(a2, a4) ! tensorial --> spectral, truncated at L=4
    ddt_nlm = matmul(dndt_ij_REG(eps), nlm) ! spectral evolution  --- d/dt nlm_i = M_ij nlm_j  
    call dndt_to_daidt(ddt_nlm, nlm(1), da2dt, da4dt) ! spectral --> tensorial 
end

function a6_CBT(a2,a4)

    ! "Clusure By Trancation" for a^(6) in terms of a^(2) and a^(4)

    implicit none
    
    real(kind=dp), intent(in) :: a2(3,3), a4(3,3,3,3)
    real(kind=dp)             :: a6_CBT(3,3,3,3,3,3)
    complex(kind=dp)          :: nlm(nlm_len), n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        
    nlm = 0.0
    nlm(1:(I_l6-1)) = a4_to_nlm(a2, a4) ! Reconstruct n_l^m

    n00 = nlm(1)
    n2m = nlm(I_l2:(I_l4-1))
    n4m = nlm(I_l4:(I_l6-1))
    n6m = 0 
    
    a6_CBT = f_ev_c6(n00, n2m, n4m, n6m) ! Calculate a^(6) given the spectral truncation n_l^m = 0 for l>4
end

