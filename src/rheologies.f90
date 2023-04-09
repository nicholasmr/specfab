! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2023

! Viscoplastic rheologies

!-----------
! Notation:
!-----------
!   eps = Strain-rate tensor (3x3)
!   tau = Deviatoric stress tensor (3x3)
!   A   = Strain-rate (flow rate) factor
!   n   = Nonlinear flow exponent
!   mi  = Material/rheological symmetry axes (m1,m2,m3)
!   Eij = Directional enhancement factors w.r.t. axes "mi"
!-----------

module rheologies

    use tensorproducts

    implicit none 

    integer, parameter, private       :: dp = 8 ! Default precision
    real(kind=dp), parameter, private :: identity(3,3)  = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    
contains      

    !---------------------------------
    ! ISOTROPIC RHEOLOGY
    !---------------------------------

    function rheo_fwd_isotropic(tau, A,n) result(eps)

        implicit none
        real(kind=dp), intent(in)     :: tau(3,3), A
        integer, intent(in)           :: n
        real(kind=dp)                 :: eps(3,3), fluidity

        fluidity = A * (doubleinner22(tau,tau))**((n-1.0d0)/2.0d0)
        eps = fluidity * tau
    end

    function rheo_rev_isotropic(eps, A,n) result(tau)

        implicit none
        real(kind=dp), intent(in)     :: eps(3,3), A
        integer, intent(in)           :: n
        real(kind=dp)                 :: tau(3,3), viscosity

        viscosity = A**(-1.0d0/n) * (doubleinner22(eps,eps))**((1.0d0-n)/(2.0d0*n))
        tau = viscosity * eps
    end
    

    !---------------------------------
    ! COMPRESSIBLE (POROUS) ISOTROPIC RHEOLOGY
    !---------------------------------

    function rheo_fwd_isotropic_compressible(sig, A,n, fa,fb) result(eps)

        implicit none
        real(kind=dp), intent(in)     :: sig(3,3), A, fa,fb
        integer, intent(in)           :: n
        real(kind=dp)                 :: eps(3,3), fluidity, sigtr, I2

        sigtr = sig(1,1)+sig(2,2)+sig(3,3)
        I2 = doubleinner22(sig,sig)
        
        fluidity = A * (fa*(I2-(sigtr**2)/3) + 2.0d0/3*fb*(sigtr/3)**2)**((n-1.0d0)/2.0d0)
        eps = fluidity * ( fa*(sig - (sigtr/3)*identity) + 2.0d0/3*fb*(sigtr/3)*(identity/3) )
    end

    function rheo_rev_isotropic_compressible(eps, A,n, fa,fb) result(sig)

        implicit none
        real(kind=dp), intent(in)     :: eps(3,3), A, fa,fb
        integer, intent(in)           :: n
        real(kind=dp)                 :: sig(3,3), viscosity, epstr, J2

        epstr = eps(1,1)+eps(2,2)+eps(3,3)
        J2 = doubleinner22(eps,eps)
        
        viscosity = A**(-1.0d0/n) * ( 1/fa*(J2 - (epstr**2)/3) + 3.0d0**3/2*1/fb*(epstr/3)**2)**((1.0d0-n)/(2.0d0*n))
        sig = viscosity * ( 1/fa*(eps - (epstr/3)*identity) + 3.0d0**3/2*1/fb*(epstr/3)*(identity/3) )
    end

    !---------------------------------
    ! TRANSVERSELY ISOTROPIC RHEOLOGY
    !---------------------------------

    function rheo_fwd_tranisotropic(tau, A,n, m,Eij) result(eps)

        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, m(3), Eij(2)
        integer, intent(in)       :: n
        real(kind=dp)             :: eps(3,3), fluidity, kI,kM,kL, I2,I4,I5, mm(3,3),tausq(3,3)
        real(kind=dp), parameter  :: d = 3.0d0
        
        call rheo_params_tranisotropic(Eij,d,n,1.0d0, kI,kM,kL)

        mm = outerprod(m,m)
        tausq = matmul(tau,tau)    
        I2 = doubleinner22(tau,tau)
        I4 = doubleinner22(tau,mm)
        I5 = doubleinner22(tausq,mm)
        
        fluidity = A*(I2 + kM*I4**2 + 2*kL*I5)**((n-1.0d0)/2.0d0)
        eps = fluidity*( tau - kI*I4*identity + kM*I4*mm + kL*(matmul(tau,mm)+matmul(mm,tau)) )
    end

    function rheo_rev_tranisotropic(eps, A,n, m,Eij) result(tau)

        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, m(3), Eij(2)
        integer, intent(in)       :: n
        real(kind=dp)             :: tau(3,3), viscosity, kI,kM,kL, I2,I4,I5, mm(3,3),epssq(3,3)
        real(kind=dp), parameter  :: d = 3.0d0

        call rheo_params_tranisotropic(Eij,d,n,-1.0d0, kI,kM,kL)

        mm = outerprod(m,m)
        epssq = matmul(eps,eps)    
        I2 = doubleinner22(eps,eps)
        I4 = doubleinner22(eps,mm)
        I5 = doubleinner22(epssq,mm)

        viscosity = A**(-1.0d0/n)*(I2 + kM*I4**2 + 2*kL*I5)**((1.0d0-n)/(2.0d0*n))
        tau = viscosity*( eps - kI*I4*identity + kM*I4*mm + kL*(matmul(eps,mm)+matmul(mm,eps)) )
    end

    subroutine rheo_params_tranisotropic(Eij,d,n,expo, kI,kM,kL)

        implicit none
        real(kind=dp), intent(in)  :: Eij(2),expo,d
        integer, intent(in)        :: n
        real(kind=dp), intent(out) :: kI,kM,kL
        real(kind=dp)              :: nexpo

        nexpo = expo * 2.0d0/(n + 1.0d0) 
        
        ! Eij := (Emm,Emt)
        kI = (Eij(1)**nexpo-1)/(d-1.)
        kM = (d*(Eij(1)**nexpo+1)-2.)/(d-1.) - 2*Eij(2)**nexpo
        kL = Eij(2)**nexpo - 1
    end

    !---------------------------------
    ! ORTHOTROPIC RHEOLOGY
    !---------------------------------

    function rheo_fwd_orthotropic(tau, A,n, m1,m2,m3, Eij) result(eps)

        implicit none
        real(kind=dp), intent(in)     :: tau(3,3), A, m1(3),m2(3),m3(3), Eij(6)
        integer, intent(in)           :: n
        real(kind=dp)                 :: eps(3,3)
        real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
        real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
        real(kind=dp)                 :: I1,I2,I3,I4,I5,I6
        real(kind=dp)                 :: fluidity

        call rheo_params_orthotropic(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)
        call rheo_structs_orthotropic(tau, m1,m2,m3, M11,M22,M33,M23,M31,M12, I1,I2,I3,I4,I5,I6)

        fluidity = A * ( &
            + lam1 * I1**2 &
            + lam2 * I2**2 &
            + lam3 * I3**2 &
            + lam4 * I4**2 &
            + lam5 * I5**2 &
            + lam6 * I6**2 &
        )**((n-1.d0)/2.d0)
        
        eps = fluidity * ( &
            + lam1 * I1 * (M22 - M33)/2.0d0 &
            + lam2 * I2 * (M33 - M11)/2.0d0 &
            + lam3 * I3 * (M11 - M22)/2.0d0 &
            + lam4 * I4 * (M23 + transpose(M23))/2.0d0 &
            + lam5 * I5 * (M31 + transpose(M31))/2.0d0 &
            + lam6 * I6 * (M12 + transpose(M12))/2.0d0 &
        )
    end

    function rheo_rev_orthotropic(eps, A,n, m1,m2,m3, Eij) result(tau)

        implicit none
        real(kind=dp), intent(in)     :: eps(3,3), A, m1(3),m2(3),m3(3), Eij(6)
        integer, intent(in)           :: n
        real(kind=dp)                 :: tau(3,3)
        real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
        real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
        real(kind=dp)                 :: J1,J2,J3,J4,J5,J6, J23,J31,J12
        real(kind=dp)                 :: viscosity

        call rheo_params_orthotropic(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)
        call rheo_structs_orthotropic(eps, m1,m2,m3, M11,M22,M33,M23,M31,M12, J1,J2,J3,J4,J5,J6)
        call rheo_auxinvars_orthotropic(J1,J2,J3, J23,J31,J12)

        viscosity = A**(-1.d0/n) * ( &
            + lam1/gam * J23**2 &
            + lam2/gam * J31**2 & 
            + lam3/gam * J12**2 &
            + 4 * (1/lam4) * J4**2 &
            + 4 * (1/lam5) * J5**2 &
            + 4 * (1/lam6) * J6**2 &
        )**((1-n)/(2.d0*n))

        tau = viscosity * ( &
            + lam1/gam * J23 * (identity - 3*M11)/2 &
            + lam2/gam * J31 * (identity - 3*M22)/2 &
            + lam3/gam * J12 * (identity - 3*M33)/2 &
            + 4 * (1/lam4) * J4 * (M23 + transpose(M23))/2 &
            + 4 * (1/lam5) * J5 * (M31 + transpose(M31))/2 &
            + 4 * (1/lam6) * J6 * (M12 + transpose(M12))/2 &
        )
    end

    function rheo_rev_orthotropic_dimless(eps, n, m1,m2,m3, Eij) result(tau)

        implicit none
        real(kind=dp), intent(in)     :: eps(3,3), m1(3),m2(3),m3(3), Eij(6)
        integer, intent(in)           :: n
        real(kind=dp)                 :: tau(3,3)
        real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
        real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
        real(kind=dp)                 :: J1,J2,J3,J4,J5,J6, J23,J31,J12

        call rheo_params_orthotropic(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)
        call rheo_structs_orthotropic(eps, m1,m2,m3, M11,M22,M33,M23,M31,M12, J1,J2,J3,J4,J5,J6)
        call rheo_auxinvars_orthotropic(J1,J2,J3, J23,J31,J12)

        tau = ( &
            + lam1/gam * J23 * (identity - 3*M11)/2 &
            + lam2/gam * J31 * (identity - 3*M22)/2 &
            + lam3/gam * J12 * (identity - 3*M33)/2 &
            + 4 * (1/lam4) * J4 * (M23 + transpose(M23))/2 &
            + 4 * (1/lam5) * J5 * (M31 + transpose(M31))/2 &
            + 4 * (1/lam6) * J6 * (M12 + transpose(M12))/2 &
        )
    end

    subroutine rheo_params_orthotropic(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)

        implicit none
        real(kind=dp), intent(in)  :: Eij(6) ! E11,E22,E33,E23(E32),E13(E31),E12(E21)
        integer, intent(in)        :: n
        real(kind=dp), intent(out) :: lam1,lam2,lam3,lam4,lam5,lam6, gam
        real(kind=dp)              :: Eexpo, E11n,E22n,E33n !,E12n,E31n,E23n

        Eexpo = 2.0d0/(n + 1.0d0) 
        
        E11n = Eij(1)**Eexpo
        E22n = Eij(2)**Eexpo
        E33n = Eij(3)**Eexpo

        lam1 = 4/3.0d0*(-E11n+E22n+E33n)
        lam2 = 4/3.0d0*(+E11n-E22n+E33n)
        lam3 = 4/3.0d0*(+E11n+E22n-E33n)
        lam4 = 2* Eij(4)**Eexpo
        lam5 = 2* Eij(5)**Eexpo
        lam6 = 2* Eij(6)**Eexpo
        
        gam = -(E11n**2 + E22n**2 + E33n**2) + 2*(E11n*E22n + E11n*E33n + E22n*E33n)
        gam = abs(gam) 
    end

    subroutine rheo_structs_orthotropic(X, m1,m2,m3, M11,M22,M33,M23,M31,M12, I1,I2,I3,I4,I5,I6)

        implicit none    
        real(kind=dp), intent(in)  :: X(3,3), m1(3),m2(3),m3(3)
        real(kind=dp), intent(out) :: M11(3,3), M22(3,3), M33(3,3), M23(3,3), M31(3,3), M12(3,3)
        real(kind=dp), intent(out) :: I1, I2, I3, I4, I5, I6
        real(kind=dp)              :: X11, X22, X33, X23, X31, X12
       
        M11 = outerprod(m1,m1)
        M22 = outerprod(m2,m2)
        M33 = outerprod(m3,m3)
        M23 = outerprod(m2,m3)
        M31 = outerprod(m3,m1)
        M12 = outerprod(m1,m2)

        X11 = doubleinner22(X, M11) 
        X22 = doubleinner22(X, M22) 
        X33 = doubleinner22(X, M33) 
        X23 = doubleinner22(X, M23) 
        X31 = doubleinner22(X, M31) 
        X12 = doubleinner22(X, M12) 
        
        I1 = (X22 - X33)/2.0d0 
        I2 = (X33 - X11)/2.0d0 
        I3 = (X11 - X22)/2.0d0 
        I4 = X23 ! I4 := (X23 + X32)/2 (but X is symmetric)
        I5 = X31 ! I5 := (X31 + X13)/2 (but X is symmetric)
        I6 = X12 ! I6 := (X12 + X21)/2 (but X is symmetric)
    end

    subroutine rheo_auxinvars_orthotropic(J1,J2,J3, J23,J31,J12)

        implicit none    
        real(kind=dp), intent(in)  :: J1,J2,J3
        real(kind=dp), intent(out) :: J23,J31,J12
       
        J23 = J2-J3 ! J23 := J2-J3 = (X22 + X33 - 2*X11)/2.0d0 = -3/2.0d0 * X11
        J31 = J3-J1 ! J31 := J3-J1 = (X11 + X33 - 2*X22)/2.0d0 = -3/2.0d0 * X22
        J12 = J1-J2 ! J12 := J1-J2 = (X11 + X22 - 2*X33)/2.0d0 = -3/2.0d0 * X33
    end

    !---------------------------------
    ! ORTHOTROPIC RHEOLOGY 
    ! ... with Pettit's hypothesis that the *fluidity* is orientation independant.
    !---------------------------------

    function rheo_fwd_orthotropic_Pettit(tau, A,n, m1,m2,m3, Eij) result(eps)

        implicit none
        real(kind=dp), intent(in)     :: tau(3,3), A, m1(3),m2(3),m3(3), Eij(6)
        integer, intent(in)           :: n
        real(kind=dp)                 :: eps(3,3)
        real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
        real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
        real(kind=dp)                 :: I1,I2,I3,I4,I5,I6
        real(kind=dp)                 :: fluidity

        call rheo_params_orthotropic(Eij, 1, lam1,lam2,lam3,lam4,lam5,lam6, gam)
        call rheo_structs_orthotropic(tau, m1,m2,m3, M11,M22,M33,M23,M31,M12, I1,I2,I3,I4,I5,I6)

        fluidity = A * (doubleinner22(tau,tau))**((n-1)/2.0d0) ! Pettit's hypothesis: fluidity = Glen's isotropic fluidity 
        
        eps = fluidity * ( &
            + lam1 * I1 * (M22 - M33)/2.0d0 &
            + lam2 * I2 * (M33 - M11)/2.0d0 &
            + lam3 * I3 * (M11 - M22)/2.0d0 &
            + lam4 * I4 * (M23 + transpose(M23))/2.0d0 &
            + lam5 * I5 * (M31 + transpose(M31))/2.0d0 &
            + lam6 * I6 * (M12 + transpose(M12))/2.0d0 &
        )
    end

    function rheo_rev_orthotropic_Pettit(eps, A,n, m1,m2,m3, Eij) result(tau)

        implicit none
        real(kind=dp), intent(in)     :: eps(3,3), A, m1(3),m2(3),m3(3), Eij(6)
        integer, intent(in)           :: n
        real(kind=dp)                 :: tau(3,3)
        real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
        real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
        real(kind=dp)                 :: J1,J2,J3,J4,J5,J6, J23,J31,J12
        real(kind=dp)                 :: viscosity 

        call rheo_params_orthotropic(Eij, 1, lam1,lam2,lam3,lam4,lam5,lam6, gam)
        call rheo_structs_orthotropic(eps, m1,m2,m3, M11,M22,M33,M23,M31,M12, J1,J2,J3,J4,J5,J6)
        call rheo_auxinvars_orthotropic(J1,J2,J3, J23,J31,J12)

        viscosity = A**(-1.d0/n) * ( &
            + 3.0d0/2 * lam1**2/gam**2 * J23**2 &
            + 3.0d0/2 * lam2**2/gam**2 * J31**2 & 
            + 3.0d0/2 * lam3**2/gam**2 * J12**2 &
            - 3.0d0/2 * lam2*lam3/gam**2 * J12*J31 &
            - 3.0d0/2 * lam3*lam1/gam**2 * J12*J23 &
            - 3.0d0/2 * lam1*lam2/gam**2 * J23*J31 &
            + 8 * (1/lam4**2) * J4**2 &
            + 8 * (1/lam5**2) * J5**2 &
            + 8 * (1/lam6**2) * J6**2 &
        )**((1-n)/(2.d0*n))

        tau = viscosity * ( &
            + lam1/gam * J23 * (identity - 3*M11)/2 &
            + lam2/gam * J31 * (identity - 3*M22)/2 &
            + lam3/gam * J12 * (identity - 3*M33)/2 &
            + 4 * (1/lam4) * J4 * (M23 + transpose(M23))/2 &
            + 4 * (1/lam5) * J5 * (M31 + transpose(M31))/2 &
            + 4 * (1/lam6) * J6 * (M12 + transpose(M12))/2 &
        )
    end

    !---------------------------------
    ! ORTHOTROPIC RHEOLOGY 
    ! ... with Martin's hypothesis that the *viscosity* is orientation independant.
    !---------------------------------

    function rheo_fwd_orthotropic_Martin(tau, A,n, m1,m2,m3, Eij) result(eps)

        implicit none
        real(kind=dp), intent(in)     :: tau(3,3), A, m1(3),m2(3),m3(3), Eij(6)
        integer, intent(in)           :: n
        real(kind=dp)                 :: eps(3,3)
        real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
        real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
        real(kind=dp)                 :: I1,I2,I3,I4,I5,I6
        real(kind=dp)                 :: fluidity 
        real(kind=dp)                 :: c1,c2,c3, fc

        call rheo_params_orthotropic_Martin(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)
        call rheo_structs_orthotropic(tau, m1,m2,m3, M11,M22,M33,M23,M31,M12, I1,I2,I3,I4,I5,I6)

        fc = 1.0d0/16 
        c1 = fc * 4*(2*lam1**2 + lam1*(lam2+lam3) - lam2*lam3)
        c2 = fc * 4*(2*lam2**2 + lam2*(lam3+lam1) - lam3*lam1)
        c3 = fc * 4*(2*lam3**2 + lam3*(lam1+lam2) - lam1*lam2)

        fluidity = A * ( &
            + c1 * I1**2 &
            + c2 * I2**2 &
            + c3 * I3**2 &
            + 0.5d0 * lam4**2 * I4**2 &
            + 0.5d0 * lam5**2 * I5**2 &
            + 0.5d0 * lam6**2 * I6**2 &
        )**((n-1.d0)/2.d0)
        
        eps = fluidity * ( &
            + lam1 * I1 * (M22 - M33)/2.0d0 &
            + lam2 * I2 * (M33 - M11)/2.0d0 &
            + lam3 * I3 * (M11 - M22)/2.0d0 &
            + lam4 * I4 * (M23 + transpose(M23))/2.0d0 &
            + lam5 * I5 * (M31 + transpose(M31))/2.0d0 &
            + lam6 * I6 * (M12 + transpose(M12))/2.0d0 &
        )
    end

    function rheo_rev_orthotropic_Martin(eps, A,n, m1,m2,m3, Eij) result(tau)

        implicit none
        real(kind=dp), intent(in)     :: eps(3,3), A, m1(3),m2(3),m3(3), Eij(6)
        integer, intent(in)           :: n
        real(kind=dp)                 :: tau(3,3)
        real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
        real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
        real(kind=dp)                 :: J1,J2,J3,J4,J5,J6, J23,J31,J12
        real(kind=dp)                 :: viscosity

        call rheo_params_orthotropic_Martin(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)
        call rheo_structs_orthotropic(eps, m1,m2,m3, M11,M22,M33,M23,M31,M12, J1,J2,J3,J4,J5,J6)
        call rheo_auxinvars_orthotropic(J1,J2,J3, J23,J31,J12)

        viscosity = A**(-1.0d0/n) * (doubleinner22(eps,eps))**((1.0d0-n)/(2.0d0*n)) ! Martin's hypothesis: viscosity = Glen's isotropic viscosity

        tau = viscosity * ( &
            + lam1/gam * J23 * (identity - 3*M11)/2 &
            + lam2/gam * J31 * (identity - 3*M22)/2 &
            + lam3/gam * J12 * (identity - 3*M33)/2 &
            + 4 * (1/lam4) * J4 * (M23 + transpose(M23))/2 &
            + 4 * (1/lam5) * J5 * (M31 + transpose(M31))/2 &
            + 4 * (1/lam6) * J6 * (M12 + transpose(M12))/2 &
        )

    end

    subroutine rheo_params_orthotropic_Martin(Eij, n, lam1,lam2,lam3,lam4,lam5,lam6, gam)

        use lambdasolver
        
        implicit none
        real(kind=dp), intent(in)  :: Eij(6)
        integer, intent(in)        :: n
        real(kind=dp), intent(out) :: lam1,lam2,lam3,lam4,lam5,lam6, gam
        real(kind=dp)              :: Eexpo, lami(3), Eii(3), f1

        ! Logitudinal 
        Eii = [Eij(1),Eij(2),Eij(3)] ! Init guess
        lami = lambdaprime(n, Eii)
        lam1 = lami(1)
        lam2 = lami(2)
        lam3 = lami(3)
        f1 = 3.0d0/16
        gam = 3 * f1 * (lam2*lam3 + lam3*lam1 + lam1*lam2)

        ! Shear
        Eexpo = 1.0d0/n
        lam4 = 2 * Eij(4)**Eexpo
        lam5 = 2 * Eij(5)**Eexpo
        lam6 = 2 * Eij(6)**Eexpo
    end

end module rheologies

