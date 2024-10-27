! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2023

! Viscoplastic rheologies

!-----------
! Notation:
!-----------
!   eps = Strain-rate tensor (3x3)
!   tau = Deviatoric stress tensor (3x3)
!   A   = Flow-rate factor
!   n   = Flow exponent (n=1 => linear, n->inf => plastic)
!   mi  = Rheological symmetry axes (m1,m2,m3)
!   Eij = Directional enhancement factors w.r.t. mi axes
!-----------

module rheologies

    use header
    use tensorproducts

    implicit none 

contains      

    function powlawexp_fwd(n) result(expo)
        implicit none
        real(kind=dp), intent(in) :: n
        real(kind=dp) :: expo
        expo = (n-1)/2
    end
    
    function powlawexp_rev(n) result(expo)
        implicit none
        real(kind=dp), intent(in) :: n
        real(kind=dp) :: expo
        expo = (1-n)/(2*n)
    end

    !---------------------------------
    ! ISOTROPIC RHEOLOGY
    !---------------------------------

    function rheo_fwd_isotropic(tau, A,n) result(eps)

        implicit none
        real(kind=dp), intent(in)     :: tau(3,3), A, n
        real(kind=dp)                 :: eps(3,3), fluidity

        fluidity = A * (doubleinner22(tau,tau))**powlawexp_fwd(n)
        eps = fluidity * tau
    end

    function rheo_rev_isotropic(eps, A,n) result(tau)

        implicit none
        real(kind=dp), intent(in)     :: eps(3,3), A, n
        real(kind=dp)                 :: tau(3,3), viscosity

        viscosity = A**(-1/n) * (doubleinner22(eps,eps))**powlawexp_rev(n)
        tau = viscosity * eps
    end
    
    !---------------------------------
    ! COMPRESSIBLE (POROUS) ISOTROPIC RHEOLOGY
    !---------------------------------

    function rheo_fwd_isotropic_compressible(sig, A,n, fa,fb) result(eps)

        implicit none
        real(kind=dp), intent(in)     :: sig(3,3), A, n, fa,fb
        real(kind=dp)                 :: eps(3,3), fluidity, I1, I2

        I1 = sig(1,1)+sig(2,2)+sig(3,3)
        I2 = doubleinner22(sig,sig)
        
        fluidity = A * (fa/2*(I2-(I1**2)/3) + fb/3*(I1/3)**2)**powlawexp_fwd(n)
        eps = fluidity * ( fa*(sig - (I1/3)*identity) + 2.0d0/3*fb*(I1/3)*(identity/3) )
    end

    function rheo_rev_isotropic_compressible(eps, A,n, fa,fb) result(sig)

        implicit none
        real(kind=dp), intent(in)     :: eps(3,3), A, n, fa,fb
        real(kind=dp)                 :: sig(3,3), viscosity, J1, J2

        J1 = eps(1,1)+eps(2,2)+eps(3,3)
        J2 = doubleinner22(eps,eps)
        
        viscosity = A**(-1/n) * (1/(2*fa)*(J2 - (J1**2)/3) + 3.0d0/(4*fb)*J1**2)**powlawexp_rev(n)
        sig = viscosity * ( 1/fa*(eps - (J1/3)*identity) + 3.0d0/(2*fb)*J1*identity )
    end

    !---------------------------------
    ! TRANSVERSELY ISOTROPIC RHEOLOGY
    !---------------------------------

    function rheo_fwd_tranisotropic(tau, A,n, m1,Eij) result(eps)

        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, n, m1(3), Eij(2)
        real(kind=dp)             :: eps(3,3), fluidity, cI,cM,cL, M(3,3), I2,I4,I5
        
        call rheo_params_tranisotropic(Eij,3,n,1, cI,cM,cL)
        call rheo_structs_tranisotropic(tau,m1, M,I2,I4,I5)
        
        fluidity = A * (I2 + cM*I4**2 + 2*cL*I5)**powlawexp_fwd(n)
        eps = fluidity * (tau - cI*I4*identity + cM*I4*M + cL*anticommutator(tau,M))
    end

    function rheo_rev_tranisotropic(eps, A,n, m1,Eij) result(tau)

        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, n, m1(3), Eij(2)
        real(kind=dp)             :: tau(3,3), viscosity, cI,cM,cL, M(3,3), J2,J4,J5

        call rheo_params_tranisotropic(Eij,3,n,-1, cI,cM,cL)
        call rheo_structs_tranisotropic(eps,m1, M,J2,J4,J5)

        viscosity = A**(-1/n) * (J2 + cM*J4**2 + 2*cL*J5)**powlawexp_rev(n)
        tau = viscosity * (eps - cI*J4*identity + cM*J4*M + cL*anticommutator(eps,M))
    end

    subroutine rheo_params_tranisotropic(Eij, d, n, ef, cI,cM,cL)

        implicit none
        real(kind=dp), intent(in)  :: Eij(2), n ! Eij = (Emm,Emt)
        integer, intent(in)        :: d, ef
        real(kind=dp), intent(out) :: cI, cM, cL
        real(kind=dp)              :: ne

        ne = ef * 2/(n+1) ! ef is exponent (pre)factor
        cI = (Eij(1)**ne - 1)/(d-1)
        cM = (d*(Eij(1)**ne+1) - 2)/(d-1) - 2*Eij(2)**ne
        cL = Eij(2)**ne - 1
    end

    subroutine rheo_structs_tranisotropic(X,m1, M,I2,I4,I5)

        implicit none    
        real(kind=dp), intent(in)  :: X(3,3), m1(3) ! X = eps or tau, m1 = rotational symmetry axis
        real(kind=dp), intent(out) :: M(3,3), I2, I4, I5

        M = outerprod(m1,m1)
        I2 = doubleinner22(X,X)
        I4 = doubleinner22(X,M)
        I5 = doubleinner22(matmul(X,X),M)
    end

    !---------------------------------
    ! ORTHOTROPIC RHEOLOGY
    !---------------------------------

    function rheo_fwd_orthotropic(tau, A, n, m1,m2,m3, Eij) result(eps)

        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: eps(3,3), lami(6), ci(6), gam, Fij(6,3,3), Ii(6), fluidity

        call rheo_params_orthotropic(Eij, n, lami, gam)
        call rheo_structs_orthotropic(tau,m1,m2,m3, 'F', Fij,Ii)
        
        ci(1:3) = 4.0d0/3 * lami(1:3)
        ci(4:6) =       2 * lami(4:6)

        fluidity = A * sum(ci*Ii**2)**powlawexp_fwd(n)
        eps = fluidity * singleinner13(ci*Ii, Fij)

    end

    function rheo_rev_orthotropic(eps, A, n, m1,m2,m3, Eij) result(tau)

        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: tau(3,3), lami(6), ci(6), gam, Fij(6,3,3), Ii(6), viscosity

        call rheo_params_orthotropic(Eij, n, lami, gam)
        call rheo_structs_orthotropic(eps,m1,m2,m3, 'R', Fij, Ii)

        ci(1:3) = 4.0d0/3 * lami(1:3)/gam
        ci(4:6) =       2 * 1/lami(4:6)

        viscosity = A**(-1/n) * sum(ci*Ii**2)**powlawexp_rev(n)
        tau = viscosity * singleinner13(ci*Ii, Fij)
    end

    subroutine rheo_params_orthotropic(Eij, n, lami, gam)

        implicit none
        real(kind=dp), intent(in)  :: Eij(6), n ! Eij = (E11,E22,E33,E23,E13,E12)
        real(kind=dp), intent(out) :: lami(6), gam
        real(kind=dp)              :: Bij(6)

        Bij = Eij**(2/(n+1))
        
        lami(1) = -Bij(1)+Bij(2)+Bij(3)
        lami(2) = +Bij(1)-Bij(2)+Bij(3)
        lami(3) = +Bij(1)+Bij(2)-Bij(3)
        lami(4:6) = Bij(4:6)
        
        gam = 2*Bij(1)*Bij(2) + 2*Bij(1)*Bij(3) + 2*Bij(2)*Bij(3) - Bij(1)**2 - Bij(2)**2 - Bij(3)**2 ! = (4*A)^2 (triangle area)
        
!        if (lami(1)<0.0d0 .or. lami(2)<0.0d0 .or. lami(3)<0.0d0) then
!            print *, 'E_ii^(2/(n+1)) for i=1,2,3: ', Bij(1:3)
!            print *, 'rheo_params_orthotropic() error: E_ii^(2/(n+1)) for i=1,2,3 does not fulfill the triangle inequality, needed to guarantee positive-valued energy dissipation.'
!        end if
    end

    subroutine rheo_structs_orthotropic(tau, m1,m2,m3, OPT, Fij, Ii)

        implicit none    
        real(kind=dp), intent(in)    :: tau(3,3), m1(3),m2(3),m3(3)
        character(len=1), intent(in) :: OPT
        real(kind=dp), intent(out)   :: Fij(6,3,3), Ii(6)
        real(kind=dp)                :: Mij(6,3,3)
        integer :: jj
       
        Mij = Mij_orthotropic(m1,m2,m3)

        if (OPT == 'F') then 
            Fij(1,:,:) = (Mij(2,:,:) - Mij(3,:,:))/2
            Fij(2,:,:) = (Mij(3,:,:) - Mij(1,:,:))/2
            Fij(3,:,:) = (Mij(1,:,:) - Mij(2,:,:))/2
        else ! OPT == 'R'
            Fij(1,:,:) = (identity - 3*Mij(1,:,:))/2
            Fij(2,:,:) = (identity - 3*Mij(2,:,:))/2
            Fij(3,:,:) = (identity - 3*Mij(3,:,:))/2
        end if 
        
        Fij(4,:,:) = symmetricpart(Mij(4,:,:))
        Fij(5,:,:) = symmetricpart(Mij(5,:,:))
        Fij(6,:,:) = symmetricpart(Mij(6,:,:))

        Ii = [(doubleinner22(tau, Fij(jj,:,:)), jj=1,6)]
    end
    
    function Mij_orthotropic(m1,m2,m3) result(Mij)

        implicit none 
        real(kind=dp), intent(in) :: m1(3),m2(3),m3(3)
        real(kind=dp)             :: Mij(6,3,3)

        Mij(1,:,:) = outerprod(m1,m1)
        Mij(2,:,:) = outerprod(m2,m2)
        Mij(3,:,:) = outerprod(m3,m3)
        Mij(4,:,:) = outerprod(m2,m3)
        Mij(5,:,:) = outerprod(m3,m1)
        Mij(6,:,:) = outerprod(m1,m2)
    end

    !---------------------------------
    ! ORTHOTROPIC RHEOLOGY 
    ! ... with Pettit's hypothesis that the *fluidity* is orientation independent
    !---------------------------------

    function rheo_fwd_orthotropic_Pettit(tau, A, n, m1,m2,m3, Eij) result(eps)

        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: eps(3,3), lami(6), ci(6), gam, Fij(6,3,3), Ii(6), fluidity

        call rheo_params_orthotropic(Eij, 1.0d0, lami, gam) ! note n = 1 for Eij exponents
        call rheo_structs_orthotropic(tau,m1,m2,m3, 'F', Fij,Ii)

        ci(1:3) = 4.0d0/3 * lami(1:3)
        ci(4:6) =       2 * lami(4:6)

        fluidity = A * (doubleinner22(tau,tau))**powlawexp_fwd(n) ! Pettit's hypothesis: fluidity = Glen's isotropic fluidity 
        eps = fluidity * singleinner13(ci*Ii, Fij)
    end

    function rheo_rev_orthotropic_Pettit(eps, A, n, m1,m2,m3, Eij) result(tau)

        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: tau(3,3), lami(6), ci(6), cvi(6), gam, Fij(6,3,3), Ii(6), viscosity, vx

        call rheo_params_orthotropic(Eij, 1.0d0, lami, gam)
        call rheo_structs_orthotropic(eps,m1,m2,m3, 'R', Fij, Ii)

        ci(1:3) = 4.0d0/3 * lami(1:3)/gam
        ci(4:6) =       2 * 1/lami(4:6)
        
        cvi(1:3) = sqrt(3.0d0/2)*ci(1:3)
        cvi(4:6) = sqrt(1.0d0/2)*ci(4:6)

        vx = cvi(2)*cvi(3)*Ii(2)*Ii(3) + cvi(3)*cvi(1)*Ii(3)*Ii(1) + cvi(1)*cvi(2)*Ii(1)*Ii(2) ! extra term appearing in viscosity expression
        viscosity = A**(-1/n) * (sum(cvi**2*Ii**2) - vx)**powlawexp_rev(n)
        tau = viscosity * singleinner13(ci*Ii, Fij)
    end

    !---------------------------------
    ! ORTHOTROPIC RHEOLOGY 
    ! ... with Martin's hypothesis that the *viscosity* is orientation independent
    !---------------------------------

    function rheo_fwd_orthotropic_Martin(tau, A, n, m1,m2,m3, Eij) result(eps)

        implicit none
        real(kind=dp), intent(in) :: tau(3,3), A, n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: eps(3,3), lami(6), ci(6), cvi(6), gam, Fij(6,3,3), Ii(6), fluidity

        call rheo_params_orthotropic_Martin(Eij, n, lami, gam) 
        call rheo_structs_orthotropic(tau,m1,m2,m3, 'F', Fij,Ii)
        
        ci(1:6) = lami(1:6)
        
        ! lambda prime, eqn. (15) in Rathmann & Lilien (2022)
        cvi(1) = 1.0d0/2 * (ci(1)**2 + ci(1)*(ci(2)+ci(3))/2 - ci(2)*ci(3)/2)
        cvi(2) = 1.0d0/2 * (ci(2)**2 + ci(2)*(ci(3)+ci(1))/2 - ci(3)*ci(1)/2)
        cvi(3) = 1.0d0/2 * (ci(3)**2 + ci(3)*(ci(1)+ci(2))/2 - ci(1)*ci(2)/2)
        cvi(4:6) = 0.5d0 * ci(4:6)**2

        fluidity = A * sum(cvi*Ii**2)**powlawexp_fwd(n)
        eps = fluidity * singleinner13(ci*Ii, Fij)
    end

    function rheo_rev_orthotropic_Martin(eps, A, n, m1,m2,m3, Eij) result(tau)

        implicit none
        real(kind=dp), intent(in) :: eps(3,3), A, n, m1(3),m2(3),m3(3), Eij(6)
        real(kind=dp)             :: tau(3,3), lami(6), ci(6), gam, Fij(6,3,3), Ii(6), viscosity

        call rheo_params_orthotropic_Martin(Eij, n, lami, gam)
        call rheo_structs_orthotropic(eps,m1,m2,m3, 'R', Fij, Ii)

        ci(1:3) = lami(1:3)/gam
        ci(4:6) = 4/lami(4:6)

        viscosity = A**(-1/n) * doubleinner22(eps,eps)**powlawexp_rev(n) ! Martin's hypothesis: viscosity = Glen's isotropic viscosity
        tau = viscosity * singleinner13(ci*Ii, Fij)
    end

    subroutine rheo_params_orthotropic_Martin(Eij, n, lami, gam)

        use lambdasolver
        
        implicit none
        real(kind=dp), intent(in)  :: Eij(6), n
        real(kind=dp), intent(out) :: lami(6), gam

        lami(1:3) = lambdaprime(n, [Eij(1),Eij(2),Eij(3)])
        lami(4:6) = 2 * Eij(4:6)**(1/n)
        gam = 9.0d0/16 * (lami(2)*lami(3) + lami(3)*lami(1) + lami(1)*lami(2))
    end

end module rheologies

