! N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

! Solid elasticities (constitutive equations)

module elasticities

    use tensorproducts

    implicit none 

    integer, parameter, private       :: dp = 8 ! Default precision
    real(kind=dp), parameter, private :: identity(3,3)  = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    
    real(kind=dp), parameter :: lam_Bennet68  = 7.15d9
    real(kind=dp), parameter :: mu_Bennet68   = 3.455d9
    real(kind=dp), parameter :: eta1_Bennet68 = 5.3e9
    real(kind=dp), parameter :: eta2_Bennet68 = -1.27e9
    real(kind=dp), parameter :: eta3_Bennet68 = -3.95e8    
    
contains      

    !---------------------------------
    ! TRANSVERSELY ISOTROPIC ELASTICITY
    !---------------------------------

    function stress_of_strain__tranisotropic(strain, lam,mu, Elam,Emu,Egam,m) result(stress)

        implicit none
        real(kind=dp), intent(in) :: strain(3,3), lam,mu, Elam,Emu,Egam,m(3)
        real(kind=dp)             :: stress(3,3)
        real(kind=dp)             :: mm(3,3),L(3,3), J1,J4
        real(kind=dp)             :: k1,k2,k3,k4,k5

        call elastic_tranisotropic_params__stress(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5)                
        call elastic_tranisotropic_tensors_and_invars(strain,m, mm,L,J1,J4)
        stress = k1*J1*identity + k2*strain + k3*(J4*identity + J1*mm) + k4*J4*mm + k5*L
    end

    function strain_of_stress__tranisotropic(stress, lam,mu, Elam,Emu,Egam,m) result(strain)

        implicit none
        real(kind=dp), intent(in) :: stress(3,3), lam,mu, Elam,Emu,Egam,m(3)
        real(kind=dp)             :: strain(3,3)
        real(kind=dp)             :: mm(3,3),L(3,3), I1,I4
        real(kind=dp)             :: k1,k2,k3,k4,k5
        
        call elastic_tranisotropic_params__strain(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5)
        call elastic_tranisotropic_tensors_and_invars(stress,m, mm,L,I1,I4)
        strain = k1*I1*identity + k2*stress + k3*(I4*identity + I1*mm) + k4*I4*mm + k5*L 
    end

    !--- AUX ---

    subroutine elastic_tranisotropic_tensors_and_invars(X,m, mm,L,I1,I4)

        implicit none
        real(kind=dp), intent(in)  :: X(3,3), m(3)
        real(kind=dp), intent(out) :: mm(3,3),L(3,3), I1,I4

        mm = outerprod(m,m)
        L  = matmul(X,mm) + matmul(mm,X)
        I1 = X(1,1) + X(2,2) + X(3,3)
        I4 = doubleinner22(X,mm)
    end
    
    subroutine elastic_tranisotropic_params__stress(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5)

        implicit none
        real(kind=dp), intent(in)  :: lam,mu, Elam,Emu,Egam
        real(kind=dp), intent(out) :: k1,k2,k3,k4,k5
        real(kind=dp)              :: gam
        
        gam = lam + 2*mu ! Isotropic P-wave modulus (gamma)

        k1 = lam
        k2 = 2*mu
        k3 = - lam*(1-Elam)
        k4 = ((1+Egam)*gam - 2*(Elam*lam + 2*Emu*mu))
        k5 = - 2*mu*(1-Emu)
    end
    
    subroutine elastic_tranisotropic_params__strain(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5)

        implicit none
        real(kind=dp), intent(in)  :: lam,mu, Elam,Emu,Egam
        real(kind=dp), intent(out) :: k1,k2,k3,k4,k5
        real(kind=dp)              :: gam, xi

        gam = lam + 2*mu ! Isotropic P-wave modulus
        xi = 2*Egam*gam*(lam+mu) - 2*(Elam*lam)**2

        k1 = ((Elam*lam)**2 - Egam*gam*lam)/(2*mu*xi)
        k2 = 1/(2*mu)
        k3 = (Egam*gam*lam - Elam*lam*(Elam*lam + 2*mu))/(2*mu*xi)
        k4 = (4*mu*lam*(1+Elam) + 4*mu**2 + Egam*gam**2 - Elam**2*lam**2)/(2*mu*xi) - 1/(Emu*mu)
        k5 = (1-Emu)/(2*Emu*mu)
    end

    subroutine elastic_tranisotropic__etai_to_Ei(lam,mu,eta1,eta2,eta3, Elam,Emu,Egam)

        implicit none
        real(kind=dp), intent(in)  :: lam,mu, eta1,eta2,eta3
        real(kind=dp), intent(out) :: Elam,Emu,Egam

        Elam = (lam+eta2)/lam
        Emu  = (mu+eta3)/mu
        Egam = (lam+2*mu + eta1 + 2*eta2 + 4*eta3)/(lam+2*mu)
    end


    subroutine elastic_tranisotropic__Ei_to_etai(lam,mu,Elam,Emu,Egam, eta1,eta2,eta3)

        implicit none
        real(kind=dp), intent(in)  :: lam,mu, Elam,Emu,Egam
        real(kind=dp), intent(out) :: eta1,eta2,eta3

        eta1 = lam + 2*mu + Egam*(lam+2*mu) - 2*Elam*lam -4*Emu*mu
        eta2 = -(1-Elam)*lam
        eta3 = -(1-Emu)*mu
    end

end module elasticities
