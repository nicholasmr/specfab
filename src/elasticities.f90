! N. M. Rathmann <rathmann@nbi.ku.dk>, 2022-2023

! Solid elasticities (constitutive equations)

module elasticities

    use tensorproducts

    implicit none 

    integer, parameter, private       :: dp = 8 ! Default precision
    real(kind=dp), parameter, private :: identity(3,3)  = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    
    ! Bennett (1968) parameters
    ! @TODO: Delete because unused? 
    real(kind=dp), parameter :: lam_B68  = 7.15d9
    real(kind=dp), parameter :: mu_B68   = 3.455d9
    real(kind=dp), parameter :: Elam_B68 = 0.8223776223776224
    real(kind=dp), parameter :: Emu_B68  = 0.8856729377713459
    real(kind=dp), parameter :: Egam_B68 = 1.0839260312944523
    
contains      

    !---------------------------------
    ! TRANSVERSELY ISOTROPIC ELASTICITY
    !---------------------------------

    function elas_rev_tranisotropic(strain, lam,mu, Elam,Emu,Egam,m) result(stress)

        implicit none
        real(kind=dp), intent(in) :: strain(3,3), lam,mu, Elam,Emu,Egam,m(3)
        real(kind=dp)             :: stress(3,3)
        real(kind=dp)             :: mm(3,3),L(3,3), J1,J4
        real(kind=dp)             :: k1,k2,k3,k4,k5

        call elas_revparams_tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5)                
        call elas_structs_tranisotropic(strain,m, mm,L,J1,J4)
        stress = k1*J1*identity + k2*strain + k3*(J4*identity + J1*mm) + k4*J4*mm + k5*L
    end

    function elas_fwd_tranisotropic(stress, lam,mu, Elam,Emu,Egam,m) result(strain)

        implicit none
        real(kind=dp), intent(in) :: stress(3,3), lam,mu, Elam,Emu,Egam,m(3)
        real(kind=dp)             :: strain(3,3)
        real(kind=dp)             :: mm(3,3),L(3,3), I1,I4
        real(kind=dp)             :: k1,k2,k3,k4,k5
        
        call elas_fwdparams_tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5)
        call elas_structs_tranisotropic(stress,m, mm,L,I1,I4)
        strain = k1*I1*identity + k2*stress + k3*(I4*identity + I1*mm) + k4*I4*mm + k5*L 
    end

    subroutine elas_structs_tranisotropic(X,m, mm,L,I1,I4)

        implicit none
        real(kind=dp), intent(in)  :: X(3,3), m(3)
        real(kind=dp), intent(out) :: mm(3,3),L(3,3), I1,I4

        mm = outerprod(m,m)
        L  = matmul(X,mm) + matmul(mm,X)
        I1 = X(1,1) + X(2,2) + X(3,3)
        I4 = doubleinner22(X,mm)
    end
    
    subroutine elas_revparams_tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5)

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
    
    subroutine elas_fwdparams_tranisotropic(lam,mu,Elam,Emu,Egam, k1,k2,k3,k4,k5)

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

    !---------------------------------
    ! ORTHOTROPIC ELASTICITY
    !---------------------------------
   
    ! TO BE IMPLEMENTED
   
    !---------------------------------
    ! AUX ROUTINES
    !--------------------------------- 
    
    subroutine Cij_to_Lame_tranisotropic(Cij, Lame)

        implicit none
        real(kind=dp), intent(in)  :: Cij(5)  ! C11,C33,C55,C12,C13
        real(kind=dp), intent(out) :: Lame(5) ! lam,mu, Elam,Emu,Egam
        real(kind=dp)              :: C66

        C66 = (Cij(1)-Cij(4))/2
        Lame(1) = Cij(4)
        Lame(2) = C66
        Lame(3) = Cij(5)/Cij(4)
        Lame(4) = Cij(3)/C66
        Lame(5) = Cij(2)/Cij(1)
    end
    
    subroutine Cij_to_Lame_orthotropic(Cij, Lame)

        implicit none
        real(kind=dp), intent(in)  :: Cij(9)  ! C11,C22,C33,C44,C55,C66, C23,C13,C12
        real(kind=dp), intent(out) :: Lame(9) ! lam11,lam22,lam33, lam23,lam13,lam12, mu1,mu2,mu3

        Lame(7) = Cij(5)+Cij(6)-Cij(4) ! mu1
        Lame(8) = Cij(4)+Cij(6)-Cij(5) ! mu2
        Lame(9) = Cij(4)+Cij(5)-Cij(6) ! mu3
        
        Lame(1) = Cij(1)-2*Lame(7)
        Lame(2) = Cij(2)-2*Lame(8)
        Lame(3) = Cij(3)-2*Lame(9)
        
        Lame(4) = Cij(7)
        Lame(5) = Cij(8)
        Lame(6) = Cij(9)
    end
    
end module elasticities
