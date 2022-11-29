! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2022

! Stran-rate enhancement factors relying on the Sachs and Taylor homogenization assumptions.

module enhancementfactors  

    use tensorproducts
    use moments
    use homogenizations

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    
contains      

    !---------------------------------
    ! VISCOUS ENHANCEMENT-FACTORS
    !---------------------------------

    function Evw(vw, tau, nlm, Ecc,Eca,alpha,nprime)

        ! Generalized enhancement factor, E_{vw}
        ! Assumes a transversely isotropic grain rheology (Ecc,Eca,alpha,nprime)

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: Ecc, Eca, alpha, vw(3,3), tau(3,3)
        integer, intent(in)          :: nprime
        real(kind=dp)                :: Evw

        Evw = (1-alpha)*Evw_sachs( vw, tau, nlm, Ecc,Eca,nprime) &
                + alpha*Evw_taylor(vw, tau, nlm, Ecc,Eca,nprime)
    end

    function Eeiej(nlm, e1,e2,e3, Ecc,Eca,alpha,nprime) result (Eij)

        ! Enhancement factors in directions (ei,ej), *not* necessarily a2 eigen directions.
        ! Returns a 3x3 symmetric matrix of enhancement factors.

        implicit none

        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), dimension(3)  :: e1,e2,e3
        real(kind=dp), intent(in)    :: Ecc, Eca, alpha
        integer, intent(in)          :: nprime
        real(kind=dp)                :: Eij(3,3)
        
        ! Longitudinal
        Eij(1,1) = Evw(outerprod(e1,e1), tau_vv(e1),     nlm, Ecc,Eca,alpha,nprime) 
        Eij(2,2) = Evw(outerprod(e2,e2), tau_vv(e2),     nlm, Ecc,Eca,alpha,nprime)
        Eij(3,3) = Evw(outerprod(e3,e3), tau_vv(e3),     nlm, Ecc,Eca,alpha,nprime)    
        
        ! Shear
        Eij(1,2) = Evw(outerprod(e1,e2), tau_vw(e1,e2),  nlm, Ecc,Eca,alpha,nprime) 
        Eij(1,3) = Evw(outerprod(e1,e3), tau_vw(e1,e3),  nlm, Ecc,Eca,alpha,nprime) 
        Eij(2,3) = Evw(outerprod(e2,e3), tau_vw(e2,e3),  nlm, Ecc,Eca,alpha,nprime)
        
        ! Symmetric matrix
        Eij(2,1) = Eij(1,2)
        Eij(3,1) = Eij(1,3) 
        Eij(3,2) = Eij(2,3)   
    end

    function Evw_sachs(vw, tau, nlm, Ecc,Eca, nprime) result (Evw)

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: Ecc, Eca, vw(3,3), tau(3,3)
        integer, intent(in)          :: nprime
        real(kind=dp)                :: Evw
        
        Evw = doubleinner22(rheo_fwd_tranisotropic_sachshomo(           tau, nlm, Ecc,Eca,nprime), vw) / &
              doubleinner22(rheo_fwd_tranisotropic_sachshomo__isotropic(tau,      Ecc,Eca,nprime), vw)
    end

    function Evw_taylor(vw, tau, nlm, Ecc,Eca, nprime) result (Evw)

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: Ecc, Eca, vw(3,3), tau(3,3)
        integer, intent(in)          :: nprime
        real(kind=dp)                :: Evw

        Evw = doubleinner22(rheo_fwd_tranisotropic_taylorhomo(           tau, nlm, Ecc,Eca,nprime), vw) / &
              doubleinner22(rheo_fwd_tranisotropic_taylorhomo__isotropic(tau,      Ecc,Eca,nprime), vw)
    end
    
    !---------------------------------
    ! ELASTIC ENHANCEMENT-FACTORS
    !---------------------------------

    ! N/A

    !---------------------------------
    ! SYNTHETIC STRESS STATES
    !---------------------------------

    function tau_vv(v) 
        ! v--v compression/extension
        implicit none
        real(kind=dp), intent(in) :: v(3)
        integer, parameter :: identity(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
        real(kind=dp) :: tau_vv(3,3)
        tau_vv = identity/3.0d0 - outerprod(v,v)
    end

    function tau_vw(v,w)
        ! v--w shear
        implicit none
        real(kind=dp), intent(in) :: v(3), w(3)
        real(kind=dp) :: tau_vw(3,3)
        tau_vw = outerprod(v,w) + outerprod(w,v)
    end

end module enhancementfactors 
