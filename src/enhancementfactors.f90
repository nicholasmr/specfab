! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2023

! Stran-rate enhancement factors assuming Sachs and Taylor homogenizations.

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

    ! Transversely isotropic grains
    
    function Eij_tranisotropic(nlm, e1,e2,e3, Eij_grain,alpha,n_grain) result (Eij)

        ! Enhancement factors w.r.t. (e1,e2,e3) axes.
        ! Note: if (e1,e2,e3) coincide with fabric symmetry axes, then these are the fabric eigenenhancements.

        implicit none

        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), dimension(3)  :: e1,e2,e3
        real(kind=dp), intent(in)    :: Eij_grain(2), alpha
        integer, intent(in)          :: n_grain
        real(kind=dp)                :: Eij(6)
        
        ! Longitudinal
        Eij(1) = Evw_tranisotropic(e1,e1, tau_vv(e1),    nlm, Eij_grain,alpha,n_grain) 
        Eij(2) = Evw_tranisotropic(e2,e2, tau_vv(e2),    nlm, Eij_grain,alpha,n_grain)
        Eij(3) = Evw_tranisotropic(e3,e3, tau_vv(e3),    nlm, Eij_grain,alpha,n_grain)    

        ! Shear
        Eij(4) = Evw_tranisotropic(e2,e3, tau_vw(e2,e3), nlm, Eij_grain,alpha,n_grain)
        Eij(5) = Evw_tranisotropic(e1,e3, tau_vw(e1,e3), nlm, Eij_grain,alpha,n_grain) 
        Eij(6) = Evw_tranisotropic(e1,e2, tau_vw(e1,e2), nlm, Eij_grain,alpha,n_grain) 
    end
    
    function Evw_tranisotropic(v,w,tau, nlm, Eij_grain,alpha,n_grain) result(Evw)

        ! Generalized enhancement factor for transversely isotropic grains and a linear Taylor--Sachs homogenization scheme.

        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(in)    :: Eij_grain(2), alpha, v(3),w(3), tau(3,3)
        integer, intent(in)          :: n_grain
        real(kind=dp)                :: vw(3,3), Evw_sachs, Evw_taylor, Evw

        vw = outerprod(v,w)    
        Evw_sachs  = 0.0d0
        Evw_taylor = 0.0d0
    
        Evw_sachs = doubleinner22(rheo_fwd_tranisotropic_sachshomo(           tau, nlm, Eij_grain,n_grain), vw) / &
                    doubleinner22(rheo_fwd_tranisotropic_sachshomo__isotropic(tau,      Eij_grain,n_grain), vw)

        Evw_taylor = doubleinner22(rheo_fwd_tranisotropic_taylorhomo(           tau, nlm, Eij_grain,n_grain), vw) / &
                     doubleinner22(rheo_fwd_tranisotropic_taylorhomo__isotropic(tau,      Eij_grain,n_grain), vw)

        Evw = (1-alpha)*Evw_Sachs + alpha*Evw_taylor       
    end

    ! Orthotropic grains
    
    function Eij_orthotropic(nlm_r1,nlm_r2,nlm_r3, e1,e2,e3, Eij_grain,alpha,n_grain)  result(Eij)

        ! Enhancement factors w.r.t. (e1,e2,e3) axes.
        ! Note: if (e1,e2,e3) coincide with fabric symmetry axes, then these are the fabric eigenenhancements.

        implicit none

        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
        real(kind=dp), dimension(3)  :: e1,e2,e3
        real(kind=dp), intent(in)    :: Eij_grain(6), alpha
        integer, intent(in)          :: n_grain
        real(kind=dp)                :: Eij(6)
        
        ! Longitudinal
        Eij(1) = Evw_orthotropic(e1,e1, tau_vv(e1),    nlm_r1,nlm_r2,nlm_r3, Eij_grain,alpha,n_grain) 
        Eij(2) = Evw_orthotropic(e2,e2, tau_vv(e2),    nlm_r1,nlm_r2,nlm_r3, Eij_grain,alpha,n_grain)
        Eij(3) = Evw_orthotropic(e3,e3, tau_vv(e3),    nlm_r1,nlm_r2,nlm_r3, Eij_grain,alpha,n_grain)    

        ! Shear
        Eij(4) = Evw_orthotropic(e2,e3, tau_vw(e2,e3), nlm_r1,nlm_r2,nlm_r3, Eij_grain,alpha,n_grain)
        Eij(5) = Evw_orthotropic(e1,e3, tau_vw(e1,e3), nlm_r1,nlm_r2,nlm_r3, Eij_grain,alpha,n_grain) 
        Eij(6) = Evw_orthotropic(e1,e2, tau_vw(e1,e2), nlm_r1,nlm_r2,nlm_r3, Eij_grain,alpha,n_grain) 
    end
    
    function Evw_orthotropic(v,w,tau, nlm_r1,nlm_r2,nlm_r3, Eij_grain,alpha,n_grain)  result(Evw)

        ! Generalized enhancement factor for orthotropic grains and a linear Taylor--Sachs homogenization scheme.

        implicit none
        
        complex(kind=dp), intent(in) :: nlm_r1(:), nlm_r2(:), nlm_r3(:)
        real(kind=dp), intent(in)    :: Eij_grain(6), alpha, v(3),w(3), tau(3,3)
        integer, intent(in)          :: n_grain
        real(kind=dp)                :: vw(3,3), Evw_sachs, Evw_taylor, Evw
        complex(kind=dp)             :: nlm_iso(size(nlm_r1))
    
        vw = outerprod(v,w)
        Evw_sachs  = 0.0d0
        Evw_taylor = 0.0d0
    
        nlm_iso(:) = 0.0d0
        nlm_iso(1) = 1.0d0
    
        Evw_sachs = doubleinner22(rheo_fwd_orthotropic_sachshomo(tau, nlm_r1,  nlm_r2,  nlm_r3,  Eij_grain,n_grain), vw) / &
                    doubleinner22(rheo_fwd_orthotropic_sachshomo(tau, nlm_iso, nlm_iso, nlm_iso, Eij_grain,n_grain), vw)
                    
        Evw_taylor = doubleinner22(rheo_fwd_orthotropic_taylorhomo(tau, nlm_r1,  nlm_r2,  nlm_r3,  Eij_grain,n_grain), vw) / &
                     doubleinner22(rheo_fwd_orthotropic_taylorhomo(tau, nlm_iso, nlm_iso, nlm_iso, Eij_grain,n_grain), vw)
                    
        Evw = (1-alpha)*Evw_Sachs + alpha*Evw_taylor                    
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
