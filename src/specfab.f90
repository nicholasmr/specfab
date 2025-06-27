! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2022

module specfab  

    ! Top level module that includes all functionalities (all submodules).

    use header
    use tensorproducts
    use mandel
    use moments
    use homogenizations
    use enhancementfactors
    use dynamics
    use rheologies
    use wavepropagation 
    use rotation
    use reducedform
    use deformationmodes
    use idealstate
    use frames

    implicit none 

    ! Distribution expansion series: n(theta,phi) = sum_{l,m}^{Lcap,:} n_l^m Y_l^m(theta,phi) 
    ! where "nlm" vector := n_l^m = (n_0^0, n_2^-2, n_2^-1, n_2^0, n_2^1, n_2^2, n_4^-4, ... ) 
    integer, private   :: Lcap    ! Truncation "L" of expansion series (internal copy of what was passed to the init routine).
    integer            :: nlm_len ! Total number of expansion coefficients (i.e. DOFs)

contains      

    !---------------------------------
    ! INIT
    !---------------------------------
           
    subroutine initspecfab(Lcap_)

        ! Must be called once before using module

        implicit none    
        integer, intent(in) :: Lcap_ ! Truncation "Lcap"
        
        Lcap = Lcap_ ! Save internal copy
        nlm_len = nlm_lenvec(Lcap) ! Number of DOFs (expansion coefficients) for full nlm vector

        call initreduced(Lcap) ! Initialize reduced nlm (rnlm) module for 2D x-z problems
        call initmoments(Lcap)
        call initdynamics(Lcap)
        call inithomogenizations(Lcap) ! Set isotropic structure tensors
        
    end

    !---------------------------------
    ! EXTERNAL ADD-ONS
    !---------------------------------

    ! Elmer/ice flow model (Lilien)
    include "elmer/specfab_elmer.f90"

end module specfab 
