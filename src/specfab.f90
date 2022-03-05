! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2022

module specfab  

    ! Top level module that includes all functionalities (all submodules).

    use tensorproducts
    use mandel
    use moments
    use homogenizations
    use enhancementfactors
    use dynamics
    use rheologies
!    use wavepropagation ! in development
    use reducedform

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    integer, private :: ll, mm ! Loop indices
    
    ! Distribution expansion series: n(theta,phi) = sum_{l,m}^{Lcap,:} n_l^m Y_l^m(theta,phi) 
    ! where "nlm" vector := n_l^m = (n_0^0, n_2^-2, n_2^-1, n_2^0, n_2^1, n_2^2, n_4^-4, ... ) 
    integer, private   :: Lcap    ! Truncation "L" of expansion series (internal copy of what was passed to the init routine).
    integer            :: nlm_len ! Total number of expansion coefficients (i.e. DOFs)
    integer, parameter :: nlm_lenvec(0:Lcap__max) = [((ll+1)*(ll+2)/2, ll=0, Lcap__max, 1)] ! nlm length for a given Lcap: sum_{l=0}^L (2*l+1) **for even l** = sum_{l=0}^{L/2} (4*l+1) = (L+1)*(L+2)/2
    
    ! (l,m) vector, etc. (nlm_lenmax set by dynamics.f90)
    integer, parameter :: lm(2,nlm_lenmax) = reshape([( (ll,mm, mm=-ll,ll), ll=0,  Lcap__max,2)], [2,nlm_lenmax]) ! These are the (l,m) values corresponding to the coefficients in "nlm".

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
        call initdynamics(Lcap) 
        call initenhancementfactors(Lcap)
        
    end

    !---------------------------------
    ! EXTERNAL, ADD-ON FAETURES
    !---------------------------------

    ! Elmer/ice flow model (Lilien)
    include "elmer/specfab_elmer.f90"

    ! JOSEF ice flow model (Rathmann and Lilien, 2021)
    include "josef/specfab_josef.f90"

end module specfab 
