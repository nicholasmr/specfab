! N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

! Constants etc.

module header  

    implicit none 

    !------------------

    integer, parameter       :: dp = 8 ! Default precision
    real(kind=dp), parameter :: Pi = 3.141592653589793d0
    !integer, parameter       :: identity(3,3)  = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    real(kind=dp), parameter :: identity(3,3) = reshape([1.0d0,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,0.0d0,1.0d0], [3,3]) ! more safe
    
    integer :: ii,jj,kk,ll,mm,nn,oo,pp ! loop indices
    
    !------------------
    
    integer, parameter :: Lcap__max  = 30 ! Hard limit
    
    ! Length of nlm state vector for a given L
    integer, parameter :: nlm_lenvec(0:Lcap__max) = [((ll+1)*(ll+2)/2, ll=0, Lcap__max, 1)] ! nlm length for a given L: sum_{l=0}^L (2*l+1) **for even l** = sum_{l=0}^{L/2} (4*l+1) = (L+1)*(L+2)/2
    integer, parameter :: nlm_lenmax = (Lcap__max+1)*(Lcap__max+2)/2
    
    !... shorthands
    integer, parameter :: L2len = nlm_lenvec(2)
    integer, parameter :: L4len = nlm_lenvec(4)
    integer, parameter :: L6len = nlm_lenvec(6)
    integer, parameter :: L8len = nlm_lenvec(8)
        
    ! Degree and order (l,m) vector
    integer, parameter :: lm(2,nlm_lenmax) = reshape([( (ll,mm, mm=-ll,ll), ll=0,  Lcap__max,2)], [2,nlm_lenmax]) ! These are the SH degree and order integers (l,m) corresponding to the coefficients in "nlm".

    ! Indices for extracting l=0,2,4,6,8 coefs of nlm
    integer, parameter :: IL(0:Lcap__max)  = [( nlm_lenvec(ll)-(2*ll+1)+1, ll=0, Lcap__max, 1)] ! Index into nlm for n_l^-l coef. given l; i.e. first coef for l=L.
    integer, parameter :: ILm(0:Lcap__max) = [( IL(ll)+ll, ll=0, Lcap__max, 1)] ! Index into nlm for n_l^0 coef. given l

    integer, parameter :: L0rng      = 1
    integer, parameter :: L2rng(5)   = [(ii, ii=IL(2), IL(4)-1, 1)] ! nlm(L2rng) = [n_2^-2, n_2^-1, n_2^0, n_2^1, n_2^2]
    integer, parameter :: L4rng(9)   = [(ii, ii=IL(4), IL(6)-1, 1)]
    integer, parameter :: L6rng(13)  = [(ii, ii=IL(6), IL(8)-1, 1)]
    integer, parameter :: L8rng(17)  = [(ii, ii=IL(8), IL(10)-1, 1)]
    integer, parameter :: L10rng(21) = [(ii, ii=IL(10), IL(12)-1, 1)]
    integer, parameter :: L12rng(25) = [(ii, ii=IL(12), IL(14)-1, 1)]

    !------------------

contains   

end module header
