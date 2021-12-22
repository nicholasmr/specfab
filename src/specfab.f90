! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2019-2022

module specfab  

    use tensorproducts
    use moments
    use gaunt
    use rheologies

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real, parameter, private    :: Pi = 3.1415927
    integer, parameter, private :: x = 1, y = 2, z = 3 ! Matrix indices
    complex(kind=dp), parameter, private :: r = (1,0), i = (0,1) ! real and imag units
    integer, private :: ii, jj, kk, ll, mm ! Loop indicies

    ! Expansion series: n(theta,phi) = sum_{l,m}^{Lcap,:} n_l^m Y_l^m(theta,phi) 
    ! where "nlm" vector := n_l^m = (n_0^0, n_2^-2, n_2^-1, n_2^0, n_2^1, n_2^2, n_4^-4, ... ) 
    integer, private   :: Lcap    ! Truncation "L" of expansion series (internal copy of what was passed to the init routine).
    integer            :: nlm_len             ! Total number of expansion coefficients (i.e. DOFs)
    integer            :: nlm_reduced_len     ! ...same but for "reduced" nlm vector (see below comment on reduced form since ODF is real-valued).
    integer            :: nlm_reduced_neg_len ! ...same but for negative part of "reduced" nlm vector
    integer, parameter :: Lcap__max                   = 60 ! Hard limit
    integer, parameter :: nlm_len_from_L(0:Lcap__max) = [((ll+1)*(ll+2)/2, ll=0, Lcap__max, 1)] ! nlm length for a given Lcap: sum_{l=0}^L (2*l+1) **for even l** = sum_{l=0}^{L/2} (4*l+1) = (L+1)*(L+2)/2
    integer, parameter :: nlm_len__max                = nlm_len_from_L(Lcap__max)
    
    ! For converting between full (nlm) and reduced (nlm_reduced) forms of nlm by noting n_l^{-m} = (-1)^m conj(n_l^m) for real-valued ODFs.
    integer, parameter :: nlm_reduced_len_from_L(0:Lcap__max)     = [( (ll+2)**2/4,     ll=0, Lcap__max, 1)] ! nlm length for a given Lcap: sum_{l=0}^L (1*l+1) **for even l** = sum_{l=0}^{L/2} (2*l+1) = (L+2)^2/4
    integer, parameter :: nlm_reduced_neg_len_from_L(0:Lcap__max) = [( ll*(ll+2)/4,     ll=0, Lcap__max, 1)] ! nlm length for a given Lcap: sum_{l=0}^L l       **for even l** = sum_{l=0}^{L/2} (2*l)   = L*(L+2)/4
    integer, parameter :: nlm_reduced_len__max                    = nlm_reduced_len_from_L(Lcap__max)
    integer, parameter :: nlm_reduced_neg_len__max                = nlm_reduced_neg_len_from_L(Lcap__max)
    integer, parameter :: I_full(nlm_len__max)                    = [( (kk+nlm_reduced_len_from_L(ll)-1*ll, kk=ll, 1, -1), (kk+nlm_reduced_len_from_L(ll)-1*ll, kk=0, ll, 1), ll=0, Lcap__max, 2)] ! nlm(1:nlm_len) = Mdiag(1:nlm_len) * nlm_reduced(I_full(1:nlm_len(L))) ...and then nlm(I_reduced_neg(1:nlm_reduced_neg_len(L))) = conj( nlm(I_reduced_neg(1:nlm_reduced_neg_len(L))) )
    integer, parameter :: I_reduced(nlm_reduced_len__max)         = [( (kk+nlm_len_from_L(ll)-ll , kk=0, ll, 1), ll=0, Lcap__max, 2)]     ! nlm_reduced     = nlm(I_reduced(    1:nlm_reduced_len(L))) = (n_0^0, n_2^0, n_2^1, n_2^2, n_4^0, ... ) 
    integer, parameter :: I_reduced_neg(nlm_reduced_neg_len__max) = [( (kk+nlm_len_from_L(ll)-2*ll , kk=0, ll-1, 1), ll=0, Lcap__max, 2)] ! nlm_reduced_neg = nlm(I_reduced_neg(1:nlm_reduced_len(L))) = (n_2^-2, n_2^-1, n_4^-4, n_4^-3, n_4^-2, n_4^-1, n_6^-6, ... ) 
    integer, parameter :: Mdiag(nlm_len__max)                     = [( ((-1)**(mm), mm=ll, 1, -1), (1, kk=0, ll, 1), ll=0, Lcap__max, 2)] ! For adjusting odd negative m values with a factor of -1 (i.e. the (-1)^m factor in n_l^{-m} = (-1)^m conj(n_l^m))

    ! (l,m) vector, etc.
    integer, parameter :: lm(2,nlm_len__max) = reshape([( (ll,mm, mm=-ll,ll), ll=0,  Lcap__max,2)], [2,nlm_len__max]) ! These are the (l,m) values corresponding to the coefficients in "nlm".
    integer, parameter :: I_l0=1, I_l2=I_l0+1, I_l4=I_l2+(2*2+1), I_l6=I_l4+(2*4+1), I_l8=I_l6+(2*6+1), I_l10=I_l8+(2*8+1) ! Indices for extracting l=0,2,4,6,8 coefs of nlm
    
    ! Dynamics
    complex(kind=dp), private, allocatable :: regmat(:,:) ! Unscaled regularization matrix
    complex(kind=dp), private, allocatable :: lapmat(:,:) ! Laplacian diffusion operator

    ! Ehancement-factor related
    real(kind=dp), private :: ev_c2_iso(3,3), ev_c4_iso(3,3,3,3), ev_c6_iso(3,3,3,3, 3,3), ev_c8_iso(3,3,3,3, 3,3,3,3) ! <c^k> for isotropic n(theta,phi)
    real(kind=dp), parameter, private :: identity(3,3)  = reshape([1,0,0, &
                                                                   0,1,0, &
                                                                   0,0,1], [3,3])
    real(kind=dp), parameter, private :: identity9(9,9) = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, &
                                                                   0, 1, 0, 0, 0, 0, 0, 0, 0, & 
                                                                   0, 0, 1, 0, 0, 0, 0, 0, 0, & 
                                                                   0, 0, 0, 1, 0, 0, 0, 0, 0, & 
                                                                   0, 0, 0, 0, 1, 0, 0, 0, 0, & 
                                                                   0, 0, 0, 0, 0, 1, 0, 0, 0, & 
                                                                   0, 0, 0, 0, 0, 0, 1, 0, 0, &
                                                                   0, 0, 0, 0, 0, 0, 0, 1, 0, &
                                                                   0, 0, 0, 0, 0, 0, 0, 0, 1], [9,9])
    
    ! Optimal n'=1 (lin) grain parameters 
    ! These are the linear mixed Taylor--Sachs best-fit parameters from Rathmann and Lilien (2021)
    real(kind=dp), parameter :: Eca_opt_lin   = 1d3
    real(kind=dp), parameter :: Ecc_opt_lin   = 1d0
    real(kind=dp), parameter :: alpha_opt_lin = 0.0125
    
    ! Optimal n'=3 (nlin) grain parameters 
    ! These are the nonlinear Sachs-only best-fit parameters (Rathmann et al., 2021) 
    real(kind=dp), parameter :: Eca_opt_nlin   = 1d4
    real(kind=dp), parameter :: Ecc_opt_nlin   = 1d0
    real(kind=dp), parameter :: alpha_opt_nlin = 0
    
contains      

!---------------------------------
! INIT
!---------------------------------
       
subroutine initspecfab(Lcap_)

    ! Needs to be called once before using the module routines.

    implicit none    
    integer, intent(in) :: Lcap_ ! Trauncation "Lcap"
    complex(kind=dp)    :: nlm_iso(nlm_len__max) = 0.0
    
    Lcap = Lcap_ ! Save internal copy
    
    nlm_len             = nlm_len_from_L(Lcap)             ! Number of DOFs (expansion coefficients) for full nlm vector
    nlm_reduced_len     = nlm_reduced_len_from_L(Lcap)     ! ...same but for reduced nlm vector
    nlm_reduced_neg_len = nlm_reduced_neg_len_from_L(Lcap) ! ...same but for negative "m" part of reduced nlm vector

    ! Set gaunt coefficients (overlap integrals involving three spherical harmonics)
    call set_gaunts()
    
    ! Calculate structure tensors for an isotropic fabric (used for calculating enhancement factors)
    nlm_iso(1) = 1/Sqrt(4*Pi) ! Normalized ODF 
    call f_ev_ck(nlm_iso, 'f', ev_c2_iso,ev_c4_iso,ev_c6_iso,ev_c8_iso) ! Sets <c^i> := a^(i) for i=2,4,6,8
    
    ! Static (constant) matrices used for spectral dynamics
    if(.not. allocated(regmat)) allocate(regmat(nlm_len,nlm_len)) ! Regularization (unscaled)
    if(.not. allocated(lapmat)) allocate(lapmat(nlm_len,nlm_len)) ! Laplacian operator
    regmat = 0.0 ! initialize
    lapmat = 0.0 ! initialize
    do ii = 1, nlm_len 
        regmat(ii,ii) = -(lm(1,ii)*(lm(1,ii)+1))**(1.5) ! if expo > 1 => hyper diffusion
        lapmat(ii,ii) = -(lm(1,ii)*(lm(1,ii)+1))
    end do    
end

!---------------------------------
! FABRIC DYNAMICS
!---------------------------------

function dndt_ij_LATROT(eps,omg, tau,Aprime,Ecc,Eca, beta)  

    ! *** LATTICE ROTATION ***

    ! Returns matrix "dndt_ij" such that
    ! 
    !       d/dt (nlm)_i = dndt_ij (nlm)_j
    ! 
    ! where "nlm" is the vector of spectral coefficient n_l^m

    implicit none

    real(kind=dp), intent(in) :: eps(3,3), omg(3,3), tau(3,3) ! strain-rate (eps), spin (omg), dev. stress (tau)
    real(kind=dp), intent(in) :: Aprime, Ecc, Eca, beta
    complex(kind=dp)          :: dndt_ij_LATROT(nlm_len,nlm_len), qe(-2:2), qt(-2:2), qo(-1:1)
    integer, parameter        :: SHI_LATROT = 6 ! Scope of harmonic interactions (in wave space) for LATROT
    complex(kind=dp), dimension(SHI_LATROT) :: g0,     gz,     gn,     gp
    complex(kind=dp), dimension(SHI_LATROT) :: g0_rot, gz_rot, gn_rot, gp_rot
    complex(kind=dp), dimension(SHI_LATROT) :: g0_Sac, gz_Sac, gn_Sac, gp_Sac
    complex(kind=dp), dimension(SHI_LATROT) :: g0_Tay, gz_Tay, gn_Tay, gp_Tay
    real(kind=dp) :: etaprime

    etaprime = Aprime*doubleinner22(tau,tau) ! Assumes eta' = A'*(I2(tau)) (i.e. n'=3).

    ! Quadric expansion coefficients
    qe = quad_rr(eps)
    qt = quad_rr(tau)
    qo = quad_tp(omg)

    ! Harmonic interaction weights

    g0_rot = [ 0*r, 0*r,0*r,0*r,0*r,0*r ]
    gz_rot = [ -i*sqrt(3.)*qo(0), 0*r,0*r,0*r,0*r,0*r ]
    gn_rot = [ -i*6/sqrt(6.)*qo(-1), 0*r,0*r,0*r,0*r,0*r ]
    gp_rot = [ +i*6/sqrt(6.)*qo(+1), 0*r,0*r,0*r,0*r,0*r ]

    g0_Tay = 3*[ 0*r, qe(-2),qe(-1),qe(0),qe(+1),qe(+2) ]
    gz_Tay = [ 0*r, -qe(-2),0*r,0*r,0*r,qe(+2) ]
    gn_Tay = [ sqrt(5./6)*qe(-1), 0*r, qe(-2), sqrt(2./3)*qe(-1), sqrt(3./2)*qe(0), 2*qe(+1) ]
    gp_Tay = [ sqrt(5./6)*qe(+1), 2*qe(-1), sqrt(3./2)*qe(0), sqrt(2./3)*qe(+1), qe(+2), 0*r ]
    
    if (beta .lt. 1.0d0 - 1d-5) then
        ! For efficiency, calculate the Sachs contribution only if beta < 1
        g0_Sac = 3*[ 0*r, Eca*qt(-2),Eca*qt(-1),Eca*qt(0),Eca*qt(+1),Eca*qt(+2) ]
        gz_Sac = [ 0*r, -Eca*qt(-2),-1./2*(Eca-1)*qt(-1),0*r,1./2*(Eca-1)*qt(+1),Eca*qt(+2) ]
        gn_Sac = [ sqrt(5./6)*qt(-1), 0*r, Eca*qt(-2), sqrt(1./6)*(3*Eca-1)*qt(-1), sqrt(3./2)*Eca*qt(0), (Eca+1)*qt(+1) ]
        gp_Sac = [ sqrt(5./6)*qt(+1), (Eca+1)*qt(-1), sqrt(3./2)*Eca*qt(0), sqrt(1./6)*(3*Eca-1)*qt(+1), Eca*qt(+2), 0*r ]

        g0 = g0_rot + (1-beta)*etaprime*g0_Sac + beta*g0_Tay
        gz = gz_rot + (1-beta)*etaprime*gz_Sac + beta*gz_Tay
        gn = gn_rot + (1-beta)*etaprime*gn_Sac + beta*gn_Tay
        gp = gp_rot + (1-beta)*etaprime*gp_Sac + beta*gp_Tay
    else
        ! If beta = 1, the ODF evolution is kinematic in the sense that c-axes rotate in response to the bulk stretching and spin.
        g0 = g0_rot + g0_Tay
        gz = gz_rot + gz_Tay
        gn = gn_rot + gn_Tay
        gp = gp_rot + gp_Tay
    end if 

    do ii = 1, nlm_len
        dndt_ij_LATROT(ii,1:nlm_len) = -1*( matmul(GC(ii,1:nlm_len,1:SHI_LATROT),g0) + matmul(GCm(ii,1:nlm_len,1:SHI_LATROT),gz) + matmul(GC_m1(ii,1:nlm_len,1:SHI_LATROT),gn) + matmul(GC_p1(ii,1:nlm_len,1:SHI_LATROT),gp) )    
    end do
    
end

function dndt_ij_DDRX(nlm, tau)

    ! *** Discontinuous dynamic recrystalization (DDRX) ***
    ! 
    ! Nucleation and migration recrystalization modelled as a decay process (Placidi et al., 2010)

    ! Returns matrix "dndt_ij" such that
    ! 
    !       d/dt (nlm)_i = dndt_ij (nlm)_j
    ! 
    ! where "nlm" is the vector of spectral coefficient n_l^m    
    
    ! **NOTICE** This is Gamma/Gamma_0. The caller must multiply by an appropriate DDRX rate factor, Gamma_0(T,tau,eps,...).

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in)    :: tau(3,3) ! Deviatoric stress tensor
    complex(kind=dp)             :: dndt_ij_DDRX(nlm_len,nlm_len)
    real(kind=dp)                :: Davg

    dndt_ij_DDRX = dndt_ij_DDRX_src(tau)

    ! Add add (nonlinear) sink term <D> to all diagonal entries such that dndt_ij = Gamma_ij/Gamma0 = (D_ij - <D>*I_ij)/(tau:tau) (where I is the identity)
    Davg = doubleinner22(matmul(tau,tau), a2(nlm)) - doubleinner22(tau,doubleinner42(a4(nlm),tau)) ! (tau.tau):a2 - tau:a4:tau 
    Davg = Davg/doubleinner22(tau,tau) ! normalize
    do ii = 1, nlm_len    
        dndt_ij_DDRX(ii,ii) = dndt_ij_DDRX(ii,ii) - Davg 
    end do
end

function dndt_ij_DDRX_src(tau)  

    ! Calculates the ODF-independent (linear) source term of DDRX 

    implicit none

    real(kind=dp), intent(in) :: tau(3,3) ! dev. stress tensor
    complex(kind=dp)          :: dndt_ij_DDRX_src(nlm_len,nlm_len), qt(-2:2)
    integer, parameter        :: SHI_DDRX = 1+5+9 ! Scope of harmonic interactions (in wave space) for DDRX
    real(kind=dp)             :: k
    complex(kind=dp)          :: g(SHI_DDRX)

    ! Quadric expansion coefficients
    qt = quad_rr(tau)

    ! Harmonic interaction weights (requires qt, sets g and k)
    include "include/DDRX__body.f90"

    ! D
    g = k*g ! common prefactor
    g = g/doubleinner22(tau,tau) ! normalize (can be done already here since it is a constant prefactor for dndt_ij_DDRX_src)
    do ii = 1, nlm_len    
        dndt_ij_DDRX_src(ii,1:nlm_len) = matmul(GC(ii,:nlm_len,:), g)
    end do
end

function dndt_ij_CDRX()

    ! *** Continuous dynamic recrystalization (CDRX) ***
    ! 
    ! Rotation recrystalization (polygonization) as a Laplacian diffusion process (Godert, 2003)

    ! Returns matrix "dndt_ij" such that
    ! 
    !       d/dt (nlm)_i = dndt_ij (nlm)_j
    ! 
    ! where "nlm" is the vector of spectral coefficient n_l^m    

    implicit none

    complex(kind=dp) :: dndt_ij_CDRX(nlm_len,nlm_len)

    dndt_ij_CDRX = lapmat 
end

!---------------------------------
! FABRIC DYNAMICS IN TENSORIAL SPACE
!---------------------------------

include "tensorialdynamics.f90"

!---------------------------------
! REGULARIZATION 
!---------------------------------

function dndt_ij_REG() 

    ! *** LAPLACIAN REGULARIZATION ***

    ! Returns matrix "dndt_ij" such that
    ! 
    !       d/dt (nlm)_i = dndt_ij (nlm)_j
    ! 
    ! where "nlm" is the vector of spectral coefficient n_l^m   

    implicit none
    complex(kind=dp) :: dndt_ij_REG(nlm_len,nlm_len)

    dndt_ij_REG = regmat ! The caller must multiply by the regularization strength, e.g. f_nu(nu0, eps, ...)
end

function f_nu(nu0, beta, eps, tau, Aprime) result (nu)
    
    ! Regularization strength
    
    implicit none
    
    real(kind=dp), intent(in) :: nu0, beta, eps(3,3),tau(3,3), Aprime
    real(kind=dp)             :: nu, etaprime

    if (beta .le. 1.0d-20) then 
        nu = f_nu_eps(nu0, eps) ! Seperate (faster) calculation if beta=0
    else 
        etaprime = Aprime*doubleinner22(tau,tau) ! Assumes eta' ~= A'*(I2(tau)) (i.e. n'=3).
        nu = beta*f_nu_eps(nu0, eps) + (1-beta)*etaprime*f_nu_tau(nu0, tau)    
    end if

end

function f_nu_eps(nu0, eps) result (nu)
    
    ! Calculates optimal nu as a function of Lcap and norm2(eps).
    
    implicit none
    
    real(kind=dp), intent(in) :: eps(3,3)
    real(kind=dp)             :: nu0, nu
    real(kind=dp), parameter  :: L0=10
    
    ! nu0=5.0e-2 ! Used for DEMO scripts
    ! nu0=0.5e-3 ! Medium loose
    ! nu0=0.1e-3 ! Used for JOSEF

    nu = nu0 * (Lcap/L0)**(-1.1) * norm2(reshape(eps,[size(eps)]))
end

function f_nu_tau(nu0, tau) result (nu)
    
    ! Calculates optimal nu as a function of Lcap and norm2(tau).
    
    implicit none
    
    real(kind=dp), intent(in) :: tau(3,3)
    real(kind=dp)             :: nu0, nu
    real(kind=dp), parameter  :: L0=10
    
    ! nu0=5e1 ! Fitted
    
    nu = nu0 * (Lcap/L0)**(-0.4) * norm2(reshape(tau,[size(tau)]))
end

!---------------------------------
! STRUCTURE TENSORS
!---------------------------------
       
subroutine f_ev_ck(nlm, opt, ev_c2,ev_c4,ev_c6,ev_c8)
    
    ! "ev_ck" are the structure tensors <c^k> := a^(k) for a given n(theta,phi) prescribed in terms of "nlm"
    
    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    character*1, intent(in)      :: opt 
    real(kind=dp), intent(inout) :: ev_c2(3,3),ev_c4(3,3,3,3),ev_c6(3,3,3,3, 3,3),ev_c8(3,3,3,3, 3,3,3,3)
    complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6), n8m(-8:8)
    
    n00 = nlm(1)
    n2m = nlm(I_l2:(I_l4-1))
    n4m = nlm(I_l4:(I_l6-1))
    ev_c2 = f_ev_c2(n00,n2m)
    ev_c4 = f_ev_c4(n00,n2m,n4m)
    
    if (opt == 'f') then
        ! Full calculation?
        n6m = nlm(I_l6:(I_l8-1))
        n8m = nlm(I_l8:(I_l10-1))
        ev_c6 = f_ev_c6(n00,n2m,n4m,n6m)
        ev_c8 = f_ev_c8(n00,n2m,n4m,n6m,n8m)
    else
        ! Reduced calculation?
        ev_c6 = 0.0d0
        ev_c8 = 0.0d0
    end if
end
      
function a2(nlm) 
    ! a^(2) := <c^2> 
    implicit none
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp)                :: a2(3,3)
    complex(kind=dp)             :: n2m(-2:2) 
    n2m = nlm(I_l2:(I_l4-1))
    a2 = f_ev_c2(nlm(1),n2m)
end

function a4(nlm) 
    ! a^(4) := <c^4> 
    implicit none
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp)                :: a4(3,3,3,3)
    complex(kind=dp)             :: n2m(-2:2), n4m(-4:4)
    n2m = nlm(I_l2:(I_l4-1))
    n4m = nlm(I_l4:(I_l6-1))
    a4 = f_ev_c4(nlm(1),n2m,n4m)
end

function a2_to_nlm(a2) result(nlm)
    ! Get n_2^m from a^(2)
    implicit none
    real(kind=dp), intent(in) :: a2(3,3)
    complex(kind=dp)          :: nlm(1+5)
    nlm = 0.0 ! init
    include "include/a2_to_nlm__body.f90"
end

function a4_to_nlm(a2, a4) result(nlm)
    ! Get n_2^m and n_4^m from a^(2) and a^(4)
    implicit none
    real(kind=dp), intent(in) :: a2(3,3), a4(3,3,3,3)
    complex(kind=dp)          :: nlm(1+5+9)
    nlm = 0.0 ! init
    include "include/a4_to_nlm__body.f90"
end

!---------------------------------
! FABRIC FRAME
!---------------------------------

subroutine frame(nlm, ftype, e1,e2,e3, eigvals)

    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    character*1, intent(in)      :: ftype ! 'x','e','p' (cartensian frame, eigen frame, 45deg-rotated eigen frame)
    integer, parameter           :: n = 3
    real(kind=dp), intent(out)   :: e1(n),e2(n),e3(n), eigvals(3)
    real(kind=dp)                :: p(n),q(n)
    ! If eigen frame
    integer            :: inf
    integer, parameter :: l=3*3-1
    real(kind=dp)      :: e_ij(n,n), work(l)
    
    eigvals = 0.0
    
    ! Cartesian frame
    if (ftype == 'x') then
        e1 = [1,0,0] 
        e2 = [0,1,0] 
        e3 = [0,0,1] 
        
    else
        ! Eigen frame    
        e_ij = a2(nlm)
        call dsyev('V','U',n,e_ij,n,eigvals,work,l,inf)
        e1 = e_ij(:,3)
        e2 = e_ij(:,2)
        e3 = e_ij(:,1)
        eigvals = [eigvals(3),eigvals(2),eigvals(1)] ! Largest first
        
        ! Rotated eigen frame
        if (ftype == 'p') then
            p = (e1+e2)/sqrt(2.0) 
            q = (e1-e2)/sqrt(2.0) 
            e1 = p
            e2 = q
            ! cross product
            e3(1) = p(2) * q(3) - p(3) * q(2)
            e3(2) = p(3) * q(1) - p(1) * q(3)
            e3(3) = p(1) * q(2) - p(2) * q(1)
        end if     
           
    end if
end

!---------------------------------
! ENHANCEMENT-FACTORS
!---------------------------------

function Evw(vw, tau, nlm, Ecc, Eca, alpha, nprime)

    ! *** GENERALIZED DIRECTIONAL ENHANCEMENT FACTOR E_{vw} ***

    ! Assumes a transversely isotropic grain rheology (Ecc,Eca,alpha,nprime)

    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in)    :: Ecc, Eca, alpha, vw(3,3), tau(3,3)
    integer, intent(in)          :: nprime
    real(kind=dp)                :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    real(kind=dp)                :: Evw

    if (nprime .eq. 1) then 
        ! Linear grain rheology (n'=1) relies only on <c^2> and <c^4>
        call f_ev_ck(nlm, 'r', ev_c2,ev_c4,ev_c6,ev_c8) ! Calculate structure tensors of orders 2,4 (6,8 are assumed zero for faster evaluation)
    else if (nprime .eq. 3) then
        ! Nonlinear grain rheology (n'=3) relies on <c^k> for k=2,4,6,8
        call f_ev_ck(nlm, 'f', ev_c2,ev_c4,ev_c6,ev_c8) ! Calculate structure tensors of orders 2,4,6,8
    end if

    Evw = (1-alpha)*Evw_Sac(vw, tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca, nprime) &
            + alpha*Evw_Tay(vw, tau, ev_c2,ev_c4,             Ecc,Eca, nprime)
end

function Eeiej(nlm, e1,e2,e3, Ecc,Eca,alpha,nprime)

    ! Enhancement factors in directions (ei,ej) 
    ! (3x3 symmetric matrix of enhancement factors)

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), dimension(3)  :: e1,e2,e3
    real(kind=dp), intent(in)    :: Ecc, Eca, alpha
    integer, intent(in)          :: nprime
    real(kind=dp)                :: Eeiej(3,3)
    
    ! Longitudinal
    Eeiej(1,1) = Evw(outerprod(e1,e1), tau_vv(e1),     nlm, Ecc,Eca,alpha,nprime) 
    Eeiej(2,2) = Evw(outerprod(e2,e2), tau_vv(e2),     nlm, Ecc,Eca,alpha,nprime)
    Eeiej(3,3) = Evw(outerprod(e3,e3), tau_vv(e3),     nlm, Ecc,Eca,alpha,nprime)    
    
    ! Shear
    Eeiej(1,2) = Evw(outerprod(e1,e2), tau_vw(e1,e2),  nlm, Ecc,Eca,alpha,nprime) 
    Eeiej(1,3) = Evw(outerprod(e1,e3), tau_vw(e1,e3),  nlm, Ecc,Eca,alpha,nprime) 
    Eeiej(2,3) = Evw(outerprod(e2,e3), tau_vw(e2,e3),  nlm, Ecc,Eca,alpha,nprime)
    
    ! Symmetric matrix
    Eeiej(2,1) = Eeiej(1,2)
    Eeiej(3,1) = Eeiej(1,3) 
    Eeiej(3,2) = Eeiej(2,3)   
end

function Evw_Sac(vw, tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca, nprime)

    implicit none
    
    real(kind=dp), intent(in) :: Ecc, Eca, vw(3,3), tau(3,3)
    integer, intent(in)       :: nprime
    real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    real(kind=dp)             :: Evw_Sac
    
    Evw_Sac = doubleinner22(ev_epsprime_Sac(tau, ev_c2,    ev_c4,    ev_c6,    ev_c8,     Ecc,Eca,nprime), vw) / &
              doubleinner22(ev_epsprime_Sac(tau, ev_c2_iso,ev_c4_iso,ev_c6_iso,ev_c8_iso, Ecc,Eca,nprime), vw)
end

function Evw_Tay(vw, tau, ev_c2,ev_c4, Ecc,Eca, nprime)

    implicit none
    
    real(kind=dp), intent(in) :: Ecc, Eca, vw(3,3), tau(3,3)
    integer, intent(in)       :: nprime
    real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3)
    real(kind=dp)             :: Evw_Tay
    
    Evw_Tay = doubleinner22(ev_epsprime_Tay(tau, ev_c2,    ev_c4,     Ecc,Eca,nprime), vw) / &
              doubleinner22(ev_epsprime_Tay(tau, ev_c2_iso,ev_c4_iso, Ecc,Eca,nprime), vw)
end

!---------------------------------
! GRAIN-AVERAGED RHEOLOGY (SACHS, TAYLOR)
!---------------------------------

function eps_of_tau__nlin_Sachs(tau, nlm, Aprime,Ecc,Eca) result(eps)

    ! Sachs grain-averaged rheology:
    !       eps = <eps'(tau)>
    ! ...where eps'(tau') is the transversely isotropic rheology for n'=3

    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in)    :: tau(3,3), Aprime, Ecc, Eca
    real(kind=dp)                :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    real(kind=dp)                :: eps(3,3)
    
    call f_ev_ck(nlm, 'f', ev_c2,ev_c4,ev_c6,ev_c8) ! Calculate structure tensors of orders 2,4,6,8
    eps = Aprime*ev_epsprime_Sac(tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca,3)
end

function eps_of_tau__lin_TaylorSachs(tau, nlm, Aprime,Ecc,Eca,alpha) result(eps)

    ! Mixed Taylor--Sachs grain-averaged rheology:
    !       eps = (1-alpha)*<eps'(tau)> + alpha*eps(<tau'>)
    ! ...where eps'(tau') is the transversely isotropic rheology for n'=1

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in)    :: tau(3,3), Aprime, Ecc, Eca, alpha
    real(kind=dp)                :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    real(kind=dp)                :: eps(3,3)
    
    ! Linear grain rheology (n'=1) relies only on <c^2> and <c^4>
    call f_ev_ck(nlm, 'r', ev_c2,ev_c4,ev_c6,ev_c8) ! Calculate structure tensors of orders 2,4 (6,8 are assumed zero for faster evaluation)

    eps = (1-alpha)*Aprime*ev_epsprime_Sac(tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca,1)  &
            + alpha*Aprime*ev_epsprime_Tay(tau, ev_c2,ev_c4,             Ecc,Eca,1) 
end

function ev_epsprime_Sac(tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca,nprime) 
    
    implicit none
    
    real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    real(kind=dp), intent(in) :: Ecc, Eca, tau(3,3)
    integer, intent(in)       :: nprime
    real(kind=dp), parameter  :: d = 3.0
    real(kind=dp)             :: ev_etac0 = 0, ev_etac2(3,3), ev_etac4(3,3,3,3)
    real(kind=dp)             :: ev_epsprime_Sac(3,3), coefA,coefB,coefC, tausq(3,3), I2

    call tranisotropic_coefs(Ecc,Eca,d,nprime,1.0d0, coefA,coefB,coefC)

    I2 = doubleinner22(tau,tau)
    tausq = matmul(tau,tau)

    ! Linear grain fludity    
    if (nprime .eq. 1) then
        ev_etac0 = 1.0
        ev_etac2 = ev_c2
        ev_etac4 = ev_c4

    ! Nonlinear orientation-dependent grain fludity (Johnson's rheology)        
    else if (nprime .eq. 3) then
        ev_etac0 = I2*  1.0 + coefB*doubleinner22(doubleinner42(ev_c4,tau),tau) + 2*coefC*doubleinner22(ev_c2,tausq)
        ev_etac2 = I2*ev_c2 + coefB*doubleinner42(doubleinner62(ev_c6,tau),tau) + 2*coefC*doubleinner42(ev_c4,tausq)
        ev_etac4 = I2*ev_c4 + coefB*doubleinner62(doubleinner82(ev_c8,tau),tau) + 2*coefC*doubleinner62(ev_c6,tausq)

    ! Same as n'=3 but *without* the orientation-dependent terms in the nonlinear grain fluidity.        
    else if (nprime .eq. -3) then 
        ev_etac0 = I2*  1.0
        ev_etac2 = I2*ev_c2
        ev_etac4 = I2*ev_c4
        
    end if

    ev_epsprime_Sac = ev_etac0*tau - coefA*doubleinner22(ev_etac2,tau)*identity + coefB*doubleinner42(ev_etac4,tau) + coefC*(matmul(tau,ev_etac2)+matmul(ev_etac2,tau))
end

function ev_epsprime_Tay(tau, ev_c2,ev_c4, Ecc,Eca,nprime) 
    
    ! Taylor model supports only n'=1. nprime is required any way for future compatibility with n'>1.
    
    implicit none
    
    real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3)
    real(kind=dp), intent(in) :: Ecc, Eca, tau(3,3)
    integer, intent(in)       :: nprime ! Dummy variable: <eps'(tau)> with Taylor hypothesis is implemented only for n' = 1.
    real(kind=dp), parameter  :: d = 3.0
    real(kind=dp)             :: P(9,9), tau_vec(9,1),  P_reg(9,9),tau_vec_reg(9,1)
    real(kind=dp)             :: ev_epsprime_Tay(3,3), coefA,coefB,coefC
    integer                   :: info
!    integer :: ipiv(9), work

    call tranisotropic_coefs(Ecc,Eca,d,nprime,-1.0d0, coefA,coefB,coefC)

    include "include/Taylor_n1.f90"

    tau_vec = reshape(tau, [9,1])
    call dposv('L', 9, 1, P, 9, tau_vec, 9, info) ! tau_vec is now "eps_vec" solution. For some reason, "U" does not work..
!    call dsysv('L', 9, 1, P, 9, ipiv, tau_vec, 9, work,8,info) ! Can be deleted

    if (info /= 0) then
        P_reg       = matmul(TRANSPOSE(P),P) + 1e-6*identity9
        tau_vec_reg = matmul(TRANSPOSE(P),tau_vec)
        call dposv('L', 9, 1, P_reg, 9, tau_vec_reg, 9, info)
        tau_vec = tau_vec_reg
        if (info /= 0) then
            stop 'specfab error: Taylor viscosity-matrix inversion failed! Please check the ODF is correct (reducing the fabric integration time-step, and/or increasing regularization, for transient problems can often help).'
        end if
    end if
    ev_epsprime_Tay = reshape(tau_vec, [3,3])
end

!---------------------------------
! SYNTHETIC STRESS STATES
!---------------------------------

function tau_vv(v) 
    ! v--v compression/extension
    implicit none
    real(kind=dp), intent(in) :: v(3)
    real(kind=dp) :: tau_vv(3,3)
    tau_vv = identity/3 - outerprod(v,v)
end

function tau_vw(v,w)
    ! v--w shear
    implicit none
    real(kind=dp), intent(in) :: v(3), w(3)
    real(kind=dp) :: tau_vw(3,3)
    tau_vw = outerprod(v,w) + outerprod(w,v)
end

!---------------------------------
! QUADRICS
!---------------------------------

function quad_rr(M) result (q2m)
    
    ! Expansion coefs of "M : rr" assuming M is symmetric and r is the radial unit vector
    
    implicit none
    
    real(kind=dp), intent(in) :: M(3,3) 
    real(kind=dp), parameter :: fsq = sqrt(2*Pi/15) 
    complex(kind=dp) :: q2m(-2:2)
    
    ! components (2,-2), (2,-1), (2,0), (2,+1), (2,+2)
    q2m = [    +fsq*   (M(x,x)-M(y,y)+2*i*M(x,y)), &
               +2*fsq* (M(x,z)+i*M(y,z)), &
               -2./3*sqrt(Pi/5)*r*(M(x,x)+M(y,y)-2*M(z,z)), &
               -2*fsq* (M(x,z)-i*M(y,z)), &
               +fsq*   (M(x,x)-M(y,y)-2*i*M(x,y))]
end
   
function quad_tp(M) result (q1m)
    
    ! Expansion coefs of "M : tp" assuming M is anti-symmetric and (t,p) are the (theta,phi) unit vectors.
    
    implicit none
    
    real(kind=dp), intent(in) :: M(3,3) 
    real(kind=dp), parameter :: fsq1 = sqrt(2*Pi/3)
    complex(kind=dp) :: q1m(-1:1)
    
    ! components (1,-1), (1,0), (1,+1)
    q1m = fsq1*[r*M(y,z)-i*M(x,z), r*sqrt(2.)*M(x,y), -r*M(y,z)-i*M(x,z)]
end

!---------------------------------
! REDUCED AND FULL FORM
!---------------------------------

! The ODF being real implies that n_l^{-m} = (-1)^m  conj(n_l^m)
! Hence, coefficients for negative m can be concontructed from the coefficients for m.

function nlm_reduced(nlm_full) 

    implicit none

    complex(kind=dp), intent(in) :: nlm_full(nlm_len)
    complex(kind=dp)             :: nlm_reduced(nlm_reduced_len)
    
    nlm_reduced = nlm_full(I_reduced(1:nlm_reduced_len))
end

function nlm_full(nlm_reduced) 

    implicit none

    complex(kind=dp), intent(in) :: nlm_reduced(nlm_reduced_len)
    complex(kind=dp)             :: nlm_full(nlm_len)

    nlm_full(1:nlm_len) = Mdiag(1:nlm_len) * nlm_reduced(I_full(1:nlm_len)) 
    nlm_full(I_reduced_neg(1:nlm_reduced_neg_len)) = conjg( nlm_full(I_reduced_neg(1:nlm_reduced_neg_len)) )

end

!---------------------------------
! EXTERNAL, NON-CORE FEATURES
!---------------------------------

! Elmer ice flow model 
include "elmer/specfab_elmer.f90"

! JOSEF ice flow model 
include "josef/specfab_josef.f90"

end module specfab 
