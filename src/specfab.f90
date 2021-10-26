! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2019-2021

module specfab  

    use tensorproducts
    use moments
    use gaunt

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real, parameter, private    :: Pi = 3.1415927
    integer, parameter, private :: x = 1, y = 2, z = 3 ! Matrix indices
    complex(kind=dp), parameter, private :: r = (1,0), i = (0,1) ! real and imag units
    integer, private :: ii, jj, kk, ll, mm ! Loop indicies

    ! Expansion series:
    !     n(theta,phi) = sum_{l,m}^{Lcap,:} n_l^m Y_l^m(theta,phi) 
    ! where "nlm" vector := n_l^m = (n_0^0, n_2^-2, n_2^-1, n_2^0, n_2^1, n_2^2, n_4^-4, ... ) 
    integer, private   :: Lcap    ! Truncation "L" of expansion series (internal copy of what was passed to the init routine).
    integer            :: nlm_len ! Total number of expansion coefficients (i.e. DOFs)
    integer, parameter :: Lcap__max          = 60 ! Hard limit
    integer, parameter :: nlm_len__max       = sum([(1+ll*2, ll=0, Lcap__max,2)])
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
    
    Lcap     = Lcap_ ! Save internal copy
    nlm_len  = sum([(1+ll*2, ll=0, Lcap,2)]) ! Number of DOFs (expansion coefficients)

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
    integer, parameter        :: lmdyn_len = 6 ! Scope of harmonic interactions is local in wave space (this is *not* a free parameter)
    complex(kind=dp), dimension(lmdyn_len) :: g0,     gz,     gn,     gp
    complex(kind=dp), dimension(lmdyn_len) :: g0_rot, gz_rot, gn_rot, gp_rot
    complex(kind=dp), dimension(lmdyn_len) :: g0_Sac, gz_Sac, gn_Sac, gp_Sac
    complex(kind=dp), dimension(lmdyn_len) :: g0_Tay, gz_Tay, gn_Tay, gp_Tay
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
    
    g0_Sac = 3*[ 0*r, Eca*qt(-2),Eca*qt(-1),Eca*qt(0),Eca*qt(+1),Eca*qt(+2) ]
    gz_Sac = [ 0*r, -Eca*qt(-2),-1./2*(Eca-1)*qt(-1),0*r,1./2*(Eca-1)*qt(+1),Eca*qt(+2) ]
    gn_Sac = [ sqrt(5./6)*qt(-1), 0*r, Eca*qt(-2), sqrt(1./6)*(3*Eca-1)*qt(-1), sqrt(3./2)*Eca*qt(0), (Eca+1)*qt(+1) ]
    gp_Sac = [ sqrt(5./6)*qt(+1), (Eca+1)*qt(-1), sqrt(3./2)*Eca*qt(0), sqrt(1./6)*(3*Eca-1)*qt(+1), Eca*qt(+2), 0*r ]

    g0 = g0_rot + (1-beta)*etaprime*g0_Sac + beta*g0_Tay
    gz = gz_rot + (1-beta)*etaprime*gz_Sac + beta*gz_Tay
    gn = gn_rot + (1-beta)*etaprime*gn_Sac + beta*gn_Tay
    gp = gp_rot + (1-beta)*etaprime*gp_Sac + beta*gp_Tay


    ! Contruct rate-of-change matrix
    do ii = 1, nlm_len    
        do jj = 1, nlm_len    
            dndt_ij_LATROT(jj,ii) = -1*sum([( &
                    GC(   kk,ii,jj)*g0(kk) + &
                    GCm(  kk,ii,jj)*gz(kk) + &
                    GC_m1(kk,ii,jj)*gn(kk) + &
                    GC_p1(kk,ii,jj)*gp(kk), kk=1,lmdyn_len)])
        end do
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

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in)    :: tau(3,3) ! Deviatoric stress tensor
    complex(kind=dp)             :: dndt_ij_DDRX(nlm_len,nlm_len)
    real(kind=dp)                :: Davg, Davgtensor(nlm_len,nlm_len)

    dndt_ij_DDRX = dndt_ij_DDRX_src(tau)
    
    ! <D> (diagonal matrix)
    Davg = doubleinner22(matmul(tau,tau), a2(nlm)) - doubleinner22(tau,doubleinner42(a4(nlm),tau)) ! (tau.tau):a2 - tau:a4:tau 
    Davgtensor = 0.0 ! initialize 
    do ii = 1, nlm_len    
        do jj = 1, nlm_len    
            if (ii .eq. jj) then
                Davgtensor(ii,ii) = Davg
            end if
        end do
    end do
    
    ! Add add (nonlinear) sink term such that Gamma/Gamma0 = (D - <D>)/(tau:tau) 
    dndt_ij_DDRX = dndt_ij_DDRX - Davgtensor/doubleinner22(tau,tau) ! This is Gamma/Gamma_0. The caller must multiply by an appropriate DDRX rate factor, Gamma_0(T,tau,eps,...).
end

function dndt_ij_DDRX_src(tau)  

    ! Calculates the ODF-independent (linear) source term of DDRX 

    implicit none

    real(kind=dp), intent(in) :: tau(3,3) ! dev. stress tensor
    complex(kind=dp)          :: dndt_ij_DDRX_src(nlm_len,nlm_len), qt(-2:2)
    integer, parameter        :: lmDDRX_len = 1+5+9 ! Scope of harmonic interactions is local in wave space (fixed parameter)
    real(kind=dp)             :: k
    complex(kind=dp)          :: g(lmDDRX_len)

    ! Quadric expansion coefficients
    qt = quad_rr(tau)

    ! Harmonic interaction weights 
    include "include/DDRX__body.f90"

    ! D
    do ii = 1, nlm_len    
        do jj = 1, nlm_len 
            dndt_ij_DDRX_src(jj,ii) = sum([ (GC(kk,ii,jj) * k*g(kk), kk=1,lmDDRX_len) ])
        end do
    end do
    
    dndt_ij_DDRX_src = dndt_ij_DDRX_src/doubleinner22(tau,tau) 
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
            stop 'Matrix is numerically singular!'
        end if
    end if
    ev_epsprime_Tay = reshape(tau_vec, [3,3])
end

!---------------------------------
! TRANSVERSELY ISOTROPIC RHEOLOGY
!---------------------------------

function eps_of_tau__tranisotropic(tau, A,n, m,Emm,Emt) result(eps)

    implicit none
    real(kind=dp), intent(in) :: tau(3,3), A, m(3), Emm, Emt
    integer, intent(in)       :: n
    real(kind=dp)             :: eps(3,3), fluidity, kI,kM,kL, I2,I4,I5, mm(3,3),tausq(3,3)
    real(kind=dp), parameter  :: d = 3.0
    real(kind=dp)             :: expo
    
    call tranisotropic_coefs(Emm,Emt,d,n,1.0d0, kI,kM,kL)

    mm = outerprod(m,m)
    tausq = matmul(tau,tau)    
    I2 = doubleinner22(tau,tau)
    I4 = doubleinner22(tau,mm)
    I5 = doubleinner22(tausq,mm)
    
    expo = (n-1)/2
    fluidity = A*(I2 + kM*I4**2 + 2*kL*I5)**(expo)
    eps = fluidity*( tau - kI*I4*identity + kM*I4*mm + kL*(matmul(tau,mm)+matmul(mm,tau)) )
end

function tau_of_eps__tranisotropic(eps, A,n, m,Emm,Emt) result(tau)

    implicit none
    real(kind=dp), intent(in) :: eps(3,3), A, m(3), Emm, Emt
    integer, intent(in)       :: n
    real(kind=dp)             :: tau(3,3), viscosity, kI,kM,kL, I2,I4,I5, mm(3,3),epssq(3,3)
    real(kind=dp), parameter  :: d = 3.0
    real(kind=dp)             :: expo

    call tranisotropic_coefs(Emm,Emt,d,n,-1.0d0, kI,kM,kL)

    mm = outerprod(m,m)
    epssq = matmul(eps,eps)    
    I2 = doubleinner22(eps,eps)
    I4 = doubleinner22(eps,mm)
    I5 = doubleinner22(epssq,mm)

    expo = (1.0d0-n)/(2.0d0*n)
    viscosity = A**(-1.0d0/n)*(I2 + kM*I4**2 + 2*kL*I5)**(expo)
    tau = viscosity*( eps - kI*I4*identity + kM*I4*mm + kL*(matmul(eps,mm)+matmul(mm,eps)) )
end

subroutine tranisotropic_coefs(Emm,Emt,d,n,expo, kI,kM,kL)

    implicit none
    real(kind=dp), intent(in)  :: Emm,Emt,expo,d
    integer, intent(in)        :: n
    real(kind=dp), intent(out) :: kI,kM,kL
    real(kind=dp)              :: nexpo

    nexpo = expo * 2.0d0/(n + 1.0d0) 
    
    kI = (Emm**nexpo-1)/(d-1.)
    kM = (d*(Emm**nexpo+1)-2.)/(d-1.) - 2*Emt**nexpo
    kL = Emt**nexpo - 1
end

!---------------------------------
! ORTHOTROPIC RHEOLOGY
!---------------------------------

function eps_of_tau__orthotropic(tau, A,n, m1,m2,m3, Eij) result(eps)

    implicit none
    real(kind=dp), intent(in)     :: tau(3,3), A, m1(3),m2(3),m3(3), Eij(3,3)
    integer, intent(in)           :: n
    real(kind=dp)                 :: eps(3,3)
    real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
    real(kind=dp)                 :: tau11,tau22,tau33,tau23,tau31,tau12
    real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
    real(kind=dp)                 :: I1,I2,I3,I4,I5,I6
    real(kind=dp)                 :: fluidity

    call orthotropic_coefs(n, Eij, lam1,lam2,lam3,lam4,lam5,lam6, gam)
    call orthotropic_tensors_and_invars(tau, m1,m2,m3, M11,M22,M33,M23,M31,M12, tau11,tau22,tau33,tau23,tau31,tau12)

    I1 = (tau22 - tau33)/2.0d0 
    I2 = (tau33 - tau11)/2.0d0 
    I3 = (tau11 - tau22)/2.0d0 
    
    ! Strictly speaking, these are sqrt(I4), sqrt(I5), and sqrt(I6)
    I4 = tau23  
    I5 = tau31
    I6 = tau12

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

function tau_of_eps__orthotropic(eps, A,n, m1,m2,m3, Eij) result(tau)

    implicit none
    real(kind=dp), intent(in)     :: eps(3,3), A, m1(3),m2(3),m3(3), Eij(3,3)
    integer, intent(in)           :: n
    real(kind=dp)                 :: tau(3,3)
    real(kind=dp), dimension(3,3) :: M11,M22,M33,M23,M31,M12
    real(kind=dp)                 :: eps11,eps22,eps33,eps23,eps31,eps12
    real(kind=dp)                 :: lam1,lam2,lam3,lam4,lam5,lam6, gam
    real(kind=dp)                 :: J1,J2,J3,J4,J5,J6, J23,J31,J12
    real(kind=dp)                 :: viscosity

    call orthotropic_coefs(n, Eij, lam1,lam2,lam3,lam4,lam5,lam6, gam)
    call orthotropic_tensors_and_invars(eps, m1,m2,m3, M11,M22,M33,M23,M31,M12, eps11,eps22,eps33,eps23,eps31,eps12)
    
    J1 = (eps22 - eps33)/2.0d0 
    J2 = (eps33 - eps11)/2.0d0 
    J3 = (eps11 - eps22)/2.0d0 
    
    ! Strictly speaking, these are sqrt(I4), sqrt(I5), and sqrt(I6)
    J4 = eps23
    J5 = eps31
    J6 = eps12
    
    J23 = -3/2.0d0 * eps11 ! J23 := J2-J3 = (eps22 + eps33 - 2*eps11)/2.0d0 = -3/2.0d0 * eps11
    J31 = -3/2.0d0 * eps22 ! J31 := J3-J1 = (eps11 + eps33 - 2*eps22)/2.0d0 = -3/2.0d0 * eps22
    J12 = -3/2.0d0 * eps33 ! J12 := J1-J2 = (eps11 + eps22 - 2*eps33)/2.0d0 = -3/2.0d0 * eps33

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

subroutine orthotropic_coefs(n, Eij, lam1,lam2,lam3,lam4,lam5,lam6, gam)

    implicit none
    real(kind=dp), intent(in)  :: Eij(3,3)
    integer, intent(in)        :: n
    real(kind=dp), intent(out) :: lam1,lam2,lam3,lam4,lam5,lam6, gam
    real(kind=dp)              :: Eexpo, E11n,E22n,E33n !,E12n,E31n,E23n

    Eexpo = 2.0d0/(n + 1.0d0) 
    
    E11n = Eij(1,1)**Eexpo
    E22n = Eij(2,2)**Eexpo
    E33n = Eij(3,3)**Eexpo

    lam1 = 4/3.0d0*(-E11n+E22n+E33n)
    lam2 = 4/3.0d0*(+E11n-E22n+E33n)
    lam3 = 4/3.0d0*(+E11n+E22n-E33n)
    lam4 = 2* Eij(2,3)**Eexpo
    lam5 = 2* Eij(3,1)**Eexpo
    lam6 = 2* Eij(1,2)**Eexpo
    
    gam = -(E11n**2 + E22n**2 + E33n**2) + 2*(E11n*E22n + E11n*E33n + E22n*E33n)
end

subroutine orthotropic_tensors_and_invars(tau, m1,m2,m3, M11,M22,M33,M23,M31,M12, tau11,tau22,tau33,tau23,tau31,tau12)

    implicit none    
    real(kind=dp), intent(in)  :: tau(3,3), m1(3),m2(3),m3(3)
    real(kind=dp), intent(out) :: M11(3,3),M22(3,3),M33(3,3),M23(3,3),M31(3,3),M12(3,3)
    real(kind=dp), intent(out) :: tau11,tau22,tau33,tau23,tau31,tau12
   
    M11 = outerprod(m1,m1)
    M22 = outerprod(m2,m2)
    M33 = outerprod(m3,m3)
    M23 = outerprod(m2,m3)
    M31 = outerprod(m3,m1)
    M12 = outerprod(m1,m2)

    tau11 = doubleinner22(tau, M11) 
    tau22 = doubleinner22(tau, M22) 
    tau33 = doubleinner22(tau, M33) 
    tau23 = doubleinner22(tau, M23) 
    tau31 = doubleinner22(tau, M31) 
    tau12 = doubleinner22(tau, M12) 
end

!---------------------------------
! ISOTROPIC (GLEN) RHEOLOGY
!---------------------------------

function eps_of_tau__isotropic(tau, A,n) result(eps)

    implicit none
    real(kind=dp), intent(in)     :: tau(3,3), A
    integer, intent(in)           :: n
    real(kind=dp)                 :: eps(3,3), fluidity

    fluidity = A * (doubleinner22(tau,tau))**((n-1)/2.0d0)
    eps = fluidity * tau
end

function tau_of_eps__isotropic(eps, A,n) result(tau)

    implicit none
    real(kind=dp), intent(in)     :: eps(3,3), A
    integer, intent(in)           :: n
    real(kind=dp)                 :: tau(3,3), viscosity

    viscosity = A**(-1./n) * (doubleinner22(eps,eps))**((1-n)/(2.0d0*n))
    tau = viscosity * eps
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
! EXTERNAL, NON-CORE FEATURES
!---------------------------------

! Elmer ice flow model 
include "elmer/specfab_elmer.f90"

! JOSEF ice flow model 
include "josef/specfab_josef.f90"

end module specfab 
