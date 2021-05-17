! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

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

    ! Expansion series --- n(theta,phi) = sum_{l,m}^{Lcap,:} n_l^m Y_l^m(theta,phi) where nlm (= n_l^m) = (n_0^0, n_2^-2, n_2^-1, n_2^0, n_2^1, n_2^2, n_4^-4, ... ) 
    integer, private :: Lcap    ! Truncation "L" of expansion series (internal copy of what was passed to the init routine).
    integer          :: nlm_len ! Total number of expansion coefficients (i.e. DOFs)
    integer, parameter :: Lcap__max          = 60
    integer, parameter :: nlm_len__max       = sum([(1+ll*2, ll=0, Lcap__max,2)])
    integer, parameter :: lm(2,nlm_len__max) = reshape([( (ll,mm, mm=-ll,ll), ll=0,  Lcap__max,2)], [2,nlm_len__max]) ! These are the (l,m) values corresponding to the coefficients in "nlm".
    integer, parameter :: I_l0=1, I_l2=I_l0+1, I_l4=I_l2+(2*2+1), I_l6=I_l4+(2*4+1), I_l8=I_l6+(2*6+1), I_l10=I_l8+(2*8+1) ! Indices for extracting l=0,2,4,6,8 coefs of nlm

    ! Dynamics
    complex(kind=dp), private, allocatable :: regmat(:,:) ! Unscaled regularization matrix

    ! Ehancement-factor related
    real(kind=dp), private :: ev_c2_iso(3,3), ev_c4_iso(3,3,3,3), ev_c6_iso(3,3,3,3, 3,3), ev_c8_iso(3,3,3,3, 3,3,3,3) ! <c^k> for isotropic n(theta,phi)
    real(kind=dp), parameter, private :: identity(3,3)  = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    real(kind=dp), parameter, private :: identity9(9,9) = reshape([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [9,9])
    
    ! Optimal n'=1 (lin) grain parameters 
    ! These are the linear mixed Taylor--Sachs best-fit parameters from Rathmann and Lilien (2021)
    real(kind=dp), parameter :: Eca_opt_lin   = 1d3
    real(kind=dp), parameter :: Ecc_opt_lin   = 1d0
    real(kind=dp), parameter :: alpha_opt_lin = 0.0125
    
    ! Optimal n'=3 (nlin) grain parameters 
    ! These are the nonlinear Sachs-only best-fit parameters (Rathmann et. al, 2021) 
    real(kind=dp), parameter :: Eca_opt_nlin   = 1d2
    real(kind=dp), parameter :: Ecc_opt_nlin   = 1d0
    real(kind=dp), parameter :: alpha_opt_nlin = 0
    
contains      

!---------------------------------
! INIT
!---------------------------------
       
subroutine initspecfab(Lcap_)

    ! Needs to be called once before using the module routines.

    implicit none    
    integer, intent(in) :: Lcap_ ! Trauncation "Lcap" to use in model. Note that regularization is calibrated for 8 <= Lcap <= 40.
    complex(kind=dp) :: nlm_iso(nlm_len__max) = 0.0
    
    Lcap     = Lcap_ ! Save internal copy
    nlm_len  = sum([(1+ll*2, ll=0, Lcap,2)]) ! Number of DOFs (coefs)

    ! Set gaunt coefficients
    call set_gaunts()
    
    ! Set isotropic structure tensors
    nlm_iso(1) = 1
    call ck_moments(nlm_iso, ev_c2_iso,ev_c4_iso,ev_c6_iso,ev_c8_iso)
    
    ! Laplacian regularization (unscaled)
    allocate(regmat(nlm_len,nlm_len))
    regmat = 0.0 ! initialize
    do ii = 1, nlm_len    
        do jj = 1, nlm_len    
            if (ii .eq. jj) then
                regmat(ii,ii) = -(lm(1,ii)*(lm(1,ii)+1))**(1.5) ! if expo > 1 => hyper diffusion
            end if
        end do
    end do
    
end

!---------------------------------
! FABRIC DYNAMICS
!---------------------------------

function dndt_ij_ROT(eps,omg, tau,Aprime,Ecc,Eca, beta)  

    ! *** LATTICE ROTATION ***

    ! Returns matrix "dndt_ij" such that
    ! 
    !       d/dt (nlm)_i = dndt_ij (nlm)_j
    ! 
    ! where "nlm" is the vector of spectral coefficient n_l^m

    implicit none

    real(kind=dp), intent(in) :: eps(3,3), omg(3,3), tau(3,3) ! strain-rate (eps), spin (omg), dev. stress (tau)
    real(kind=dp), intent(in) :: Aprime, Ecc, Eca, beta
    complex(kind=dp) :: dndt_ij_ROT(nlm_len,nlm_len), qe(-2:2), qt(-2:2), qo(-1:1)
    integer, parameter :: lmdyn_len = 6 ! Scope of harmonic interactions is local in wave space (this is *not* a free parameter)
    complex(kind=dp), dimension(lmdyn_len) :: g0,     gz,     gn,     gp
    complex(kind=dp), dimension(lmdyn_len) :: g0_rot, gz_rot, gn_rot, gp_rot
    complex(kind=dp), dimension(lmdyn_len) :: g0_Sac, gz_Sac, gn_Sac, gp_Sac
    complex(kind=dp), dimension(lmdyn_len) :: g0_Tay, gz_Tay, gn_Tay, gp_Tay
    real(kind=dp) :: etaprime

    etaprime = Aprime*doubleinner22(tau,tau) ! Assumes eta' = A'*(I2(tau)) (i.e. n'=3).

    ! Quadric expansion coefficients
    qe = rrquad(eps)
    qt = rrquad(tau)
    qo = tpquad(omg)

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
            dndt_ij_ROT(jj,ii) = -1*sum([( &
                    GC(   kk,ii,jj)*g0(kk) + &
                    GCm(  kk,ii,jj)*gz(kk) + &
                    GC_m1(kk,ii,jj)*gn(kk) + &
                    GC_p1(kk,ii,jj)*gp(kk), kk=1,lmdyn_len)])
        end do
    end do
    
end

function dndt_ij_DRX(nlm, tau)

    ! *** Dynamic recrystalization (DRX) ***
    ! 
    ! Nucleation and migration recrystalization modelled as a decay process (Placidi et al., 2010)

    ! Returns matrix "dndt_ij" such that
    ! 
    !       d/dt (nlm)_i = dndt_ij (nlm)_j
    ! 
    ! where "nlm" is the vector of spectral coefficient n_l^m    

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in) :: tau(3,3) ! Deviatoric stress tensor
    complex(kind=dp) :: dndt_ij_DRX(nlm_len,nlm_len)
    real(kind=dp) :: Davg, Davgtensor(nlm_len,nlm_len)

    dndt_ij_DRX = dndt_ij_DRX_src(tau)
    
    ! <D> (diagonal matrix)
    Davg = doubleinner22(matmul(tau,tau), a2_ij(nlm)) - doubleinner22(tau,doubleinner42(a4_ijkl(nlm),tau)) ! (tau.tau):a2 - tau:a4:tau 
    Davgtensor = 0.0 ! initialize 
    do ii = 1, nlm_len    
        do jj = 1, nlm_len    
            if (ii .eq. jj) then
                Davgtensor(ii,ii) = Davg
            end if
        end do
    end do
    
    ! Add add (nonlinear) sink term such that Gamma/Gamma0 = (D - <D>)/(tau:tau) 
    dndt_ij_DRX = dndt_ij_DRX - Davgtensor/doubleinner22(tau,tau) ! This is Gamma/Gamma_0. The caller must multiply by an appropriate DRX rate factor, Gamma_0(T,tau,eps,...).
end

function dndt_ij_DRX_src(tau)  

    ! Calculates the ODF-independent (linear) source term of DRX 

    implicit none

    real(kind=dp), intent(in) :: tau(3,3) ! DRX rate contant, dev. stress tensor
    complex(kind=dp) :: dndt_ij_DRX_src(nlm_len,nlm_len), qt(-2:2)
    integer, parameter :: lmDRX_len = 1+5+9 ! Scope of harmonic interactions is local in wave space (this is NOT a free parameter)
    real(kind=dp) :: k
    complex(kind=dp), dimension(lmDRX_len) :: g

    ! Quadric expansion coefficients
    qt = rrquad(tau)

    ! Harmonic interaction weights 
    include "include/DRX__body.f90"

    ! D
    do ii = 1, nlm_len    
        do jj = 1, nlm_len 
            dndt_ij_DRX_src(jj,ii) = sum([ (GC(kk,ii,jj) * k*g(kk), kk=1,lmDRX_len) ])
        end do
    end do
    
    dndt_ij_DRX_src = dndt_ij_DRX_src/doubleinner22(tau,tau) 
end

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

    dndt_ij_REG = regmat ! The caller must multiply by the regularization strength f_nu(nu0, eps, ...)
end

function f_nu(nu0, beta, eps, tau, Aprime) result (nu)
    
    ! Regularization strength
    
    implicit none
    
    real(kind=dp), intent(in) :: nu0, beta, eps(3,3),tau(3,3), Aprime
    real(kind=dp) :: nu, etaprime

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
    real(kind=dp) :: nu0, nu
    real(kind=dp), parameter :: L0=10
    
    ! nu0=5.0e-2 ! Used for DEMO scripts
    ! nu0=0.5e-3 ! Medium loose
    ! nu0=0.1e-3 ! Used for JOSEF

    nu = nu0 * (Lcap/L0)**(-1.1) * norm2(reshape(eps,[size(eps)]))
end

function f_nu_tau(nu0, tau) result (nu)
    
    ! Calculates optimal nu as a function of Lcap and norm2(tau).
    
    implicit none
    
    real(kind=dp), intent(in) :: tau(3,3)
    real(kind=dp) :: nu0, nu
    real(kind=dp), parameter :: L0=10
    
    ! nu0=5e1 ! Fitted
    
    nu = nu0 * (Lcap/L0)**(-0.4) * norm2(reshape(tau,[size(tau)]))
end


!---------------------------------
! QUADRICS
!---------------------------------

function rrquad(M) result (q2m)
    
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
   
function tpquad(M) result (q1m)
    
    ! Expansion coefs of "M : tp" assuming M is anti-symmetric and (t,p) are the (theta,phi) unit vectors.
    
    implicit none
    
    real(kind=dp), intent(in) :: M(3,3) 
    real(kind=dp), parameter :: fsq1 = sqrt(2*Pi/3)
    complex(kind=dp) :: q1m(-1:1)
    
    ! components (1,-1), (1,0), (1,+1)
    q1m = fsq1*[r*M(y,z)-i*M(x,z), r*sqrt(2.)*M(x,y), -r*M(y,z)-i*M(x,z)]
end

!---------------------------------
! FABRIC FRAME
!---------------------------------

subroutine eigenframe(nlm, e1,e2,e3, eigvals)

    ! Determine eigen vectors e1,e2,e3 and corresponding eigen values (eigvals(1),eigvals(2),eigvals(3)). 
    ! Ordered by largest eigen value first.
    
    implicit none
    
    integer, parameter :: n=3, l=3*3-1
    complex(kind=dp), intent(in)  :: nlm(nlm_len)
    real(kind=dp), dimension(n), intent(out) :: e1,e2,e3, eigvals
    integer inf
    real(kind=dp) :: e_ij(n,n), work(l)
    
    e_ij = a2_ij(nlm)
    call dsyev('V','U',n,e_ij,n,eigvals,work,l,inf)
    e1 = e_ij(:,3)
    e2 = e_ij(:,2)
    e3 = e_ij(:,1)
    eigvals = [eigvals(3),eigvals(2),eigvals(1)] ! Largest first
end 

subroutine pqframe(nlm, p23,p12,p13, q23,q12,q13)

    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), dimension(3), intent(out) :: p23,p12,p13, q23,q12,q13
    real(kind=dp), dimension(3):: e1,e2,e3, eigvals
    
    call eigenframe(nlm, e1,e2,e3, eigvals)
    
    p23 = (e2+e3)/sqrt(2.0) ! p_{e2,e3}
    q23 = (e2-e3)/sqrt(2.0) ! q_{e2,e3}
    
    p12 = (e1+e2)/sqrt(2.0) ! p_{e1,e2}
    q12 = (e1-e2)/sqrt(2.0) ! q_{e1,e2}
    
    p13 = (e1+e3)/sqrt(2.0) ! p_{e1,e3}
    q13 = (e1-e3)/sqrt(2.0) ! q_{e1,e3} 
end

! Collection of routines used for numerical ice-flow model in FEniCS
include "frame.f90"

!---------------------------------
! ORIENTATION (STRUCTURE) TENSORS
!---------------------------------
       
subroutine ck_moments(nlm, ev_c2,ev_c4,ev_c6,ev_c8)
    
    ! "ev_ck" are the tensors <c^k> for a given n(theta,phi) in terms of the expansion coefficients "nlm"
    
    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(inout) :: ev_c2(3,3),ev_c4(3,3,3,3),ev_c6(3,3,3,3, 3,3),ev_c8(3,3,3,3, 3,3,3,3)
    complex(kind=dp) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6), n8m(-8:8)
    
    n00 = nlm(1)
    n2m = nlm(I_l2:(I_l4-1))
    n4m = nlm(I_l4:(I_l6-1))
    n6m = nlm(I_l6:(I_l8-1))
    n8m = nlm(I_l8:(I_l10-1))
    
    ev_c2 = f_ev_c2(n00,n2m)
    ev_c4 = f_ev_c4(n00,n2m,n4m)
    ev_c6 = f_ev_c6(n00,n2m,n4m,n6m)
    ev_c8 = f_ev_c8(n00,n2m,n4m,n6m,n8m)
end
      
! Just some shortcuts for computational efficiency...
      
function a2_ij(nlm) 
    
    ! a^(2) := <c^2> 
    
    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp) :: a2_ij(3,3)
    complex(kind=dp) :: n2m(-2:2) 

    n2m = nlm(I_l2:(I_l4-1))
    a2_ij = f_ev_c2(nlm(1),n2m)
end

function a4_ijkl(nlm) 
    
    ! a^(4) := <c^4> 
    
    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp) :: a4_ijkl(3,3,3,3)
    complex(kind=dp) :: n2m(-2:2), n4m(-4:4)
    
    n2m = nlm(I_l2:(I_l4-1))
    n4m = nlm(I_l4:(I_l6-1))
    a4_ijkl = f_ev_c4(nlm(1),n2m,n4m)
end

!---------------------------------
! SPECTRAL TO TENSORIAL CONVERSION
!---------------------------------

function a2_to_nlm(a2) result(nlm)
    
    ! Given a^(2), returns the equivelent spectral coefs assuming the true ODF is truncated at L=2 (higher order modes vanish)

    implicit none
    
    real(kind=dp), intent(in) :: a2(3,3)
    complex(kind=dp) :: nlm(nlm_len)
    
    nlm = 0.0 ! init
    include "include/a2_to_nlm__body.f90"
end

function a4_to_nlm(a2, a4) result(nlm)

    ! Given a^(2) and a^(4), returns the equivelent spectral coefs assuming the true ODF is truncated at L=4 (higher order modes vanish)

    implicit none
    
    real(kind=dp), intent(in) :: a2(3,3), a4(3,3,3,3)
    complex(kind=dp) :: nlm(nlm_len)
    
    nlm = 0.0 ! init
    include "include/a4_to_nlm__body.f90"
end

function da2dt_DRX(tau, a2, a4)
    
    ! Returns DRX contribution to d/dt a^(2) --- useful routine for external tensorially-based fabric models 
    ! 
    !  *** The caller must multiply with the appropriate rate factor Gamma_0(T, tau, ...) ***
    
    implicit none

    real(kind=dp), intent(in) :: tau(3,3), a2(3,3), a4(3,3,3,3)
    real(kind=dp) :: da2dt_DRX(3,3)
    complex(kind=dp) :: nlm(nlm_len), ddt_nlm(nlm_len)
    
    ! (1) tensorial --> spectral
    nlm = a4_to_nlm(a2, a4)
    
    ! (2) Calculate spectral evolution due to DRX
    ddt_nlm = matmul(dndt_ij_DRX(nlm, tau), nlm) ! d/dt(nlm) = M_ij nlm_j 
    
    ! (3) spectral --> tensorial 
    da2dt_DRX = f_ev_c2( (1d0,0)/Sqrt(4*Pi), nlm(I_l2:(I_l4-1)) ) ! ...Assumes normalized ODF: n00 = 1/sqrt(4*pi) 
    da2dt_DRX = da2dt_DRX - identity/3.0 ! Remove time-constant isotropic (monopole) part due to calculating <c^2> using f_ev_c2() 
end

!---------------------------------
! ENHANCEMENT-FACTORS
!---------------------------------

function Eeiej(nlm, Ecc,Eca,alpha,nprime)

    ! Enhancement-factors in eigen directions: E_{e_i e_j}

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in) :: Ecc, Eca, alpha
    integer, intent(in) :: nprime    
    real(kind=dp), dimension(3) :: e1,e2,e3, eigvals
    real(kind=dp) :: Eeiej(3,3)
    
    call eigenframe(nlm, e1,e2,e3, eigvals)
    
    ! Stored column-wise
    
    Eeiej(1,1) = Evw(outerprod(e1,e1), tau_vv(e1),     nlm, Ecc,Eca,alpha,nprime) ! Longitidinal
    Eeiej(2,1) = Evw(outerprod(e1,e2), tau_vw(e1,e2),  nlm, Ecc,Eca,alpha,nprime) ! Shear
    Eeiej(3,1) = Evw(outerprod(e1,e3), tau_vw(e1,e3),  nlm, Ecc,Eca,alpha,nprime) ! Shear

    Eeiej(1,2) = Evw(outerprod(e2,e1), tau_vw(e2,e1),  nlm, Ecc,Eca,alpha,nprime)
    Eeiej(2,2) = Evw(outerprod(e2,e2), tau_vv(e2),     nlm, Ecc,Eca,alpha,nprime)
    Eeiej(3,2) = Evw(outerprod(e2,e3), tau_vw(e2,e3),  nlm, Ecc,Eca,alpha,nprime)

    Eeiej(1,3) = Evw(outerprod(e3,e1), tau_vw(e3,e1),  nlm, Ecc,Eca,alpha,nprime)
    Eeiej(2,3) = Evw(outerprod(e3,e2), tau_vw(e3,e2),  nlm, Ecc,Eca,alpha,nprime)
    Eeiej(3,3) = Evw(outerprod(e3,e3), tau_vv(e3),     nlm, Ecc,Eca,alpha,nprime)
end

function Exixj(nlm, Ecc,Eca,alpha,nprime) 

    ! Enhancement-factors cartesian directions: E_{x_i x_j}

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in) :: Ecc, Eca, alpha
    integer, intent(in) :: nprime    
    real(kind=dp), dimension(3), parameter :: e1=[1,0,0],e2=[0,1,0],e3=[0,0,1]
    real(kind=dp) :: Exixj(3,3)

    ! Stored column-wise
    
    Exixj(1,1) = Evw(outerprod(e1,e1), tau_vv(e1),     nlm, Ecc,Eca,alpha,nprime) ! Longitidinal
    Exixj(2,1) = Evw(outerprod(e1,e2), tau_vw(e1,e2),  nlm, Ecc,Eca,alpha,nprime) ! Shear
    Exixj(3,1) = Evw(outerprod(e1,e3), tau_vw(e1,e3),  nlm, Ecc,Eca,alpha,nprime) ! Shear

    Exixj(1,2) = Evw(outerprod(e2,e1), tau_vw(e2,e1),  nlm, Ecc,Eca,alpha,nprime)
    Exixj(2,2) = Evw(outerprod(e2,e2), tau_vv(e2),     nlm, Ecc,Eca,alpha,nprime)
    Exixj(3,2) = Evw(outerprod(e2,e3), tau_vw(e2,e3),  nlm, Ecc,Eca,alpha,nprime)

    Exixj(1,3) = Evw(outerprod(e3,e1), tau_vw(e3,e1),  nlm, Ecc,Eca,alpha,nprime)
    Exixj(2,3) = Evw(outerprod(e3,e2), tau_vw(e3,e2),  nlm, Ecc,Eca,alpha,nprime)
    Exixj(3,3) = Evw(outerprod(e3,e3), tau_vv(e3),     nlm, Ecc,Eca,alpha,nprime)
end

function Epijqij(nlm, Ecc,Eca,alpha,nprime)

    ! Enhancement-factors in rotated eigen-directions: E_{p_ij q_ij}

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in) :: Ecc, Eca, alpha
    integer, intent(in) :: nprime    
    real(kind=dp), dimension(3) :: e1,e2,e3, eigvals
    real(kind=dp), dimension(3) :: p23,p12,p13, q23,q12,q13
    real(kind=dp), dimension(3) :: Epijqij

    call eigenframe(nlm, e1,e2,e3, eigvals)    
    call pqframe(nlm, p23,p12,p13, q23,q12,q13)
    
    Epijqij(1) = Evw(outerprod(p23,q23), tau_vw(p23,q23), nlm, Ecc,Eca,alpha,nprime) 
    Epijqij(2) = Evw(outerprod(p12,q12), tau_vw(p12,q12), nlm, Ecc,Eca,alpha,nprime) 
!    Epijqij(2) = Evw(outerprod(e1,e3), tau_vv(p13), nlm, Ecc,Eca,alpha,nprime) ! Override with E^{dagger}_{mt}(tau_0[I/3-pp])
    Epijqij(3) = Evw(outerprod(p13,q13), tau_vw(p13,q13), nlm, Ecc,Eca,alpha,nprime)
end

function Evw(vw, tau, nlm, Ecc, Eca, alpha, nprime)

    ! This is THE generalized directional enhancement factor

    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in) :: Ecc, Eca, alpha, vw(3,3), tau(3,3)
    integer, intent(in) :: nprime
    real(kind=dp) :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    real(kind=dp) :: Evw
    
    call ck_moments(nlm, ev_c2,ev_c4,ev_c6,ev_c8)
    if (alpha < 1e-10) then
        Evw = Evw_Sac(vw, tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca, nprime) ! Faster execution of only Sachs case is of interest
    else
        Evw = (1-alpha)*Evw_Sac(vw, tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca, nprime) + alpha*Evw_Tay(vw, tau, ev_c2,ev_c4, Ecc,Eca, nprime)
    end if
end

function Evw_Sac(vw, tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca, nprime)

    implicit none
    
    real(kind=dp), intent(in) :: Ecc, Eca, vw(3,3), tau(3,3)
    integer, intent(in) :: nprime
    real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    real(kind=dp) :: Evw_Sac
    
    Evw_Sac = doubleinner22(epsprime_mean_Sac(tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca,nprime), vw) / doubleinner22(epsprime_mean_Sac(tau, ev_c2_iso,ev_c4_iso,ev_c6_iso,ev_c8_iso, Ecc,Eca,nprime), vw)
end

function Evw_Tay(vw, tau, ev_c2,ev_c4, Ecc,Eca, nprime)

    implicit none
    
    real(kind=dp), intent(in) :: Ecc, Eca, vw(3,3), tau(3,3)
    integer, intent(in) :: nprime
    real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3)
    real(kind=dp) :: Evw_Tay
    
    Evw_Tay = doubleinner22(epsprime_mean_Tay(tau, ev_c2,ev_c4, Ecc,Eca,nprime), vw) / doubleinner22(epsprime_mean_Tay(tau, ev_c2_iso,ev_c4_iso, Ecc,Eca,nprime), vw)
end

!---------------------------------
! GRAIN-AVERAGED RHEOLOGY (SACHS, TAYLOR)
!---------------------------------

function epsprime_mean_Sac(tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca,nprime) 
    
    implicit none
    
    real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    real(kind=dp), intent(in) :: Ecc, Eca, tau(3,3)
    integer, intent(in) :: nprime
    real(kind=dp), parameter :: d = 3.0
    real(kind=dp) :: ev_etac0 = 0, ev_etac2(3,3), ev_etac4(3,3,3,3)
    real(kind=dp) :: epsprime_mean_Sac(3,3), coefA,coefB,coefC, tausq(3,3), I2

    coefA = (Ecc-1)/(d-1.)
    coefB = (d*(Ecc+1)-2.)/(d-1.) - 2*Eca
    coefC = Eca - 1

    I2 = doubleinner22(tau,tau)
    tausq = matmul(tau,tau)
    
    if (nprime .eq. 1) then
        ev_etac0 = 1.0
        ev_etac2 = ev_c2
        ev_etac4 = ev_c4
    else if (nprime .eq. 3) then
        ev_etac0 = I2*  1.0 + coefB*doubleinner22(doubleinner42(ev_c4,tau),tau) + 2*coefC*doubleinner22(ev_c2,tausq)
        ev_etac2 = I2*ev_c2 + coefB*doubleinner42(doubleinner62(ev_c6,tau),tau) + 2*coefC*doubleinner42(ev_c4,tausq)
        ev_etac4 = I2*ev_c4 + coefB*doubleinner62(doubleinner82(ev_c8,tau),tau) + 2*coefC*doubleinner62(ev_c6,tausq)
    else if (nprime .eq. -3) then ! n'=3 but without the orientation-dependent terms in the grain fluidity.
        ev_etac0 = I2*  1.0
        ev_etac2 = I2*ev_c2
        ev_etac4 = I2*ev_c4
    end if

    epsprime_mean_Sac = ev_etac0*tau - coefA*doubleinner22(ev_etac2,tau)*identity + coefB*doubleinner42(ev_etac4,tau) + coefC*(matmul(tau,ev_etac2)+matmul(ev_etac2,tau))
end

function epsprime_mean_Tay(tau, ev_c2,ev_c4, Ecc,Eca,nprime) 
    
    implicit none
    
    real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3)
    real(kind=dp), intent(in) :: Ecc, Eca, tau(3,3)
    integer, intent(in) :: nprime ! Dummy variable: <eps'(tau)> with Taylor hypothesis is implemented only for n' = 1.
    real(kind=dp), parameter :: d = 3.0
    real(kind=dp) :: P(9,9), tau_vec(9,1),  P_reg(9,9),tau_vec_reg(9,1)
    real(kind=dp) :: epsprime_mean_Tay(3,3), coefA,coefB,coefC
    integer :: info
!    integer :: ipiv(9), work

    coefA = (Ecc**(-1)-1)/(d-1.)
    coefB = (d*(Ecc**(-1)+1)-2.)/(d-1.) - 2*Eca**(-1)
    coefC = Eca**(-1) - 1

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
    epsprime_mean_Tay = reshape(tau_vec, [3,3])
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
! TRANSVERSELY ISOTROPIC RHEOLOGY
!---------------------------------

function eps_of_tau(tau,m,A,Emm,Emt,n) result(eps)

    implicit none
    real(kind=dp), intent(in) :: tau(3,3), m(3), A, Emm, Emt
    integer, intent(in) :: n
    real(kind=dp) :: eps(3,3), fluidity, kI,kM,kL, I2,I4,I5, mm(3,3),tausq(3,3)
    real(kind=dp), parameter :: d = 3.0
    
    call tranisocoefs(Emm,Emt,d,1.0d0, kI,kM,kL)

    mm = outerprod(m,m)
    tausq = matmul(tau,tau)    
    I2 = doubleinner22(tau,tau)
    I4 = doubleinner22(tau,mm)
    I5 = doubleinner22(tausq,mm)
    
    if (n .lt. 0) then
        fluidity = A*(I2)**((n-1)/2)
    else
        fluidity = A*(I2 + kM*I4**2 + 2*kL*I5)**((n-1)/2)
    end if
    eps = fluidity*( tau - kI*I4*identity + kM*I4*mm + kL*(matmul(tau,mm)+matmul(mm,tau)) )
end

function tau_of_eps(eps,m,A,Emm,Emt,n) result(tau)

    implicit none
    real(kind=dp), intent(in) :: eps(3,3), m(3), A, Emm, Emt
    integer, intent(in) :: n
    real(kind=dp) :: tau(3,3), viscosity, kI,kM,kL, I2,I4,I5, mm(3,3),epssq(3,3)
    real(kind=dp), parameter :: d = 3.0
    real(kind=dp) :: expo

    expo = (1.0d0-n)/(2.0d0*n)
    call tranisocoefs(Emm,Emt,d,-1.0d0, kI,kM,kL)

    mm = outerprod(m,m)
    epssq = matmul(eps,eps)    
    I2 = doubleinner22(eps,eps)
    I4 = doubleinner22(eps,mm)
    I5 = doubleinner22(epssq,mm)
    
    if (n .lt. 0) then
        viscosity = A**(-1.0d0/n)*(I2)**(expo)
    else
        viscosity = A**(-1.0d0/n)*(I2 + kM*I4**2 + 2*kL*I5)**(expo)
    end if

    tau = viscosity*( eps - kI*I4*identity + kM*I4*mm + kL*(matmul(eps,mm)+matmul(mm,eps)) )
end

subroutine tranisocoefs(Emm,Emt,d,expo, kI,kM,kL)

    implicit none
    real(kind=dp), intent(in) :: Emm,Emt,expo,d
    real(kind=dp), intent(out) :: kI,kM,kL

    kI = (Emm**expo-1)/(d-1.)
    kM = (d*(Emm**expo+1)-2.)/(d-1.) - 2*Emt**expo
    kL = Emt**expo - 1
end

end module specfab 
