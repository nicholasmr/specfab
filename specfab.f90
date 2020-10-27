! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020

module specfab  

    use tensorproducts
    use moments
    use gaunt

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real, parameter, private    :: Pi = 3.1415927
    integer, parameter, private :: x = 1, y = 2, z = 3 ! Matrix indices
    complex(kind=dp), parameter, private :: r = (1,0), i = (0,1) ! shortcuts for real and imag units
    integer, private :: ii, jj, kk, ll, mm ! Loop indicies

    ! Expansion series --- n(theta,phi) = sum_{l,m}^{Lcap,:} n_l^m Y_l^m(theta,phi) where nlm (= n_l^m) = (n_0^0, n_2^-2, n_2^-1, n_2^0, n_2^1, n_2^2, n_4^-4, ... ) 
    integer, private :: Lcap    ! Truncation "L" of expansion series (internal copy of what was passed to the init routine).
    integer          :: nlm_len ! Total number of expansion coefficients (i.e. DOFs)
    integer, parameter :: Lcap__max          = 60
    integer, parameter :: nlm_len__max       = sum([(1+ll*2, ll=0, Lcap__max,2)])
    integer, parameter :: lm(2,nlm_len__max) = reshape([( (ll,mm, mm=-ll,ll), ll=0,  Lcap__max,2)], [2,nlm_len__max]) ! These are the (l,m) values corresponding to the coefficients in "nlm".
    integer, parameter :: I_l0=1, I_l2=I_l0+1, I_l4=I_l2+(2*2+1), I_l6=I_l4+(2*4+1), I_l8=I_l6+(2*6+1), I_l10=I_l8+(2*8+1) ! Indices for extracting l=0,2,4,6,8 coefs of nlm

    ! Ehancement-factor related
    real(kind=dp), private :: ev_c2_iso(3,3), ev_c4_iso(3,3,3,3), ev_c6_iso(3,3,3,3, 3,3), ev_c8_iso(3,3,3,3, 3,3,3,3) ! <c^k> for isotropic n(theta,phi)
    real(kind=dp), parameter, private :: identity(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
    
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
end

!---------------------------------
! FABRIC TIME-EVOLUTION
!---------------------------------

function dndt_ij(eps,omg)  

    ! d/dt (nlm)_i = dndt_ij (nlm)_j

    implicit none

    real(kind=dp), intent(in) :: eps(3,3), omg(3,3) ! strain-rate (eps) and spin (omg)
    complex(kind=dp) :: dndt_ij(nlm_len,nlm_len), reg_ij(nlm_len,nlm_len), qe(-2:2), qo(-1:1)
    integer, parameter :: lmdyn_len = 6 ! Scope of harmonic interactions is local in wave space (this is NOT a free parameter)
    complex(kind=dp), dimension(lmdyn_len) :: w, wm, w_m1, w_p1
    real(kind=dp) :: nu

    ! Quadric expansion coefficients
    qe = rrquad(eps)
    qo = tpquad(omg)

    ! Harmonic interaction weights
    w    = 3*[ 0*r, qe(-2),qe(-1),qe(0),qe(+1),qe(+2) ]
    wm   = [ -i*sqrt(3.)*qo(0), -qe(-2),0*r,0*r,0*r,qe(+2) ]
    w_m1 = [ -i*6/sqrt(6.)*qo(-1) + sqrt(5./6)*qe(-1), 0*r, qe(-2), sqrt(2./3)*qe(-1), sqrt(3./2)*qe(0), 2*qe(+1) ]
    w_p1 = [ +i*6/sqrt(6.)*qo(+1) + sqrt(5./6)*qe(+1), 2*qe(-1), sqrt(3./2)*qe(0), sqrt(2./3)*qe(+1), qe(+2), 0*r ]

    ! Contruct rate-of-change matrix
    do ii = 1, nlm_len    
        do jj = 1, nlm_len    
            dndt_ij(jj,ii) = -1*sum([( &
                    GC(   kk,ii,jj)*w(kk) + &
                    GCm(  kk,ii,jj)*wm(kk) + &
                    GC_m1(kk,ii,jj)*w_m1(kk) + &
                    GC_p1(kk,ii,jj)*w_p1(kk), kk=1,lmdyn_len)])
        end do
    end do
    
    ! Laplacian regularization (diagonal matrix)
    reg_ij = 0.0
    nu = f_nu(eps)
    do ii = 1, nlm_len    
        do jj = 1, nlm_len    
            if (ii .eq. jj) then
                reg_ij(ii,ii) = -nu * (lm(1,ii)*(lm(1,ii)+1))**(1.0) ! **(1.5)
            end if
        end do
    end do
    
    dndt_ij = dndt_ij + reg_ij
end

function f_nu(eps) result (nu)
    
    ! Regularization diffusion coefficient.
    ! Calculates optimal nu as a function of Lcap and norm2(eps).
    
    implicit none
    
    real(kind=dp), intent(in) :: eps(3,3)
    real(kind=dp) :: nu
    real(kind=dp), parameter :: L0=10, nu0=5e-2 !  Diffusion exponent = 1.0 --- nu0>=10e-2 seems very safe, 8e-2 is OK.
!    real(kind=dp), parameter :: L0=10, nu0=4e-3 ! Diffusion exponent = 1.5 ---

    if (nu0 .le. 1.0d-20) then ! For debugging the regularization
        nu = 3e-3
    else 
        nu = nu0*(Lcap/L0)**(-1.1)
    end if
    
    nu = nu*norm2(reshape(eps,[size(eps)]))
end

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
! FABRIC STATE
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
       
function a2_ij(nlm) 
    
    ! a^(2) = <c^2> 
    
    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp) :: a2_ij(3,3)
    complex(kind=dp) :: n2m(-2:2) 

    n2m = nlm(I_l2:(I_l4-1))
    a2_ij = f_ev_c2(nlm(1),n2m)
end

!---------------------------------
! ENHANCEMENT-FACTORS
!---------------------------------

function Eeiej(nlm, Ecc,Eca,nprime)

    ! Enhancement-factors in eigen directions: E_{e_i e_j}

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in) :: Ecc, Eca
    integer, intent(in) :: nprime    
    real(kind=dp) :: Eeiej(3,3), e1(3),e2(3),e3(3), eigvals(3)
    
    call eigenframe(nlm, e1,e2,e3, eigvals)
    
    ! Stored column-wise
    
    Eeiej(1,1) = Evw(outerprod(e1,e1), tau_vv(e1),     nlm, Ecc,Eca,nprime) ! Longitidinal
    Eeiej(2,1) = Evw(outerprod(e1,e2), tau_vw(e1,e2),  nlm, Ecc,Eca,nprime) ! Shear
    Eeiej(3,1) = Evw(outerprod(e1,e3), tau_vw(e1,e3),  nlm, Ecc,Eca,nprime) ! Shear

    Eeiej(1,2) = Evw(outerprod(e2,e1), tau_vw(e2,e1),  nlm, Ecc,Eca,nprime)
    Eeiej(2,2) = Evw(outerprod(e2,e2), tau_vv(e2),     nlm, Ecc,Eca,nprime)
    Eeiej(3,2) = Evw(outerprod(e2,e3), tau_vw(e2,e3),  nlm, Ecc,Eca,nprime)

    Eeiej(1,3) = Evw(outerprod(e3,e1), tau_vw(e3,e1),  nlm, Ecc,Eca,nprime)
    Eeiej(2,3) = Evw(outerprod(e3,e2), tau_vw(e3,e2),  nlm, Ecc,Eca,nprime)
    Eeiej(3,3) = Evw(outerprod(e3,e3), tau_vv(e3),     nlm, Ecc,Eca,nprime)

end

function Epijqij(nlm, Ecc,Eca,nprime)

    ! Enhancement-factors in rotated eigen-directions: E_{p_ij q_ij}

    implicit none

    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in) :: Ecc, Eca
    integer, intent(in) :: nprime    
    real(kind=dp), dimension(3) :: Epijqij, p23,p12,p13, q23,q12,q13
    
    call pqframe(nlm, p23,p12,p13, q23,q12,q13)
    
    Epijqij(1) = Evw(outerprod(p23,q23), tau_vw(p23,q23), nlm, Ecc,Eca,nprime) 
    Epijqij(2) = Evw(outerprod(p12,q12), tau_vw(p12,q12), nlm, Ecc,Eca,nprime) 
    Epijqij(3) = Evw(outerprod(p13,q13), tau_vw(p13,q13), nlm, Ecc,Eca,nprime)
end

function Evw(vw, tau, nlm, Ecc, Eca, nprime)

    implicit none
    
    complex(kind=dp), intent(in) :: nlm(nlm_len)
    real(kind=dp), intent(in) :: Ecc, Eca, vw(3,3), tau(3,3)
    integer, intent(in) :: nprime
    real(kind=dp) :: Evw, ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    
    call ck_moments(nlm, ev_c2,ev_c4,ev_c6,ev_c8)
    Evw = doubleinner22(epsprime_mean(tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca,nprime), vw) / doubleinner22(epsprime_mean(tau, ev_c2_iso,ev_c4_iso,ev_c6_iso,ev_c8_iso, Ecc,Eca,nprime), vw)
end

function epsprime_mean(tau, ev_c2,ev_c4,ev_c6,ev_c8, Ecc,Eca,nprime) 
    
    implicit none
    
    real(kind=dp), intent(in) :: ev_c2(3,3), ev_c4(3,3,3,3), ev_c6(3,3,3,3, 3,3), ev_c8(3,3,3,3, 3,3,3,3)
    real(kind=dp), intent(in) :: Ecc, Eca, tau(3,3)
    integer, intent(in) :: nprime
    real(kind=dp), parameter :: d = 3.0
    real(kind=dp) :: ev_etac0 = 0, ev_etac2(3,3), ev_etac4(3,3,3,3)
    real(kind=dp) :: epsprime_mean(3,3), coefA,coefB,coefC, tausq(3,3), I2

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
    end if

    epsprime_mean = ev_etac0*tau - coefA*doubleinner22(ev_etac2,tau)*identity + coefB*doubleinner42(ev_etac4,tau) + coefC*(matmul(tau,ev_etac2)+matmul(ev_etac2,tau))
end

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

end module specfab 
