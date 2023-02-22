! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2022

! Fabric (ODF) dynamics in spectral space.

module dynamics  

    use tensorproducts
    use moments ! used by tensorial dynamics routines
    use gaunt

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real,    parameter, private :: Pi = 3.141592653589793
    integer, parameter, private :: x = 1, y = 2, z = 3 ! Matrix indices (used for readability)
    complex(kind=dp), parameter, private :: r = (1,0), i = (0,1) ! real and imag units
    integer, private :: ii, ll, mm ! Loop indices

    integer, private   :: Lcap    ! Truncation "L" of expansion series (internal copy of what was passed to the init routine).
    integer, private   :: nlm_len ! Total number of expansion coefficients (i.e. DOFs)

    integer, parameter :: Lcap__max  = 30 ! Hard limit
    integer, parameter :: nlm_lenmax = (Lcap__max+1)*(Lcap__max+2)/2
    
    ! Static (constant) matrices used for spectral dynamics
    real(kind=dp), parameter :: Ldiag(nlm_lenmax) = [( (-ll*(ll+1),mm=-ll,ll), ll=0,  Lcap__max,2)] ! Diagonal entries of Laplacian diffusion operator.
    
contains      

    !---------------------------------
    ! INIT
    !---------------------------------

    subroutine initdynamics(Lcap_)

        ! Needs to be called once before using the module routines.

        implicit none    
        integer, intent(in) :: Lcap_ ! Truncation "Lcap"
        
        ! Save internal copy
        Lcap    = Lcap_ 
        nlm_len = (Lcap+1)*(Lcap+2)/2 

        ! Set gaunt coefficients (overlap integrals involving three spherical harmonics)
        call set_gaunts()
    end

    !---------------------------------
    ! FABRIC DYNAMICS
    !---------------------------------

    function M_LROT(eps, omg, iota, zeta)  

        !------------------
        ! Lattice rotation
        !------------------

        ! Returns matrix M such that d/dt (nlm)_i = M_ij (nlm)_j

        implicit none

        real(kind=dp), intent(in) :: eps(3,3), omg(3,3), iota, zeta ! strain-rate (eps), spin (omg), eps^1 coefficient (iota), eps^2 coefficient (zeta)
        complex(kind=dp)          :: M_LROT(nlm_len,nlm_len), qe(-2:2), qo(-1:1)
        integer, parameter        :: SHI_LATROT = 6 ! Scope of harmonic interactions (in wave space) for LATROT
        complex(kind=dp), dimension(SHI_LATROT) :: g0,     gz,     gn,     gp
        complex(kind=dp), dimension(SHI_LATROT) :: g0_rot, gz_rot, gn_rot, gp_rot
        complex(kind=dp), dimension(SHI_LATROT) :: g0_Tay, gz_Tay, gn_Tay, gp_Tay
        real(kind=dp) :: zetanorm, epssq(3,3)
        
        ! Quadric expansion coefficients
        epssq = matmul(eps,eps)
        zetanorm = zeta/sqrt(epssq(1,1)+epssq(2,2)+epssq(3,3)) ! Normalize zeta so that evolution depends only on total accumulated strain
        qe = quad_rr(iota*eps + zetanorm*epssq)
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
        
        g0 = g0_rot + g0_Tay
        gz = gz_rot + gz_Tay
        gn = gn_rot + gn_Tay
        gp = gp_rot + gp_Tay


        do ii = 1, nlm_len
            M_LROT(ii,1:nlm_len) = -1*( matmul(GC(ii,1:nlm_len,1:SHI_LATROT),g0) + matmul(GCm(ii,1:nlm_len,1:SHI_LATROT),gz) + matmul(GC_m1(ii,1:nlm_len,1:SHI_LATROT),gn) + matmul(GC_p1(ii,1:nlm_len,1:SHI_LATROT),gp) )    
        end do
    end

    function M_DDRX(nlm, tau)

        !-----------------------------------------------
        ! Discontinuous dynamic recrystallization (DDRX)
        !----------------------------------------------- 
        
        ! Nucleation and migration recrystalization modeled as a decay process (Placidi et al., 2010).
        ! Returns matrix M such that d/dt (nlm)_i = M_ij (nlm)_j
        ! NOTICE: This is Gamma/Gamma_0. The caller must multiply by an appropriate DDRX rate factor, Gamma_0(T,tau,eps,...).
        
        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        real(kind=dp), intent(in)    :: tau(3,3) ! Deviatoric stress tensor
        complex(kind=dp)             :: M_DDRX(nlm_len,nlm_len)
        real(kind=dp)                :: Davg
        real(kind=dp)                :: tausq(3,3), trtausq ! tau.tau, tau:tau
        real(kind=dp)                :: a4v(6,6), a2v(6), tauv(6), a4tauv(6)

        ! Get linear source term
        M_DDRX = M_DDRX_src(tau)

        ! Calculate nonlinear sink term, <D>
        tausq = matmul(tau,tau) ! = tau.tau
        trtausq = tausq(1,1) + tausq(2,2) + tausq(3,3) ! = tr(tau.tau) = tau:tau
        call f_ev_ck_Mandel(nlm, a2v, a4v) ! Structure tensors in Mandel notation
        tauv = mat_to_vec(tau) ! 6x1 Mandel vector        
        a4tauv = matmul(a4v,tauv) ! = a4:tau in Mandel notation (6x1 vector)
        Davg = dot_product(mat_to_vec(tausq),a2v) - dot_product(tauv,a4tauv) ! = (tau.tau):a2 - tau:a4:tau 
        Davg = Davg/trtausq ! normalize by tau:tau
        
        ! Add source and sink terms
        ! I.e. add <D> to all diagonal entries such that M_ij := Gamma_ij/Gamma0 = (D_ij - <D>*I_ij)/(tau:tau) (where I is the identity)
        do ii = 1, nlm_len    
            M_DDRX(ii,ii) = M_DDRX(ii,ii) - Davg 
        end do
    end

    function M_DDRX_src(tau)  

        !------------------
        ! DDRX source term
        !------------------

        implicit none

        real(kind=dp), intent(in) :: tau(3,3) ! dev. stress tensor
        complex(kind=dp)          :: M_DDRX_src(nlm_len,nlm_len), qt(-2:2)
        integer, parameter        :: SHI_DDRX = 1+5+9 ! Scope of harmonic interactions (in wave space) for DDRX
        real(kind=dp)             :: k
        complex(kind=dp)          :: g(SHI_DDRX)

        ! Quadric expansion coefficients
        qt = quad_rr(tau)

        ! Harmonic interaction weights (requires qt, sets g and k)
        include "include/ddrx-coupling-weights.f90"

        ! D
        g = k*g ! common prefactor
        g = g/doubleinner22(tau,tau) ! normalize (can be done already here since it is a constant prefactor for M_DDRX_src)
        do ii = 1, nlm_len    
            M_DDRX_src(ii,1:nlm_len) = matmul(GC(ii,:nlm_len,1:SHI_DDRX), g) ! len(g)=SHI_DDRX
        end do
    end
    
    function Gamma0(eps, T, A, Q)
        
        !-----------------------------------------------
        ! DDRX rate factor modelled as an Arrhenius activation process
        !----------------------------------------------- 
        
        ! eps = strain rate tensor 
        ! T   = temperature (deg. K)
        ! A   = rate factor magnitude
        ! Q   = activation energy
        
        implicit none

        real(kind=dp), intent(in) :: eps(3,3) ! strain-rate tensor
        real(kind=dp), intent(in) :: T, A, Q
        real(kind=dp)             :: Gamma0, epssq(3,3), I2
        real(kind=dp), parameter  :: gasconst = 8.31446261815324
        
        epssq = matmul(eps,eps)
        I2 = sqrt(0.5d0*(epssq(1,1)+epssq(2,2)+epssq(3,3))) ! = sqrt(0.5* eps:eps)
        Gamma0 = A * I2 * exp(-Q/(gasconst*T))
    end

    function M_CDRX()

        !--------------------------------------------
        ! Continuous dynamic recrystalization (CDRX)
        !--------------------------------------------
        
        ! Rotation recrystalization (polygonization) as a Laplacian diffusion process (Godert, 2003).
        ! Returns matrix M such that d/dt (nlm)_i = M_ij (nlm)_j
        ! NOTICE: This gives the unscaled effect of CDRX. The caller must multiply by an appropriate CDRX rate factor (scale) that should depend on temperature, stress, etc.

        implicit none
        real(kind=dp) :: M_CDRX(nlm_len,nlm_len)

        M_CDRX = 0.0
        do ii = 1, nlm_len 
            M_CDRX(ii,ii) = Ldiag(ii) ! Laplacian
        end do  
    end

    function M_REG(eps) 

        !----------------
        ! Regularization 
        !----------------
        
        ! Returns matrix M such that d/dt (nlm)_i = M_ij (nlm)_j
        ! NOTICE: Calibrated for 4 <= L <= 8 and L=20. For different L, the caller must specify an appropriate scaling.
        ! Calibration is provided by the script in tests/calibrate-regularization.
        
        implicit none
        real(kind=dp), intent(in) :: eps(3,3)
        real(kind=dp)             :: M_REG(nlm_len,nlm_len)
        real(kind=dp)             :: nu, expo, scalefac
        
        if (Lcap == 4) then
            expo = 1.55d0
            nu   = 2.9d0
        end if 

        if (Lcap == 6) then
!            expo = 2.000
!            nu   = 3.923652e+00
            expo = 1.25d0
            nu   = 3.6d0
        end if 
        
        if (Lcap == 8) then
            expo = 3.000
            nu   = 3.755041e+00
        end if 

!        if (Lcap == 8) then
!            expo = 1.500
!            nu   = 6.500000e+00
!        end if 

        if (Lcap == 10) then
            expo = 1.700
            nu   = 9e+0
        end if 
        
        if (Lcap == 20) then
            expo = 3.000
            nu   = 1.495197e+01
        end if 

        if ((Lcap .gt. 12) .and. (Lcap .lt. 20)) then
    !        print *, 'specfab error: returning the unscaled (but normalized) Laplacian matrix for you to scale yourself.'
            expo = 1
            scalefac = -1
        else
            scalefac = -nu * norm2(reshape(eps,[size(eps)])) 
        end if
        
        M_REG = 0.0
        do ii = 1, nlm_len 
            M_REG(ii,ii) = scalefac * abs( Ldiag(ii)/(Lcap*(Lcap+1)) )**expo 
        end do
    end

    !---------------------------------
    ! FABRIC DYNAMICS IN TENSORIAL SPACE
    !---------------------------------

    include "tensorialdynamics.f90"

    !---------------------------------
    ! BOUNDS AND POWER SPECTRUM
    !---------------------------------

    function apply_bounds(nlm) result (nlm_bounded)
    
        ! The entries of nlm are bounded in the sense that the corresponding angular power spectrum, S(l), must not exceed that of the delta function (perfect single maximum). 
        ! This function adjusts nlm *if needed* such that nlm is consistent with delta function bounds for components l=2,4 (for all m).
        ! For transient problems, this usually (combined with regularization) gives numerically stable and well-behaved fabric evolution.
        
        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        complex(kind=dp)             :: nlm_bounded(nlm_len)
        real(kind=dp)                :: S0, S2_rel, S4_rel 
        
        nlm_bounded = nlm
            
        S0 = real(nlm(1))**2 ! m=0 components are real by definition
        S2_rel = Sl(nlm,2)/S0
        S4_rel = Sl(nlm,4)/S0
        
        if (S2_rel > 1.0d0) then
            ! if S(2) > S_delta(2), then scale n2m such that S(2) = S_delta(2) (ensures a2 eigen values are correctly bounded: 0 <= lam1,lam2,lam3 <= 1)
            nlm_bounded(I_l2:(I_l4-1)) = nlm(I_l2:(I_l4-1)) / sqrt(S2_rel) 
        end if
        
        if (S4_rel > 1.0d0) then
            ! if S(4) > S_delta(4), then scale n4m such that S(4) = S_delta(4)
            nlm_bounded(I_l4:(I_l6-1)) = nlm(I_l4:(I_l6-1)) / sqrt(S4_rel)
        end if
    end

    function Sl(nlm, l)
    
        ! Angular power spectrum at degree l
        
        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        integer, intent(in)          :: l
        real(kind=dp)                :: Sl
        complex(kind=dp)             :: nlm_sub(2*l+1)
        integer                      :: I0, I1
        
        I1 = (l+1)*(l+2)/2      ! number of coefs if L=l
        I0 = I1 - (2*l+1) + 1   ! 2*l+1 coefs (possible m) for a given l
        nlm_sub(:) = nlm(I0:I1) ! nlm subrange with all "m" components for "l"
        Sl = 1.0d0/(2*l+1) * sum(abs(nlm_sub)**2) 
    end

    !---------------------------------
    ! QUADRICS
    !---------------------------------

    function quad_rr(M) result (q2m)
        
        ! Expansion coefficients of "M : rr" assuming M is symmetric and r is the radial unit vector
        
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
        
        ! Expansion coefficients of "M : tp" assuming M is anti-symmetric and (t,p) are the (theta,phi) unit vectors.
        
        implicit none
        
        real(kind=dp), intent(in) :: M(3,3) 
        real(kind=dp), parameter :: fsq1 = sqrt(2*Pi/3)
        complex(kind=dp) :: q1m(-1:1)
        
        ! components (1,-1), (1,0), (1,+1)
        q1m = fsq1*[r*M(y,z)-i*M(x,z), r*sqrt(2.)*M(x,y), -r*M(y,z)-i*M(x,z)]
    end
    
end module dynamics
