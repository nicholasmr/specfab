! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2024

! CPO dynamics in spectral space.

module dynamics  

    use header
    use tensorproducts
    use moments ! used by tensorial dynamics routines
    use gaunt

    implicit none 

    integer, parameter, private :: x = 1, y = 2, z = 3 ! Matrix indices (used for readability)
    complex(kind=dp), parameter, private :: r = (1,0), i = (0,1) ! real and imag units

    integer, private :: Lcap    ! Truncation "L" of expansion series (internal copy of what was passed to the init routine).
    integer, private :: nlm_len ! Total number of expansion coefficients (i.e. DOFs)

    ! Static (constant) matrices used for spectral dynamics
    real(kind=dp), parameter :: Ldiag(nlm_lenmax) = [( (-ll*(ll+1),mm=-ll,ll), ll=0,  Lcap__max,2)] ! Diagonal entries of Laplacian diffusion operator.
    
    integer, parameter :: SHI_LATROT=1+5
!    integer, parameter :: SHI_DDRX=1+5+9
    
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

        ! Lattice rotation modelled using plastic spin theory.
        ! Returns matrix M such that d/dt (nlm)_i = M_ij (nlm)_j

        implicit none

        real(kind=dp), intent(in) :: eps(3,3), omg(3,3), iota, zeta ! strain-rate (eps), spin (omg), eps^1 coefficient (iota), eps^2 coefficient (zeta)
        complex(kind=dp)          :: M_LROT(nlm_len,nlm_len), qe(-2:2), qo(-1:1)
        integer, parameter        :: SHI_LATROT = 1+5 ! Scope of harmonic interactions for LATROT
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
    
    function nlm_LROT(nlm0, dt, Nt, D,W, iota) result(nlm)
        ! Euler integrator of lattice rotation
        implicit none
        complex(kind=dp), intent(in) :: nlm0(nlm_len)
        integer, intent(in)          :: Nt
        real(kind=dp), intent(in)    :: dt, D(:,:,:), W(:,:,:), iota
        complex(kind=dp)             :: nlm(Nt,nlm_len)
        nlm(1,:) = nlm0 ! initial state
        do jj = 1,Nt-1 ! time step
            nlm(jj+1,:) = nlm(jj,:) + dt*matmul(M_LROT(D(jj,:,:), W(jj,:,:), iota, 0.0d0), nlm(jj,:)) ! Euler step
        end do
    end
    
    function dri_LROT(ri, D,W, iota) result(dri)
        ! Discrete version of lattice rotation
        implicit none
        real(kind=dp), intent(in) :: ri(:,:) ! axis no, xyz
        real(kind=dp), intent(in) :: D(3,3), W(3,3), iota
        real(kind=dp)             :: mm(3,3), Wp(3,3), dri(size(ri,1),size(ri,2))
        do jj = 1,size(ri,1) ! grain no
            mm = outerprod(ri(jj,:), ri(jj,:))
            Wp = iota*(matmul(mm,D) - matmul(D,mm))
            dri(jj,:) = matmul(W+Wp, ri(jj,:))
        end do
    end
    
    function ri_LROT(ri0, dt, Nt, D,W, iota) result(ri)
        ! Euler integrator of discrete version of lattice rotation
        implicit none
        real(kind=dp), intent(in) :: ri0(:,:) ! nn (crystallographic axis of nn'th grain), xyz
        integer, intent(in)       :: Nt
        real(kind=dp), intent(in) :: dt, D(:,:,:), W(:,:,:), iota
        real(kind=dp)             :: ri(Nt,size(ri0,1),size(ri0,2)), ri_new(size(ri0,1),size(ri0,2))
        ri(1,:,:) = ri0 ! initial state
        do ii = 1,Nt-1 ! time step
            ri_new = ri(ii,:,:) + dt*dri_LROT(ri(ii,:,:), D(ii,:,:), W(ii,:,:), iota) ! Euler step
            ri(ii+1,:,:) = ri_new/spread(norm2(ri_new,2), 2, 3) ! renormalize
        end do
    end
    
    !---- DEBUG NEW LATROT FORMULATION ----
    
    function plm_LROT(D,W, iota) result(plm)
        implicit none
        complex(kind=dp) :: plm(-2:2) ! entries correspond to lm pairs {(2,-2),(2,-1),(2,0),(2,1),(2,2)}
        real(kind=dp)    :: D(3,3),W(3,3), iota
        complex(kind=dp) :: WD(3,3), V(3,3)
        integer          :: ism(-2:2)
        real(kind=dp), parameter :: sr = sqrt(2.0d0/3)
        
        plm(:) = 0.0d0
        WD = W-iota*D
        
        do mm = -2,2 ! loop over m index
            ism(:) = 0
            if (mm == -2) ism(-2) = 1
            if (mm == -1) ism(-1) = 1
            if (mm == 0)  ism(0)  = 1
            if (mm == 1)  ism(1)  = 1
            if (mm == 2)  ism(2)  = 1
            
            V(1,1) =   -sr*ism(0) +ism(-2)+ism(2) ! xx
            V(2,2) =   -sr*ism(0) -ism(-2)-ism(2) ! yy
            V(3,3) = +2*sr*ism(0) ! zz
            V(2,3) = i*(ism(-1)+ism(1)) ! yz
            V(1,3) =   (ism(-1)-ism(1)) ! xz
            V(1,2) = i*(ism(-2)-ism(2)) ! xy
            V(3,2) = V(2,3)
            V(3,1) = V(1,3)
            V(2,1) = V(1,2)
            
            plm(mm) = sqrt(Pi/30.0d0) * doubleinner22complex(WD, V)
        end do
    end

    function qlm_LROT(D,W, iota) result(qlm)
        implicit none
        complex(kind=dp) :: qlm(-1:1) ! entries correspond to lm pairs {(1,-1), (1,0), (1,1)}
        real(kind=dp)    :: D(3,3),W(3,3), iota
        complex(kind=dp) :: WD(3,3), V(3,3)
        integer          :: ism(-1:1)
        
        qlm(:) = 0.0d0
        WD = W-iota*D
        
        do mm = -1,1 ! loop over m index
            ism(:) = 0
            if (mm == -1) ism(-1) = 1
            if (mm == 0)  ism(0)  = 1
            if (mm == 1)  ism(1)  = 1
            
            V(:,:) = 0.0d0
            V(2,3) =    -ism(-1)+ism(1)  ! yz
            V(1,3) = i*(+ism(-1)+ism(1)) ! xz
            V(1,2) = -sqrt(2.0d0)*ism(0) ! xy
            V(3,2) = -V(2,3)
            V(3,1) = -V(1,3)
            V(2,1) = -V(1,2)

            qlm(mm) = sqrt(Pi/6.0d0) * doubleinner22complex(WD, V)
        end do
    end
    
    function M_LROT_new(D, W, iota) result(M_LROT)
        implicit none
        real(kind=dp), intent(in) :: D(3,3), W(3,3), iota
        complex(kind=dp)          :: M_LROT(nlm_len,nlm_len)
        integer,parameter         :: SHI_A(5) = [4,5,6,7,8]
        integer,parameter         :: SHI_B(3) = [1,2,3]
        do ii = 1, nlm_len
!            M_LROT(ii,1:nlm_len) = matmul(GC_LROT_A(ii,1:nlm_len,1:5), plm_LROT(D,W,iota)) - i * matmul(GC_LROT_B(ii,1:nlm_len,1:3), qlm_LROT(D,W,iota))
            M_LROT(ii,1:nlm_len) = matmul(GC_LROT_A(ii,1:nlm_len,SHI_A), plm_LROT(D,W,iota)) - i * matmul(GC_LROT_B(ii,1:nlm_len,SHI_B), qlm_LROT(D,W,iota))
        end do
    end
    
    !---------------
    
    function plm_DDRX_C23(c0) result(plm)
        implicit none
        complex(kind=dp) :: plm(-1:1) ! entries correspond to lm pairs {(1,-1),(1,0),(1,1), (2,-2),(2,-1),(2,0),(2,1),(2,2)}
        real(kind=dp)    :: c0(3)
        complex(kind=dp) :: V(3)
        integer          :: ism(-1:1)
        
        plm(:) = 0.0d0
        
        do mm = -1,1 ! loop over m index
            ism(:) = 0
            if (mm == -1) ism(-1) = 1
            if (mm == 0)  ism(0)  = 1
            if (mm == 1)  ism(1)  = 1
            
            V(1) =    ism(-1) - ism(1)
            V(2) = i*(ism(-1) + ism(1))
            V(3) = sqrt(2.0d0)*ism(0)
            
            plm(mm) = sqrt(Pi*2.0d0/3) * dot_product(c0,V)
!            print *, plm
        end do
    end
    
    function M_DDRX_C23(c0) result(M_DDRX)
        implicit none
        real(kind=dp), intent(in) :: c0(3)
        complex(kind=dp)          :: M_DDRX(nlm_len,nlm_len)
        do ii = 1, nlm_len
            M_DDRX(ii,1:nlm_len) = matmul(GC_LROT_A(ii,1:nlm_len,1:3), plm_DDRX_C23(c0))  ! - i * matmul(GC_LROT_B(ii,1:nlm_len,1:3), qlm_LROT(D,W,iota))
        end do
    end
    
    !-----------------------------------------------

    function M_DDRX(nlm, tau)

        !-----------------------------------------------
        ! Discontinuous dynamic recrystallization (DDRX)
        !----------------------------------------------- 
        
        ! Nucleation and migration recrystalization modeled as a decay process (Placidi et al., 2010).
        ! Returns matrix M such that d/dt (nlm)_i = M_ij (nlm)_j
        
        ! NOTICE: This is Gamma/Gamma0. The caller must multiply by an appropriate DDRX rate factor, Gamma0(T,tau,eps,...).
        
        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        real(kind=dp), intent(in)    :: tau(3,3) ! Stress tensor
        complex(kind=dp)             :: M_DDRX(nlm_len,nlm_len)
        real(kind=dp)                :: Davg

        M_DDRX = M_DDRX_src(tau) ! Linear source term (D in matrix form)
        Davg   = ev_D2(nlm, tau) ! Nonlinear sink term (<D>)
        
        do ii = 1, nlm_len    
            M_DDRX(ii,ii) = M_DDRX(ii,ii) - Davg ! Add <D> to diagonal such that M = Gamma/Gamma0 = D - <D>*I
        end do
    end

    function M_DDRX_src(tau)  

        !------------------
        ! DDRX source term
        !------------------

        implicit none

        real(kind=dp), intent(in) :: tau(3,3) ! Stress tensor
        complex(kind=dp)          :: M_DDRX_src(nlm_len,nlm_len), qt(-2:2)
        integer, parameter        :: SHI_DDRX = 1+5+9 ! Scope of harmonic interactions for DDRX
        real(kind=dp)             :: k
        complex(kind=dp)          :: g(SHI_DDRX)

        qt = quad_rr(tau) ! Quadric expansion coefficients
        include "include/ddrx-coupling-weights.f90" ! Harmonic expansion coefficients of D (requires qt, sets g and k)
        g = k*g * 5/doubleinner22(tau,tau) ! Normalize by tau:tau/5
        
        do ii = 1, nlm_len    
            M_DDRX_src(ii,1:nlm_len) = matmul(GC(ii,:nlm_len,1:SHI_DDRX), g) ! len(g)=SHI_DDRX
        end do
    end
    
    function ev_D(nlm, tau, pow) result(D)
    
        !------------------
        ! Average basal-plane RSS to power "pow"
        !------------------
    
        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        real(kind=dp), intent(in)    :: tau(3,3) ! deviatoric stress tensor
        integer, intent(in)          :: pow
        real(kind=dp)                :: D

        D = 0.0d0
        
        if (pow == 2) then 
            D = ev_D2(nlm, tau)
        else if (pow == 4) then
            D = ev_D4(nlm, tau)
        end if
    end
    
    function ev_D2(nlm, tau) result(D)
    
        !------------------
        ! Average basal-plane RSS to *second* power (deformability "D", Placidi et al, 2010)
        !------------------
    
        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        real(kind=dp), intent(in)    :: tau(3,3) ! deviatoric stress tensor
        real(kind=dp)                :: D
        real(kind=dp)                :: tauv(6), tausq(3,3), norm, a4v(6,6), a2v(6)

        tauv = mat_to_vec(tau) ! 6x1 Mandel vector
        tausq = matmul(tau,tau) ! = tau.tau
        norm = tausq(1,1) + tausq(2,2) + tausq(3,3) ! tr(tau.tau) = tau:tau
        
        call f_ev_ck_Mandel(nlm, a2v, a4v) ! structure tensors in Mandel notation
        D = dot_product(mat_to_vec(tausq),a2v) - dot_product(tauv,matmul(a4v,tauv)) ! (tau.tau):a2 - tau:a4:tau
        D = 5 * D/norm ! normalize by basal RSS for isotropic polycrystal (tau:tau/5)
    end

    function ev_D4(nlm, tau) result(D)
    
        !------------------
        ! Average basal-plane RSS to the *forth* power
        !------------------
    
        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        real(kind=dp), intent(in)    :: tau(3,3) ! deviatoric stress tensor
        real(kind=dp)                :: D
        real(kind=dp)                :: norm, normsq, tausq(3,3), tautau(3,3,3,3), tausqtau(3,3,3,3)
        real(kind=dp)                :: a2mat(3,3), a4mat(3,3,3,3), a6mat(3,3,3,3, 3,3), a8mat(3,3,3,3, 3,3,3,3)

        tausq    = matmul(tau,tau)          ! = tau.tau
        tautau   = outerprodmat2(tau,tau)   ! = tau otimes tau
        tausqtau = outerprodmat2(tausq,tau) ! = (tau.tau) otimes tau
        norm = tausq(1,1) + tausq(2,2) + tausq(3,3) ! tr(tau.tau) = tau:tau
        normsq = norm**2
        
        call f_ev_ck(nlm, 'f', a2mat, a4mat, a6mat, a8mat) ! structure tensors 
        D = doubleinner22(tausq,doubleinner42(a4mat,tausq)) 
        D = D + doubleinner44(tautau,doubleinner84(a8mat,tautau))
        D = D - 2*doubleinner44(tausqtau,doubleinner62(a6mat,tau))
        D = 35/2.0d0 * D/normsq ! normalize by basal RSS for isotropic polycrystal
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
        I2 = sqrt(0.5d0*(epssq(1,1)+epssq(2,2)+epssq(3,3))) ! = sqrt(0.5* tr(eps.eps))
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

    function M_REG(D) 

        !----------------
        ! Regularization 
        !----------------
        
        ! Returns matrix M such that d/dt (nlm)_i = M_ij (nlm)_j
        ! Calibration constants (nu,expo) are calculated in tests/calibrate-regularization
        
        implicit none
        real(kind=dp), intent(in) :: D(3,3)
        real(kind=dp)             :: M_REG(nlm_len,nlm_len), M0_REG(nlm_len,nlm_len)
        real(kind=dp)             :: nu, expo, ratemag
        
        include "include/regcalib.f90"

        ! M0_REG is the matrix diag(~0, less small value, ..., 1) 
        M0_REG = 0.0
        do ii = 1, nlm_len 
            M0_REG(ii,ii) = abs( Ldiag(ii)/(Lcap*(Lcap+1)) )**expo ! Normalize by L*(L+1) to ensure =1 at modes l=L (note: this simply implies absorbing a constant into nu)
        end do
        
        ratemag = nu*norm2(reshape(D,[size(D)])) ! nu * ||D||_2
        M_REG = -ratemag*M0_REG
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
            nlm_bounded(L2rng) = nlm(L2rng) / sqrt(S2_rel) 
        end if
        
        if (S4_rel > 1.0d0) then
            ! if S(4) > S_delta(4), then scale n4m such that S(4) = S_delta(4)
            nlm_bounded(L4rng) = nlm(L4rng) / sqrt(S4_rel)
        end if
    end

    function Sl(nlm, l)
    
        ! Angular power spectrum at degree l
        
        implicit none

        complex(kind=dp), intent(in) :: nlm(nlm_len)
        integer, intent(in)          :: l
        real(kind=dp)                :: Sl
        complex(kind=dp)             :: nlm_sub(2*l+1)
        
        nlm_sub(:) = nlm(IL(l):(IL(l+2)-1)) ! [n_l^-l, ... n_l^l]
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
    
    !---------------------------------
    ! CORRELATIONS
    !---------------------------------
    
    function nhat40_empcorr_ice(nhat20) result(nhat40)
    
        ! Emperical correlation between nhat40=n40/n00 and nhat20=n20/n00 in a vertically symmetric frame, based on ice cores.
    
        implicit none
        
        real(kind=dp), intent(in) :: nhat20(:)
        real(kind=dp)             :: nhat40(size(nhat20))
        real(kind=dp), parameter  :: a=0.07283201, b=0.47585027, c=-0.27165286, d=0.14299372
        real(kind=dp), parameter  :: eps = 1.0d-4, lolim = -sqrt(5.0d0)/2, uplim = +sqrt(5.0d0)
        
        if (any(nhat20 > uplim+eps)) stop 'nhat40_empcorr_ice() error: nhat20 > sqrt(5)'
        if (any(nhat20 < lolim-eps)) stop 'nhat40_empcorr_ice() error: nhat20 < -sqrt(5)/2'

        do ii = 1, size(nhat20)
            if (nhat20(ii) > (uplim-eps)) then
                nhat40(ii) = 3.0 ! delta value
            else if (nhat20(ii) < (lolim+eps)) then
                nhat40(ii) = 1.1249998952803144 ! delta girdle value
            else
                nhat40(ii) = a*nhat20(ii) + b*nhat20(ii)**2 + c*nhat20(ii)**3 + d*nhat20(ii)**4
            end if
        end do

    end
    
end module dynamics
