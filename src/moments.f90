! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-

! Normalized vector moments (structure tensors) given the spectral expansion coefficients of CPO.

module moments 

    use header
    use mandel
    use tensorproducts

    implicit none 
    
    ! Set by init()
    integer, private :: Lcap ! Truncation "L" of expansion series 
    integer, private :: nlm_len
    
    real(kind=dp), parameter :: a2_iso(3,3)     = identity/3
    real(kind=dp), parameter :: a4_iso(3,3,3,3) = reshape([ ( ( ( ( (identity(ii,jj)*identity(kk,ll)+identity(ii,kk)*identity(jj,ll)+identity(ii,ll)*identity(jj,kk))/15, ii=1,3), jj=1,3), kk=1,3), ll=1,3) ], [3,3,3,3]) 
    
contains      

    subroutine initmoments(Lcap_) 

        ! Needs to be called once before using the module routines.

        implicit none    
        integer, intent(in) :: Lcap_ ! Truncation "Lcap"
        
        Lcap = Lcap_ ! Save internal copy
        nlm_len = nlm_lenvec(Lcap)
    end

    !---------------------------------
    ! Conversion between nlm and <c^k> 
    !---------------------------------
              
    function a2(nlm) 
        ! a^(2) := <c^2> 
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a2(3,3)

        a2 = f_ev_c2(nlm(L0rng),nlm(L2rng))
    end

    function a4(nlm) 
        ! a^(4) := <c^4> 
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a4(3,3,3,3)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)

        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        a4 = f_ev_c4(n00,n2m,n4m)
    end
    
    function a6(nlm) 
        ! a^(6) := <c^6> 
        implicit none
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: a6(3,3,3,3,3,3)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        a6 = f_ev_c6(n00,n2m,n4m,n6m)
    end
    
    function a2_to_nlm(a2) result(nlm)
        implicit none
        real(kind=dp), intent(in) :: a2(3,3)
        complex(kind=dp)          :: nlm(nlm_lenvec(2)) 
        nlm = (0.0d0, 0.0d0) ! init
        include "include/a2_to_nlm__body.f90"
    end

    function a4_to_nlm(a4) result(nlm)
        implicit none
        real(kind=dp), intent(in) :: a4(3,3,3,3)
        complex(kind=dp)          :: nlm(nlm_lenvec(4)) 
        real(kind=dp)             :: a4Mandel(6,6)        
        nlm = 0.0d0 ! init
        a4Mandel = a4_to_mat(a4) ! 6x6 Mandel matrix of a4
        include "include/a4_to_nlm__body.f90"
    end

    function a6_to_nlm(a6) result(nlm)
        implicit none
        real(kind=dp), intent(in) :: a6(3,3,3,3,3,3)
        complex(kind=dp)          :: nlm(nlm_lenvec(6)) 
        nlm = 0.0d0 ! init
        include "include/a6_to_nlm__body.f90"
    end

    !---------------------------------
    ! Conversion from ri to nlm
    !---------------------------------

    function ri_to_nlm(ri, wi, L) result(nlm)
    
        implicit none
        real(kind=dp), intent(in)     :: ri(:,:), wi(:) ! (grain no, xyz), (grain no)
        integer, intent(in)           :: L
        integer                       :: N
        real(kind=dp)                 :: a2mat(3,3), a4mat(3,3,3,3), a6mat(3,3,3,3,3,3)
        complex(kind=dp)              :: nlm(nlm_lenvec(L)) 
    
        N = size(ri,1) ! number of grains
        nlm = 0.0d0
    
        if (L == 2) then
            a2mat = 0.0d0
            do ii=1,N
                a2mat = a2mat + wi(ii)*outerprod(ri(ii,:),ri(ii,:))
            end do
            nlm(1:L2len) = a2_to_nlm(a2mat)
            
        elseif (L == 4) then
            a4mat = 0.0d0
            do ii=1,N
                a4mat = a4mat + wi(ii)*outerprod4(ri(ii,:),ri(ii,:),ri(ii,:),ri(ii,:))
            end do
            nlm(1:L4len) = a4_to_nlm(a4mat)
            
        else
            a6mat = 0.0d0
            do ii=1,N
                a6mat = a6mat + wi(ii)*outerprod6(ri(ii,:),ri(ii,:),ri(ii,:),ri(ii,:),ri(ii,:),ri(ii,:))
            end do
            nlm(1:L6len) = a6_to_nlm(a6mat)
        end if
    end

    !---------------------------------
    ! All moments <c^k> (k = 2,4,6,8)
    !---------------------------------
    
    subroutine f_ev_ck(nlm, opt, ev_c2,ev_c4,ev_c6,ev_c8)
        
        ! "ev_ck" are the structure tensors <c^k> := a^(k) for a given n(theta,phi) prescribed in terms of "nlm" expansion coefficients.
        
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        character*1, intent(in)      :: opt 
        real(kind=dp), intent(inout) :: ev_c2(3,3),ev_c4(3,3,3,3),ev_c6(3,3,3,3, 3,3),ev_c8(3,3,3,3, 3,3,3,3)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6), n8m(-8:8)
        
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        ev_c2 = f_ev_c2(n00,n2m)
        ev_c4 = f_ev_c4(n00,n2m,n4m)
        
        if (opt == 'f') then
            ! Full calculation?
            n8m = nlm(L8rng)
            ev_c6 = f_ev_c6(n00,n2m,n4m,n6m)
            ev_c8 = f_ev_c8(n00,n2m,n4m,n6m,n8m)
        else
            ! Reduced calculation? Higher-order contributions not needed by caller.
            ev_c6 = 0.0d0
            ev_c8 = 0.0d0
        end if
    end
    
    subroutine f_ev_ck_Mandel(nlm, ev_c2_Mandel, ev_c4_Mandel)
            
        ! a2 and a4 in Mandel's notation. 
        ! Mandel's notation should be used wherever possible for performance.
            
        implicit none
        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(out)   :: ev_c2_Mandel(6),ev_c4_Mandel(6,6)
        complex(kind=dp)             :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        ev_c2_Mandel = mat_to_vec(f_ev_c2(n00,n2m)) ! 6x1 Mandel vector of a2
        ev_c4_Mandel = f_ev_c4_Mandel(n00,n2m,n4m)  ! 6x6 Mandel matrix of a4
    end

    !---------------------------------
    ! Individual <c^k> moments
    !---------------------------------

    function f_ev_c0(n00) result(ev)
        ! Integral over orientation distribution = total number of c-axes = N
        complex(kind=dp), intent(in) :: n00
        real(kind=dp) :: ev
        ev = REAL(sqrt(4*Pi)*n00)
    end

    function f_ev_c2(n00,n2m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2)
        complex(kind=dp) :: n2mhat(0:2)
!        real(kind=dp) :: k = 0.0
        real(kind=dp) :: ev(3,3)
        ev = 0.0d0
        include "include/ev_c2__body.f90"
!        ev = ev * k/f_ev_c0(n00) ! not needed, <c^2> already normalized
    end
    
    function f_ev_c4(n00,n2m,n4m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4)
        real(kind=dp) :: k = 0.0d0, ev(3,3, 3,3)
        ev = 0.0d0
        include "include/ev_c4__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c4_Mandel(n00,n2m,n4m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4)
        real(kind=dp) :: k = 0.0d0, ev(6,6)
        ev = 0.0d0
        include "include/ev_c4_Mandel__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c6(n00,n2m,n4m,n6m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)
        real(kind=dp) :: k = 0.0d0, ev(3,3, 3,3, 3,3)
        ev = 0.0d0
        include "include/ev_c6__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end
    
    function f_ev_c8(n00,n2m,n4m,n6m,n8m) result(ev)
        implicit none
        complex(kind=dp), intent(in) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6), n8m(-8:8)
        real(kind=dp) :: k = 0.0d0, ev(3,3, 3,3, 3,3, 3,3)
        ev = 0.0d0
        include "include/ev_c8__body.f90"
        ev = ev * k/f_ev_c0(n00)
    end

    !---------------------------------
    ! Orthogonal distribution
    !---------------------------------
    
    function a2_orth(blm,nlm) result(ev)
    
        ! Calculates <v^2> where v = b x n
        ! Identically includes xi=sin^2(alpha) kernel due to cross product

        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: k=0.0, norm=0.0, ev(3,3)
        complex(kind=dp)             :: b00,n00, b2m(-2:2),n2m(-2:2), b4m(-4:4),n4m(-4:4), b6m(-6:6),n6m(-6:6)
        
        call decompose_nlm(blm, b00,b2m,b4m,b6m)
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        
        ev = 0.0
        include "include/ev_v2__body.f90"
        ev = ev * k/norm
    end
    
    function a4_orth(blm,nlm) result(ev)
    
        ! Calculates <v^4> where v = b x n
        ! Identically includes xi=sin^4(alpha) kernel due to cross product

        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: k=0.0, norm=0.0, ev(3,3,3,3)
        complex(kind=dp)             :: b00,n00, b2m(-2:2),n2m(-2:2), b4m(-4:4),n4m(-4:4), b6m(-6:6),n6m(-6:6)

        call decompose_nlm(blm, b00,b2m,b4m,b6m)
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        
        ev = 0.0
        include "include/ev_v4__body.f90"
        ev = ev * k/norm 
    end
  
    function a4_joint(blm,nlm) result(ev)
    
        ! Calculates < b^2 n^2 xi(r,r') > where kernel xi=sin^2(alpha) where alpha = angle between n and b vectors

        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: k=0.0, norm=0.0, ev(3,3,3,3)
        complex(kind=dp)             :: b00,n00, b2m(-2:2),n2m(-2:2), b4m(-4:4),n4m(-4:4), b6m(-6:6),n6m(-6:6)

        call decompose_nlm(blm, b00,b2m,b4m,b6m)
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        
        ev = 0.0
        include "include/ev_c2b2__body.f90"
        ev = ev * k/norm 
    end
    
    function a4_jointcross(blm,nlm) result(ev)
    
        ! Calculates <n^2 v^2> where v = b x n 
        ! Identically includes xi=sin^2(alpha) kernel due to cross product

        implicit none
        complex(kind=dp), intent(in) :: blm(:), nlm(:)
        real(kind=dp)                :: k=0.0, norm=0.0, ev(3,3,3,3)
        complex(kind=dp)             :: b00,n00, b2m(-2:2),n2m(-2:2), b4m(-4:4),n4m(-4:4), b6m(-6:6),n6m(-6:6)

        call decompose_nlm(blm, b00,b2m,b4m,b6m)
        call decompose_nlm(nlm, n00,n2m,n4m,n6m)
        
        ev = 0.0
        include "include/ev_c2v2__body.f90"
        ev = ev * k/norm 
    end
    
    !---------------------------------
    ! Irreducible parts
    !---------------------------------
    
    function a2_irpart(a2_)
        ! Irreducible part of a2 (Zheng and Zou, 2001)
        implicit none
        real(kind=dp), intent(in) :: a2_(3,3)
        real(kind=dp)             :: a2_irpart(3,3)
        
        a2_irpart = a2_ - a2_iso ! a2 - I/3
    end
    
    function a4_irpart(a2_,a4_)
        ! Irreducible part of a4 (Zheng and Zou, 2001)
        implicit none
        real(kind=dp), intent(in) :: a2_(3,3), a4_(3,3,3,3)
        real(kind=dp)             :: a4_irpart(3,3,3,3), a2I(3,3,3,3)
        
        a2I = reshape([ ( ( ( ( identity(ii,jj)*a2_(kk,ll) + identity(ii,kk)*a2_(jj,ll) + identity(ii,ll)*a2_(jj,kk) &
                              + identity(jj,kk)*a2_(ii,ll) + identity(jj,ll)*a2_(ii,kk) + identity(kk,ll)*a2_(ii,jj) , ii=1,3), jj=1,3), kk=1,3), ll=1,3) ], [3,3,3,3]) 
        a4_irpart = a4_ - 1.0d0/7*a2I + 15.0d0/35*a4_iso 
    end
    
    !---------------------------------
    ! AUX
    !---------------------------------
    
    subroutine decompose_nlm(nlm, n00,n2m,n4m,n6m)

        implicit none
        
        complex(kind=dp), intent(in)  :: nlm(:)
        complex(kind=dp), intent(out) :: n00, n2m(-2:2), n4m(-4:4), n6m(-6:6)

        n4m = 0.0d0
        n6m = 0.0d0
        
        n00 = nlm(L0rng)
        n2m = nlm(L2rng)
        if (size(nlm) >= nlm_lenvec(4)) n4m = nlm(L4rng)
        if (size(nlm) >= nlm_lenvec(6)) n6m = nlm(L6rng)
    end
    
    subroutine ai_orthotropic(qlm_1,qlm_2,qlm_3, a2_i, a4_ii, a4_jk)
    
        implicit none
        
        complex(kind=dp), intent(in) :: qlm_1(:), qlm_2(:), qlm_3(:) ! distributions of of slip-system axes qi = (b,n,v), respectively
        real(kind=dp), intent(out)   :: a2_i(3, 3,3), a4_ii(3, 3,3,3,3), a4_jk(3, 3,3,3,3) ! a4_jk = <m_j m_j m_k m_k>

        a2_i(1, :,:)      = a2(qlm_1) ! <q1^2>
        a2_i(2, :,:)      = a2(qlm_2) ! <q2^2>
        a4_ii(1, :,:,:,:) = a4(qlm_1) ! <q1^4>
        a4_ii(2, :,:,:,:) = a4(qlm_2) ! <q2^4>
        a4_jk(3, :,:,:,:) = a4_joint(qlm_1,qlm_2) ! <q1^2 q2^2>

        ! Is qlm_3 given? Use it
        if (real(qlm_3(1)) > 1e-8) then
            a2_i(3, :,:)      = a2(qlm_3) ! <q3^2>
            a4_ii(3, :,:,:,:) = a4(qlm_3) ! <q3^4>
            a4_jk(1, :,:,:,:) = a4_joint(qlm_2,qlm_3) ! <q2^2 q3^2>
            a4_jk(2, :,:,:,:) = a4_joint(qlm_1,qlm_3) ! <q1^2 q3^2>
            
        ! else derive qlm_3 from qlm_1, qlm_2 (i.e. v from n,b distributions)
        else
            a2_i(3, :,:)      = a2_orth(qlm_1, qlm_2) ! <q3^2> 
            a4_ii(3, :,:,:,:) = a4_orth(qlm_1, qlm_2) ! <q3^4>
            a4_jk(1, :,:,:,:) = a4_jointcross(qlm_2,qlm_1) ! <q2^2 q3^2> 
            a4_jk(2, :,:,:,:) = a4_jointcross(qlm_1,qlm_2) ! <q1^2 q3^2> 
        end if
    end
    
    subroutine ai_orthotropic_discrete(qi, a2_i, a4_ii, a4_jk)

        implicit none
        
        real(kind=dp), intent(in)  :: qi(:,:,:) ! (3,3,N) = (m'_i, xyz comp., grain no.) 
        real(kind=dp), intent(out) :: a2_i(3, 3,3), a4_ii(3, 3,3,3,3), a4_jk(3, 3,3,3,3)
        real(kind=dp) :: q1(3), q2(3), q3(3)
        integer       :: ii, nn, N

        a2_i  = 0.0d0
        a4_ii = 0.0d0
        a4_jk = 0.0d0
        
        N = size(qi, dim=3) ! number of grains
        
        do nn=1,N
        
            do ii=1,3
                q1 = qi(ii,:,nn) ! use q1 as dummy
                a2_i(ii, :,:)      = a2_i(ii, :,:)      + outerprod(q1,q1)/N
                a4_ii(ii, :,:,:,:) = a4_ii(ii, :,:,:,:) + outerprod4(q1,q1,q1,q1)/N
            end do 
            
            q1 = qi(1,:,nn)
            q2 = qi(2,:,nn)
            q3 = qi(3,:,nn)
            
            a4_jk(1, :,:,:,:) = a4_jk(1, :,:,:,:) + outerprod4(q2,q2, q3,q3)/N
            a4_jk(2, :,:,:,:) = a4_jk(2, :,:,:,:) + outerprod4(q1,q1, q3,q3)/N
            a4_jk(3, :,:,:,:) = a4_jk(3, :,:,:,:) + outerprod4(q1,q1, q2,q2)/N
        end do
    end
    
end module moments
