! N. M. Rathmann <rathmann@nbi.ku.dk>

! Fabric eigenframe etc.

module frames

    use header
    use moments

    implicit none 
    
contains

    subroutine eig(nlm, ei,lami)

        implicit none
        integer, parameter           :: n = 3
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp), intent(out)   :: ei(n,n), lami(n) ! (i,xyz), (i) for eigenpair i=1,2,3

        call eigframe(a2(nlm), 'ij', ei,lami) ! shorthand for eigenframe of nlm
    end

    subroutine eigframe(M,plane, ei,lami)

        implicit none
        integer, parameter         :: n = 3
        real(kind=dp), intent(in)  :: M(n,n)
        character*2, intent(in)    :: plane ! ['ij'|'xy'|'xz']
        real(kind=dp), intent(out) :: ei(n,n), lami(n) ! (i,xyz), (i) for eigenpair i=1,2,3
        
        real(kind=dp)              :: e1(n),e2(n),e3(n)
        integer                    :: k, sort(n), Imax
        
        call eig3(M, e1,e2,e3, lami) ! largest eigenvalue pair first
        ei(1,:) = e1
        ei(2,:) = e2
        ei(3,:) = e3
        
        if (plane == 'ij') then
            sort = [1,2,3] ! no sorting required
        else
            if (plane == 'xy') then
                k = 3 ! z is the out-of-plane direction 
            else if (plane == 'xz') then
                k = 2 ! y is the out-of-plane direction
            else
                stop 'specfab error: eigframe(M, plane) argument "plane" was not any of ij,xy,xz'
            end if

            Imax = maxloc(abs(ei(:,k)), 1)
            if (Imax==1) sort = [2,3,1] ! e1 is out-of-plane
            if (Imax==2) sort = [1,3,2] ! e2 is out-of-plane
            if (Imax==3) sort = [1,2,3] ! e3 is out-of-plane
        end if
        
        ! Apply sorting
        ei(:,:) = ei(sort,:)
        lami(:) = lami(sort)
    end
    
    subroutine eig3(M, e1,e2,e3, eigvals)
        
        ! Eigenvalue pair of symmetric real-valued 3x3 matrix

        implicit none
        integer, parameter           :: n = 3
        real(kind=dp), intent(in)    :: M(n,n)
        real(kind=dp), intent(out)   :: e1(n),e2(n),e3(n), eigvals(3)

        integer            :: inf
        integer, parameter :: l=3*3-1
        real(kind=dp)      :: work(l)
        
        call dsyev('V','U',n,M,n,eigvals,work,l,inf)
        e1 = M(:,3)
        e2 = M(:,2)
        e3 = M(:,1)
        eigvals = [eigvals(3),eigvals(2),eigvals(1)] ! Largest first        
    end
    
!    subroutine frame(nlm, ftype, e1,e2,e3, eigvals)

!        implicit none
!        complex(kind=dp), intent(in) :: nlm(:)
!        character*1, intent(in)      :: ftype ! 'x','e','p' (Cartesian frame, a2 eigenframe, 45deg rotated eigenframe)
!        real(kind=dp), intent(out)   :: e1(3),e2(3),e3(3), eigvals(3)
!        real(kind=dp)                :: p(3),q(3)
!        
!        ! Cartesian frame
!        if (ftype == 'x') then
!            e1 = [1,0,0] 
!            e2 = [0,1,0] 
!            e3 = [0,0,1] 
!            
!        else
!            ! a2 eigen frame    
!            call eig3(a2(nlm), e1,e2,e3, eigvals)
!            
!            ! Rotated a2 eigen frame
!            if (ftype == 'p') then
!                p = (e1+e2)/sqrt(2.0) 
!                q = (e1-e2)/sqrt(2.0) 
!                e1 = p
!                e2 = q
!                ! cross product
!                e3(1) = p(2) * q(3) - p(3) * q(2)
!                e3(2) = p(3) * q(1) - p(1) * q(3)
!                e3(3) = p(1) * q(2) - p(2) * q(1)
!            end if     
!               
!        end if
!    end
    
    function pqframe(ei) result(pq)

        ! Rotated a2 eigenframe

        implicit none
        real(kind=dp), intent(in) :: ei(3,3)
        real(kind=dp)             :: pq(3,3)

        pq(1,:) = (ei(1,:)+ei(2,:))/sqrt(2.0) 
        pq(2,:) = (ei(1,:)-ei(2,:))/sqrt(2.0) 
        ! cross product
        pq(3,1) = pq(1,2) * pq(2,3) - pq(1,3) * pq(2,2)
        pq(3,2) = pq(1,3) * pq(2,1) - pq(1,1) * pq(2,3)
        pq(3,3) = pq(1,1) * pq(2,2) - pq(1,2) * pq(2,1)
    end

    subroutine a4_eigentensors(nlm, Q1,Q2,Q3,Q4,Q5,Q6, eigvals6)  

        implicit none        
        complex(kind=dp), intent(in) :: nlm(:)
        real(kind=dp)                :: M(6,6)
        real(kind=dp), intent(out)   :: Q1(3,3),Q2(3,3),Q3(3,3),Q4(3,3),Q5(3,3),Q6(3,3)
        
        integer, parameter :: n = 6
        real(kind=dp), intent(out) :: eigvals6(n)
        integer            :: inf
        integer, parameter :: l6=n*n-1
        real(kind=dp)      :: work6(l6)
        
        M = a4_to_mat(a4(nlm)) ! 6x6
        call dsyev('V','U',n,M,n,eigvals6,work6,l6,inf)
        Q1 = vec_to_mat(M(:,1)) 
        Q2 = vec_to_mat(M(:,2)) 
        Q3 = vec_to_mat(M(:,3)) 
        Q4 = vec_to_mat(M(:,4))
        Q5 = vec_to_mat(M(:,5))
        Q6 = vec_to_mat(M(:,6))
    end
    
end module frames
