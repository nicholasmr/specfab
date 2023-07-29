! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien, 2019-2023

! Tensor products

module tensorproducts 

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real, parameter, private    :: Pi = 3.1415927
    integer, private            :: ii,jj,kk,ll,mm,nn,oo,pp ! loop indices

contains      

    !------------------------
    ! GENERAL
    !------------------------

    function symmetricpart(A) result(Y)
        ! (A + A^T)/2
        implicit none
        real(kind=dp), intent(in) :: A(3,3)
        real(kind=dp) :: Y(3,3)
        Y = (A + transpose(A))/2
    end
    
    function antisymmetricpart(A) result(Y)
        ! (A - A^T)/2
        implicit none
        real(kind=dp), intent(in) :: A(3,3)
        real(kind=dp) :: Y(3,3)
        Y = (A - transpose(A))/2
    end
    
    function anticommutator(A,B) result(Y)
        ! A.B + B.A
        implicit none
        real(kind=dp), intent(in) :: A(3,3), B(3,3)
        real(kind=dp) :: Y(3,3)
        Y = matmul(A,B) + matmul(B,A)
    end    
    
    function diag(v) result(Y)
        ! create diagonal matrix Y(:,:) with diagonal elements v(:)
        implicit none
        real(kind=dp), intent(in) :: v(:)
        real(kind=dp)             :: Y(size(v),size(v))
        Y = 0.0d0
        do ii=1,size(v)
            Y(ii,ii) = v(ii)
        end do
    end

    !------------------------
    ! OUTER PRODUCTS
    !------------------------

    function outerprod(a,b) 
        ! a_i b_j = rank-2 tensor
        implicit none
        real(kind=dp), intent(in) :: a(3), b(3)
        real(kind=dp) :: outerprod(3,3)
        outerprod = reshape( [( [( a(ii)*b(jj), ii=1,3)], jj=1,3)] , [3,3])
    end
    
    function outerprod4(a,b,c,d) 
        ! a_i b_j c_k d_l = rank-4 tensor
        implicit none
        real(kind=dp), intent(in) :: a(3), b(3), c(3), d(3)
        real(kind=dp) :: outerprod4(3,3,3,3)
        outerprod4 = reshape( [( [( [( [( a(ii)*b(jj)*c(kk)*d(ll), ii=1,3)], jj=1,3)], kk=1,3)], ll=1,3)] , [3,3,3,3])
    end
    
    function outerprod6(a,b) 
        ! a_i b_j = rank-2 tensor
        implicit none
        real(kind=dp), intent(in) :: a(6), b(6)
        real(kind=dp) :: outerprod6(6,6)
        outerprod6 = reshape( [( [( a(ii)*b(jj), ii=1,6)], jj=1,6)] , [6,6])
    end
    
    function outerprod9(a,b) 
        ! a_i b_j = rank-2 tensor
        implicit none
        real(kind=dp), intent(in) :: a(9), b(9)
        real(kind=dp) :: outerprod9(9,9)
        outerprod9 = reshape( [( [( a(ii)*b(jj), ii=1,9)], jj=1,9)] , [9,9])
    end

    function outerprod22(A,B) result(C)
        ! A_ij B_kl = rank-4 tensor
        implicit none
        real(kind=dp), intent(in) :: A(3,3), B(3,3)
        real(kind=dp) :: C(3,3,3,3)
        C = reshape( [( [( [( [( A(ii,jj)*B(kk,ll), ii=1,3)], jj=1,3)], kk=1,3)], ll=1,3)] , [3,3,3,3])
    end    

    !------------------------
    ! INNER PRODUCTS
    !------------------------

    function singleinner13(x, F) result(Y)
        ! Y_jk = x_i F_ijk
        implicit none
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(in) :: F(size(x),3,3)
        real(kind=dp) :: Y(3,3)
        Y = 0.0d0
        do ii=1,size(x) 
            Y = Y + x(ii)*F(ii,:,:)
        end do
    end

    function doubleinner22(A,B) 
        ! A_ij B_ji = scalar
        implicit none
        real(kind=dp), intent(in) :: A(3,3), B(3,3)
        real(kind=dp) :: doubleinner22
        doubleinner22 = sum([( sum([(A(ii,jj)*B(jj,ii),jj=1,3)]) ,ii=1,3)])
    end

    function doubleinner42(A,B) 
        ! A_lkij B_ji = rank-2 tensor
        implicit none
        real(kind=dp), intent(in) :: A(3,3,3,3), B(3,3)
        real(kind=dp) :: doubleinner42(3,3)
        do ll=1,3
            do kk=1,3
                doubleinner42(ll,kk) = doubleinner22(A(ll,kk,:,:),B) 
            end do
        end do
    end

    function doubleinner24(B,A) result(C)
        ! B_ij A_jilk = rank-2 tensor
        implicit none
        real(kind=dp), intent(in) :: A(3,3,3,3), B(3,3)
        real(kind=dp) :: C(3,3)
        do ll=1,3
            do kk=1,3
                C(ll,kk) = doubleinner22(B,A(:,:,ll,kk)) 
            end do
        end do
    end

    function doubleinner42_firstlast(A,B) result(C)
        ! A_lkij B_lj = rank-2 tensor
        implicit none
        real(kind=dp), intent(in) :: A(3,3,3,3), B(3,3)
        real(kind=dp) :: C(3,3)
        do kk=1,3
            do ii=1,3
                C(kk,ii) = doubleinner22(A(:,kk,ii,:),B)
            end do
        end do
    end

    function doubleinner62(A,B) 
        ! A_nmlkij B_ji = rank-4 tensor
        implicit none
        real(kind=dp), intent(in) :: A(3,3,3,3, 3,3), B(3,3)
        real(kind=dp) :: doubleinner62(3,3,3,3)
        do nn=1,3
            do mm=1,3
                do ll=1,3
                    do kk=1,3
                        doubleinner62(nn,mm,ll,kk) = doubleinner22(A(nn,mm,ll,kk,:,:),B) 
                    end do
                end do
            end do
        end do
    end

    function doubleinner82(A,B) 
        ! A_ponmlkij B_ji = rank-6 tensor
        implicit none
        real(kind=dp), intent(in) :: A(3,3,3,3, 3,3,3,3), B(3,3)
        real(kind=dp) :: doubleinner82(3,3,3,3, 3,3)
        do pp=1,3
            do oo=1,3
                do nn=1,3
                    do mm=1,3
                        do ll=1,3
                            do kk=1,3
                                doubleinner82(pp,oo,nn,mm,ll,kk) = doubleinner22(A(pp,oo,nn,mm,ll,kk,:,:),B) 
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end
    
end module tensorproducts
