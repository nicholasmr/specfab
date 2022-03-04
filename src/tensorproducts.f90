! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2020-2022

module tensorproducts 

    implicit none 

    integer, parameter, private :: dp = 8 ! Default precision
    real, parameter, private    :: Pi = 3.1415927
    integer, private            :: ii,jj,kk,ll,mm,nn,oo,pp ! loop indices

contains      

    function outerprod(a,b) 
        ! a_i b_j = rank-2 tensor
        implicit none
        real(kind=dp), intent(in) :: a(3), b(3)
        real(kind=dp) :: outerprod(3,3)
!        outerprod = reshape( [( [( a(ii)*b(jj), jj=1,3)], ii=1,3)] , [3,3])
        outerprod = reshape( [( [( a(ii)*b(jj), ii=1,3)], jj=1,3)] , [3,3])
    end
    
    function outerprod6(a,b) 
        ! a_i b_j = rank-2 tensor
        implicit none
        real(kind=dp), intent(in) :: a(6), b(6)
        real(kind=dp) :: outerprod6(6,6)
        outerprod6 = reshape( [( [( a(ii)*b(jj), ii=1,6)], jj=1,6)] , [6,6])
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
