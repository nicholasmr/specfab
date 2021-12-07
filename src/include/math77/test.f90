program drdnqsol

    use lambdasolver

    integer :: nglen
    double precision :: Eii(3), lami(3)
    
    nglen = 3
    Eii = [5,10,20]
    
    lami = lambdaprime(nglen, Eii)
    print *, lami
    
end program
