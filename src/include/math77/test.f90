program drdnqsol

    use lambdasolver

    double precision :: nglen, Eii(3), lami(3)
    
    nglen = 3.0d0
    Eii = [5,10,20]
    
    lami = lambdaprime(nglen, Eii)
    print *, lami
    
end program
