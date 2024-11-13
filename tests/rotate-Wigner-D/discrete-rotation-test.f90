! N. M. Rathmann <rathmann@nbi.ku.dk>, 2022

program demo

    use specfab    
    
    implicit none

!    integer, parameter :: dp = 8
    
    real(kind=dp), dimension(3), parameter :: x1 = [1,0,0], x2 = [0,1,0], x3 = [0,0,1] ! x,y,z dir.
    real(kind=dp) :: theta = 0.0d0, phi = 0.0d0
    
    call initspecfab(4) 
   
    theta = 1*(+pi/2)
    phi   = 0*(+pi/2)

    print *, 'x1=', x1
    print *, rotate_vector(x1, theta, phi)
   
    print *, 'x2=', x2
    print *, rotate_vector(x2, theta, phi)
    
end program
