n2mhat(0:2) = n2m(0:2)/n00

ev(1,1) = -0.5d0*sqrt(2.0d0/3)*real(n2mhat(0)) + real(n2mhat(2)) 
ev(2,2) = -0.5d0*sqrt(2.0d0/3)*real(n2mhat(0)) - real(n2mhat(2)) 
ev(3,3) =       +sqrt(2.0d0/3)*real(n2mhat(0))

ev(1,2) = -aimag(n2mhat(2))
ev(2,1) = ev(1,2)

ev(1,3) = -real(n2mhat(1))
ev(3,1) = ev(1,3)

ev(2,3) = +aimag(n2mhat(1))
ev(3,2) = ev(2,3)

!k = f_ev_c0(n00) ! already normalized => k/f_ev_c0(n00) = 1
ev = sqrt(2/15.0d0)*ev + identity/3.0d0
