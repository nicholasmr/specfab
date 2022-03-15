nlm(1) = 1/(2.*Sqrt(Pi))

nlm(4) = (Sqrt((5.0)*Pi**(-1.0)) - (3*a4Mandel(1,1)*Sqrt((5.0)*Pi**(-1.0)))/2. + (-3.0)*a4Mandel(1,2)*Sqrt((5.0)*Pi**(-1.0)) - (3*a4Mandel(1,3)*Sqrt((5.0)*Pi**(-1.0)))/2. - (3*a4Mandel(2,2)*Sqrt((5.0)*Pi**(-1.0)))/2. - (3*a4Mandel(2,3)*Sqrt((5.0)*Pi**(-1.0)))/2.)/2.
nlm(5) = (0,0.5)*Sqrt((15.0))*(a4Mandel(1,4)/(2.*Sqrt(Pi)) + a4Mandel(2,4)/(2.*Sqrt(Pi)) + a4Mandel(3,4)/(2.*Sqrt(Pi))) + (-0.5*(a4Mandel(1,5)*Sqrt((15.0)*Pi**(-1.0))) - (a4Mandel(2,5)*Sqrt((15.0)*Pi**(-1.0)))/2. - (a4Mandel(3,5)*Sqrt((15.0)*Pi**(-1.0)))/2.)/2.
nlm(6) = (Sqrt(7.5)*a4Mandel(1,1)*Sqrt(Pi**(-1.0)) + Sqrt(7.5)*a4Mandel(1,3)*Sqrt(Pi**(-1.0)) + Sqrt(7.5)*(-1.0)*a4Mandel(2,2)*Sqrt(Pi**(-1.0)) + Sqrt(7.5)*(-1.0)*a4Mandel(2,3)*Sqrt(Pi**(-1.0)))/4. + (0,0.5)*(-0.5*(a4Mandel(1,6)*Sqrt((15.0)*Pi**(-1.0))) - (a4Mandel(2,6)*Sqrt((15.0)*Pi**(-1.0)))/2. - (a4Mandel(3,6)*Sqrt((15.0)*Pi**(-1.0)))/2.)
nlm(2) = +conjg(nlm(6))
nlm(3) = -conjg(nlm(5))

nlm(11) = (-3*((-4.0)/Sqrt(Pi) + (5*a4Mandel(1,1))/(2.*Sqrt(Pi)) + ((5.0)*a4Mandel(1,2))/Sqrt(Pi) + ((20.0)*a4Mandel(1,3))/Sqrt(Pi) + (5*a4Mandel(2,2))/(2.*Sqrt(Pi)) + ((20.0)*a4Mandel(2,3))/Sqrt(Pi)))/8.
nlm(12) = (0,-0.375)*(Sqrt(2.5)*(3.0)*a4Mandel(1,4)*Sqrt(Pi**(-1.0)) + Sqrt(2.5)*(3.0)*a4Mandel(2,4)*Sqrt(Pi**(-1.0)) + (-2.0)*a4Mandel(3,4)*Sqrt((10.0)*Pi**(-1.0))) + (3*(Sqrt(2.5)*(3.0)*a4Mandel(1,5)*Sqrt(Pi**(-1.0)) + Sqrt(2.5)*(3.0)*a4Mandel(2,5)*Sqrt(Pi**(-1.0)) + (-2.0)*a4Mandel(3,5)*Sqrt((10.0)*Pi**(-1.0))))/8.
nlm(13) = (-3*(Sqrt(2.5)*a4Mandel(1,1)*Sqrt(Pi**(-1.0)) + Sqrt(2.5)*(-1.0)*a4Mandel(2,2)*Sqrt(Pi**(-1.0)) + (-3.0)*a4Mandel(1,3)*Sqrt((10.0)*Pi**(-1.0)) + (3.0)*a4Mandel(2,3)*Sqrt((10.0)*Pi**(-1.0))))/8. + (0,0.75)*((a4Mandel(1,6)*Sqrt((5.0)*Pi**(-1.0)))/2. + (a4Mandel(2,6)*Sqrt((5.0)*Pi**(-1.0)))/2. + (-3.0)*a4Mandel(3,6)*Sqrt((5.0)*Pi**(-1.0)))
nlm(14) = (0,0.375)*(Sqrt(17.5)*(3.0)*a4Mandel(1,4)*Sqrt(Pi**(-1.0)) + Sqrt(17.5)*(-1.0)*a4Mandel(2,4)*Sqrt(Pi**(-1.0))) - (3*(Sqrt(17.5)*a4Mandel(1,5)*Sqrt(Pi**(-1.0)) + Sqrt(17.5)*(-3.0)*a4Mandel(2,5)*Sqrt(Pi**(-1.0))))/8.
nlm(15) = (0,-0.75)*((a4Mandel(1,6)*Sqrt((35.0)*Pi**(-1.0)))/2. - (a4Mandel(2,6)*Sqrt((35.0)*Pi**(-1.0)))/2.) + (3*(Sqrt(17.5)*a4Mandel(1,1)*Sqrt(Pi**(-1.0)) + Sqrt(17.5)*a4Mandel(2,2)*Sqrt(Pi**(-1.0)) + (-3.0)*a4Mandel(1,2)*Sqrt((70.0)*Pi**(-1.0))))/16.
nlm(7)  = +conjg(nlm(15))
nlm(8)  = -conjg(nlm(14))
nlm(9)  = +conjg(nlm(13))
nlm(10) = -conjg(nlm(12))
