function Y00(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=1/(2.*Sqrt(Pi))
end

function Y20(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(Sqrt((5.0d0)*Pi**(-1.0d0))*((-1.0d0) + (3.0d0)*Cos(theta)**(2.0d0)))/4.
end

function Y40(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(3*((3.0d0) + (-30.0d0)*Cos(theta)**(2.0d0) + (35.0d0)*Cos(theta)**(4.0d0)))/(16.*Sqrt(Pi))
end

function Y60(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(Sqrt((13.0d0)*Pi**(-1.0d0))*((-5.0d0) + (105.0d0)*Cos(theta)**(2.0d0) + (-315.0d0)*Cos(theta)**(4.0d0) + (231.0d0)*Cos(theta)**(6.0d0)))/32.
end

function Y80(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(Sqrt((17.0d0)*Pi**(-1.0d0))*((35.0d0) + (-1260.0d0)*Cos(theta)**(2.0d0) + (6930.0d0)*Cos(theta)**(4.0d0) + (-12012.0d0)*Cos(theta)**(6.0d0) + (6435.0d0)*Cos(theta)**(8.0d0)))/256.
end

function Y100(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(Sqrt((21.0d0)*Pi**(-1.0d0))*((-63.0d0) + (46189.0d0)*Cos(theta)**(10.0d0) + (3465.0d0)*Cos(theta)**(2.0d0) + (-30030.0d0)*Cos(theta)**(4.0d0) + (90090.0d0)*Cos(theta)**(6.0d0) + (-109395.0d0)*Cos(theta)**(8.0d0)))/512.
end

function Y120(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(5*((231.0d0) + (-1939938.0d0)*Cos(theta)**(10.0d0) + (676039.0d0)*Cos(theta)**(12.0d0) + (-18018.0d0)*Cos(theta)**(2.0d0) + (225225.0d0)*Cos(theta)**(4.0d0) + (-1021020.0d0)*Cos(theta)**(6.0d0) + (2078505.0d0)*Cos(theta)**(8.0d0)))/(2048.*Sqrt(Pi))
end

function Y140(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(Sqrt((29.0d0)*Pi**(-1.0d0))*((-429.0d0) + (22309287.0d0)*Cos(theta)**(10.0d0) + (-16900975.0d0)*Cos(theta)**(12.0d0) + (5014575.0d0)*Cos(theta)**(14.0d0) + (45045.0d0)*Cos(theta)**(2.0d0) + (-765765.0d0)*Cos(theta)**(4.0d0) + (4849845.0d0)*Cos(theta)**(6.0d0) + (-14549535.0d0)*Cos(theta)**(8.0d0)))/4096.
end

function Y160(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(Sqrt((33.0d0)*Pi**(-1.0d0))*((6435.0d0) + (-1487285800.0d0)*Cos(theta)**(10.0d0) + (1825305300.0d0)*Cos(theta)**(12.0d0) + (-1163381400.0d0)*Cos(theta)**(14.0d0) + (300540195.0d0)*Cos(theta)**(16.0d0) + (-875160.0d0)*Cos(theta)**(2.0d0) + (19399380.0d0)*Cos(theta)**(4.0d0) + (-162954792.0d0)*Cos(theta)**(6.0d0) + (669278610.0d0)*Cos(theta)**(8.0d0)))/65536.
end

function Y180(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(Sqrt((37.0d0)*Pi**(-1.0d0))*((-12155.0d0) + (10039179150.0d0)*Cos(theta)**(10.0d0) + (-17644617900.0d0)*Cos(theta)**(12.0d0) + (18032411700.0d0)*Cos(theta)**(14.0d0) + (-9917826435.0d0)*Cos(theta)**(16.0d0) + (2268783825.0d0)*Cos(theta)**(18.0d0) + (2078505.0d0)*Cos(theta)**(2.0d0) + (-58198140.0d0)*Cos(theta)**(4.0d0) + (624660036.0d0)*Cos(theta)**(6.0d0) + (-3346393050.0d0)*Cos(theta)**(8.0d0)))/131072.
end

function Y200(theta) result (Y)
implicit none
real(kind=dp), intent(in) :: theta ! latitude
real(kind=dp)             :: Y
Y=(Sqrt((41.0d0)*Pi**(-1.0d0))*((46189.0d0) + (-116454478140.0d0)*Cos(theta)**(10.0d0) + (273491577450.0d0)*Cos(theta)**(12.0d0) + (-396713057400.0d0)*Cos(theta)**(14.0d0) + (347123925225.0d0)*Cos(theta)**(16.0d0) + (-167890003050.0d0)*Cos(theta)**(18.0d0) + (34461632205.0d0)*Cos(theta)**(20.0d0) + (-9699690.0d0)*Cos(theta)**(2.0d0) + (334639305.0d0)*Cos(theta)**(4.0d0) + (-4461857400.0d0)*Cos(theta)**(6.0d0) + (30117537450.0d0)*Cos(theta)**(8.0d0)))/524288.
end

