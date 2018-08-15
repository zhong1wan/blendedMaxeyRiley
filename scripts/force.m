function [frc1,frc2]=force(famp)
global y

%(me)
flag = randn;
% (me)
a1 = (1.+0.0*flag)*famp*sin(2*y);
% Farazmand forcing
%a1 = famp*sin(4*y);
a2 = 0*a1;
frc1 = fft2(a1);
frc2 = fft2(a2);