function [frc1,frc2]=force_new(t)
global y kx ky x n1 n2 dt phase1 phase2 kappa_min kappa_max

%kappa_f = 3;
%calculate energy of modes within the cell
%u1_f = zeros(n1,n2); % velocity corresponding to modes within cell in the x-direction
%u2_f = zeros(n1,n2); % velocity corresponding to modes within cell in the y-direction
%u1_f = u1_f +fu1(1,1)*exp(1i*0*x)/double(n1*n2);
%u2_f = u2_f +fu2(1,1)*exp(1i*0*x)/double(n1*n2);
%for ii=1:n1
%   for jj=1:n2
%       if (kx(ii,jj)^2+ky(ii,jj)^2 < kappa_f^2 +1)
%       u1_f = u1_f +real(fu1(ii,jj)*exp(1i*(kx(ii,jj)*x +ky(ii,jj)*y) )/double(n1*n2));
%       u2_f = u2_f +real(fu2(ii,jj)*exp(1i*(kx(ii,jj)*x +ky(ii,jj)*y) )/double(n1*n2));
%       end
%   end
%end
%Calculate energy as integral over space of the velocity of the modes
%within the cell
%Energy_f = trapz(y(:,1),trapz(x(1,:), 0.5*(u1_f.*u1_f+u2_f.*u2_f) ,2) );
% Calculate forcing
%P = 5;
%a1 = P*u1_f/(2*Energy_f); a2 = P*u2_f/(2*Energy_f);


% Forcing acting only on certain wavenumbers
force_amp = .6;
kappa_max = 16;
kappa_min = 15;
a1 = zeros(n1,n2);
a2 = zeros(n1,n2);

for ii=1:n1
   for jj=1:n2
       if (kx(ii,jj)^2+ky(ii,jj)^2 < kappa_max^2 +1 && kx(ii,jj)^2+ky(ii,jj)^2 > kappa_min^2 -1)
       a1 = a1 +real(force_amp*exp(1i*(kx(ii,jj)*x +ky(ii,jj)*y +2*pi*phase1(ii,jj)) ));
       a2 = a2 +real(force_amp*exp(1i*(kx(ii,jj)*x +ky(ii,jj)*y +2*pi*phase2(ii,jj)) ));
       end
   end
end
%disp(t);


% Farazmand forcing
%a1 = famp*sin(4*y);
%a2 = 0*a1;
frc1 = fft2(a1);
frc2 = fft2(a2);