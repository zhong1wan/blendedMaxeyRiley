function [rhs1, rhs2]=rhs(u1,u2,t)
global kx ky nu n1 n2

fu1=fft2(u1); fu2=fft2(u2);

ksq = kx.^2+ky.^2;

% calculate u cross w
w=vort(u1,u2);
rhs1 =  u2.*w; rhs2 = -u1.*w;
frhs1 = fft2(rhs1); frhs2=fft2(rhs2);

% calculate -grad(.5*|u|^2)
usq = -.5*(u1.^2 + u2.^2);
fusq = fft2(usq);
frhs1 = frhs1 + 1i*kx.*fusq;
frhs2 = frhs2 + 1i*ky.*fusq;

% add force + drag
% (me) this function I need to change to 
% (me) include random gaussian forcing and not Kolmogorov forcing.
% (me) The forcing should be transformed to Fourier space.
[frc1, frc2]=force_new(t);
frhs1=frhs1+frc1;
frhs2=frhs2+frc2;

% dissipation
% (me) add dissipation term (1/Re)*grad^2(f) => (1/Re)*k^2 *f.
frhs1 = frhs1- nu*ksq.*fu1; 
frhs2 = frhs2- nu*ksq.*fu2;

% Attempt to add hypoviscosity
for ii=1:n1
   for jj = 1:n2
      if (ksq(ii,jj) > 0.0)
         ksq_inv(ii,jj) = 1/ksq(ii,jj); 
      end
   end
end
frhs1 = frhs1+ 0.0*ksq_inv.*fu1; 
frhs2 = frhs2+ 0.0*ksq_inv.*fu2;

% dealiase + project to div free fields
% (me) Question: Why isn't the aliasing done after projection
% (me) to the divergence-free subspace?
[frhs1, frhs2] = dealiase(frhs1,frhs2);
% (me) project the function f to a divfree subspace
% (me) u = P(f), via grad(p).
[frhs1, frhs2] = divfree(frhs1,frhs2);

rhs1 = ifft2(frhs1,'symmetric');
rhs2 = ifft2(frhs2,'symmetric');