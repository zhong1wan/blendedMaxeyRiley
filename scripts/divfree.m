function [f1 f2]=divfree(f1,f2);
global kx ky

% (me) f = u + grad(p)
tmp1=f1; tmp2=f2;

ksq = kx.^2+ky.^2;
index = find(ksq==0);
ksq(index)=1;

% (me) In Fourier space we project the f-solution into 
% (me) a subspace where div(u) = 0 (incompressibility condition
% (me) via the grad(p). We then get u from f.
f1 = (ky.^2.*tmp1 - kx.*ky.*tmp2)./ksq;
f2 = (-kx.*ky.*tmp1 + kx.^2.*tmp2)./ksq;