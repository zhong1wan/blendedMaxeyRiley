clear;
close all

L1=2*pi; n1=192; % Length of domain in the x1-direction and number of points.
L2=2*pi; n2=192; % Length of domain in the x2-direction and number of points. 

[x,y,kx,ky]=gvars(n1,n2,L1,L2);

force_amp = 0.2;
kappa_max = 5;
kappa_min = 4;
a1 = zeros(n1,n2);
a2 = zeros(n1,n2);
phase1 = rand(n1,n2);
phase2 = rand(n1,n2);
for ii=1:n1
   for jj=1:n2
       if (kx(ii,jj)^2+ky(ii,jj)^2 < kappa_max^2 +1 && kx(ii,jj)^2+ky(ii,jj)^2 > kappa_min^2 -1)
       a1 = a1 +real(force_amp*exp(1i*(kx(ii,jj)*x +ky(ii,jj)*y +2*pi*rand(1))));
       a2 = a2 +real(force_amp*exp(1i*(kx(ii,jj)*x +ky(ii,jj)*y +2*pi*rand(1)) ));
       end
   end
end

%%
subplot(2,1,1)
hf1 = pcolor(x,y,a1);
set(hf1,'EdgeColor','none')
subplot(2,1,2)
hf2 = pcolor(x,y,a2);
set(hf2,'EdgeColor','none')
