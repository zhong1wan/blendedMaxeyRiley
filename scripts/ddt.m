function [du1dt,du2dt] = ddt(u1,u2,t)
global kx ky

[du1dt,du2dt] = rhs(u1,u2,t);

%%% calculate -u dot grad(u) %%%
% calculate u cross w
w=vort(u1,u2);
adv1 =  u2.*w; adv2 = -u1.*w;
fadv1 = fft2(adv1); fadv2=fft2(adv2);

% calculate -grad(.5*|u|^2)
usq = -.5*(u1.^2 + u2.^2);
fusq = fft2(usq);
fadv1 = fadv1 + 1i*kx.*fusq;
fadv2 = fadv2 + 1i*ky.*fusq;

% dealiasing
[fadv1,fadv2] = dealiase(fadv1,fadv2);

adv1 = ifft2(fadv1,'symmetric');
adv2 = ifft2(fadv2,'symmetric');

%%% DDt = ddt + adv %%% 
du1dt = du1dt - adv1;
du2dt = du2dt - adv2;
