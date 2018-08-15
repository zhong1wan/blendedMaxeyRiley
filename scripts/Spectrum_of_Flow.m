clear all
close all

global kx ky

%load(['C:\Users\15693\Desktop\512x512_Re200_k8_Hypo\turb_u_' sprintf('%4.4i',1)]);
%load(['velocity_results/turb_u_' sprintf('%4.4i',1)]);
load(['C:\Users\15693\Desktop\Other_Cases\Navier_Stokes\velocity_results\turb_u_' sprintf('%4.4i',1)]);
kx1 = [0:n1/2-1 -n1/2:-1];
ky1 = [0:n2/2-1 -n2/2:-1];

[x,y,kx,ky]=gvars(n1,n2,x_right,y_right);
w=vort(u1,u2);


xx = x(1,:);
yy = y(:,1);

kk_total = 31;
kx_flag = zeros(1,kk_total);
ky_flag = zeros(1,kk_total);
for ii=1:kk_total
   kx_flag(ii) = -kk_total+(kk_total-1)/2 +ii ;
   %kx_flag(ii) = 8.*cos(2.*pi*double(ii-1)/double(kk_total));
   %ky_flag(ii) = 8.*sin(2.*pi*double(ii-1)/double(kk_total));
end
%kx_flag = 10.*kx_flag/max(kx_flag);

times = 65;
Ymn = zeros(kk_total,kk_total,times);
Mmn = zeros(kk_total,kk_total);
Zmn = zeros(kk_total,kk_total);
for tt=1:1:times
%load(['C:\Users\15693\Desktop\512x512_Re200_k8_Hypo\turb_u_' sprintf('%4.4i',tt)]);
%load(['velocity_results/turb_u_' sprintf('%4.4i',tt)]);
load(['C:\Users\15693\Desktop\Other_Cases\Navier_Stokes\velocity_results\turb_u_' sprintf('%4.4i',tt)]);
w=vort(u1,u2);
    
for ii=1:kk_total
    for jj=1:kk_total
        Ymn(ii,jj,tt) = trapz(yy, trapz(xx, exp(1i*(kx_flag(ii)*x+kx_flag(jj)*y)).*w ,2));
    end
end
clc;
disp(tt);
end

%%
for ii=1:kk_total
   for jj=1:kk_total
      Mmn(ii,jj) = sum(real(Ymn(ii,jj,:)))/double(times); 
   end
end

for ii=1:kk_total
   for jj=1:kk_total
      Zmn(ii,jj) = sum( (real(Ymn(ii,jj,:)) -Mmn(ii,jj)).^2 )/double(times); 
   end
end


%save('velocity_results/Spectrum_Data','Ymn','Mmn','Zmn','kx_flag','kk_total');

%%
figure(1)
subplot(3,1,1)
h1 = pcolor(kx_flag,kx_flag,real(Ymn(:,:,11)));
set(h1,'EdgeColor','none');
subplot(3,1,2)
h2 = pcolor(kx_flag,kx_flag,real(Mmn));
set(h2,'EdgeColor','none');
subplot(3,1,3)
h3 = pcolor(kx_flag,kx_flag,real(Zmn));
set(h3,'EdgeColor','none');