clear all
close all

%load('C:\Users\15693\Desktop\NS_384x384_Re_200_k4\Spectrum_Data')
load('C:\Users\15693\Desktop\New_kn_forcing\Spectrum_Data')

%%

eta_video = VideoWriter('C:\Users\15693\Desktop\eta_history.mp4','MPEG-4');
eta_video.FrameRate = 24;
open(eta_video);
figeta = figure(1);
scrsz = get(groot,'ScreenSize');
figeta.Position = [1 scrsz(4)/15 scrsz(3)/1.5 scrsz(4)/1.5];
timestart = 1;
for tt=timestart:10:size(Ymn,3)-size(Ymn,3)+1

    if (tt == timestart)
       hf1 = pcolor(kx_flag,kx_flag,real(Ymn(:,:,tt)));
       set(hf1,'EdgeColor','none')
    end
    set(hf1,'CData',real(Ymn(:,:,tt)));
    xlabel('$k_x$','interpreter','latex','fontsize',16,'rot',0)
    ylabel('$k_y$','interpreter','latex','fontsize',16,'rot',0)
    colorbar
    caxis([-50,50])
    c1 = colorbar;
    ylabel(c1,'$S(k_x,k_y;t)$','interpreter','latex','fontsize',14,'rot',0)
    ylabh = get(c1,'YLabel');
    set(ylabh,'Position',[3.0 0.85 0])
    frame = getframe(figeta);
    writeVideo(eta_video,frame);
    disp(tt);
end

close(eta_video);

%%
TIME = linspace(0,80,800);
Ymnnew = Ymn(:,:,1:10:800);
TIMEnew = TIME(1:10:800);

for ii=1:kk_total
   for jj=1:kk_total 
       Mmn(ii,jj) = trapz(TIMEnew,real(Ymnnew(ii,jj,:)))/(TIMEnew(size(TIMEnew,2))-TIMEnew(1));
       Zmn(ii,jj) = trapz(TIMEnew, (real(Ymnnew(ii,jj,:)) -Mmn(ii,jj)).^2 )/(TIME(size(TIMEnew,2))-TIMEnew(1));
   end
end

%%
figeta = figure(2);
scrsz = get(groot,'ScreenSize');
figeta.Position = [1 scrsz(4)/15 scrsz(3)/1.5 scrsz(4)/1.5];
h2 = pcolor(kx_flag,kx_flag,real(Mmn));
set(h2,'EdgeColor','none')
caxis([-10,10])
xlabel('$k_x$','interpreter','latex','fontsize',16,'rot',0)
ylabel('$k_y$','interpreter','latex','fontsize',16,'rot',0)
colorbar
c1 = colorbar;
ylabel(c1,'$M(k_x,k_y)$','interpreter','latex','fontsize',14,'rot',0)
ylabh = get(c1,'YLabel');
set(ylabh,'Position',[3.0 0.85 0])

%%
figeta = figure(3);
scrsz = get(groot,'ScreenSize');
figeta.Position = [1 scrsz(4)/15 scrsz(3)/1.5 scrsz(4)/1.5];
h3 = pcolor(kx_flag,kx_flag,log10(real(Zmn)));
set(h3,'EdgeColor','none')
%caxis([0,1.7*10^7])
xlabel('$k_x$','interpreter','latex','fontsize',16,'rot',0)
ylabel('$k_y$','interpreter','latex','fontsize',16,'rot',0)
colorbar
c1 = colorbar;
ylabel(c1,'$log_{10}(Z(k_x,k_y))$','interpreter','latex','fontsize',12,'rot',0)
ylabh = get(c1,'YLabel');
set(ylabh,'Position',[3.2, 0.9 0])



% figure(1)
% subplot(3,1,1)
% h1 = pcolor(kx_flag,kx_flag,real(Ymn(:,:,10)));
% caxis([-2500,2500])
% set(h1,'EdgeColor','none');
% subplot(3,1,2)
% h2 = pcolor(kx_flag,kx_flag,real(Mmn));
% caxis([-800,600])
% set(h2,'EdgeColor','none');
% subplot(3,1,3)
% h3 = pcolor(kx_flag,kx_flag,real(Zmn));
% caxis([0,1.7*10^7])
% set(h3,'EdgeColor','none');
