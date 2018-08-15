clear all
close all
global kx ky

%load(['C:\Users\15693\Desktop\Other_Cases\Navier_Stokes\velocity_results\turb_u_' sprintf('%4.4i',1)]);
load(['C:\Users\15693\Desktop\512x512_Re5000_k12\turb_u_' sprintf('%4.4i',1)]);

kx1 = [0:n1/2-1 -n1/2:-1];
ky1 = [0:n2/2-1 -n2/2:-1];

[kx,ky] = meshgrid(kx1,ky1);

eta_video = VideoWriter('C:\Users\15693\Desktop\eta_history.mp4','MPEG-4');
eta_video.FrameRate = 24;
open(eta_video);
figeta = figure(1);
scrsz = get(groot,'ScreenSize');
figeta.Position = [1 scrsz(4)/15 scrsz(3)/1.5 scrsz(4)/1.5];
timestart = 1;
for tt=timestart:3:801
    %load(['C:\Users\15693\Desktop\Other_Cases\Navier_Stokes\velocity_results\turb_u_' sprintf('%4.4i',tt)]);
    load(['C:\Users\15693\Desktop\512x512_Re5000_k12\turb_u_' sprintf('%4.4i',tt)]);
    %u1(:,:) = uu1(tt,:,:);
    %u2(:,:) = uu2(tt,:,:);
    %vel = u1.^2+u2.^2;
    subplot(1,2,1)
    w=vort(u1,u2);
    if (tt == timestart)
       hf1 = pcolor(x,y,w(:,:));
       set(hf1,'EdgeColor','none')
    end
    set(hf1,'CData',w(:,:));
    xlim([x_left,x_right]);
    ylim([x_left,x_right]);
    title('No hypoviscosity','interpreter','latex','fontsize',16,'rot',0)
    xlabel('$x$','interpreter','latex','fontsize',16,'rot',0)
    ylabel('$y$','interpreter','latex','fontsize',16,'rot',0)
    colorbar
    caxis([-65,65])
    c1 = colorbar;
    ylabel(c1,'$\omega_z$','interpreter','latex','fontsize',16,'rot',0)
    %ylabh = get(c1,'YLabel');
    %set(ylabh,'Position',[3.0 0.85 0])
    
    
    load(['C:\Users\15693\Desktop\384x384_R200_Hypo\turb_u_' sprintf('%4.4i',tt)]);
    subplot(1,2,2)
    w=vort(u1,u2);
    if (tt == timestart)
       hf2 = pcolor(x,y,w(:,:));
       set(hf2,'EdgeColor','none')
    end
    set(hf2,'CData',w(:,:));
    xlim([x_left,x_right]);
    ylim([x_left,x_right]);
    title('Hypoviscosity','interpreter','latex','fontsize',16,'rot',0)
    xlabel('$x$','interpreter','latex','fontsize',16,'rot',0)
    ylabel('$y$','interpreter','latex','fontsize',16,'rot',0)
    colorbar
    caxis([-65,65])
    c1 = colorbar;
    ylabel(c1,'$\omega_z$','interpreter','latex','fontsize',16,'rot',0)
    %ylabh = get(c1,'YLabel');
    %set(ylabh,'Position',[3.0 0.85 0])
    
    
    
    
    
    frame = getframe(figeta);
    writeVideo(eta_video,frame);
    disp(tt);
end

close(eta_video);

