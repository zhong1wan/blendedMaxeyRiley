clear all
close all
global kx ky

%load(['C:\Users\15693\Desktop\Other_Cases\Navier_Stokes\velocity_results\turb_u_' sprintf('%4.4i',1)]);
%load(['C:\Users\15693\Desktop\512x512_Re200_k8_Hypo\turb_u_' sprintf('%4.4i',1)]);
%load(['C:\Users\15693\Desktop\New_kn_forcing\turb_u_' sprintf('%4.4i',1)]);
folder = 'Re250_force1_k4';
load(['../' folder '/turb_u_' sprintf('%4.4i',100)]);

kx1 = [0:n1/2-1 -n1/2:-1];
ky1 = [0:n2/2-1 -n2/2:-1];

[kx,ky] = meshgrid(kx1,ky1);

[x,y,kx,ky]=gvars(n1,n2,x_right,y_right);

% eta_video = VideoWriter('C:\Users\15693\Desktop\eta_history.mp4','MPEG-4');
eta_video = VideoWriter(['../' folder '/eta_history.mp4'],'MPEG-4');
eta_video.FrameRate = 10;
open(eta_video);
figeta = figure(1);
timestart = 1;
times = 100;
for tt=timestart:1:times
    load(['../' folder '/turb_u_' sprintf('%4.4i',tt)]);
    %u1(:,:) = uu1(tt,:,:);
    %u2(:,:) = uu2(tt,:,:);
    %vel = u1.^2+u2.^2;
    w=vort(u1,u2);
    if (tt == timestart)
       hf1 = pcolor(x,y,w(:,:));
       set(hf1,'EdgeColor','none')
    end
    set(hf1,'CData',w(:,:));
    set(gcf,'color','w');
    xlim([x_left,x_right]);
    ylim([x_left,x_right]);
    xlabel('$x$','interpreter','latex','fontsize',16,'rot',0)
    ylabel('$y$','interpreter','latex','fontsize',16,'rot',0)
    colorbar
    caxis([-70, 70])
    c1 = colorbar;
    ylabel(c1,'$\omega_z$','interpreter','latex','fontsize',16,'rot',0)
    %ylabh = get(c1,'YLabel');
    %set(ylabh,'Position',[3.0 0.85 0])
    frame = getframe(figeta);
    writeVideo(eta_video,frame);
    disp(tt);
end

close(eta_video);

