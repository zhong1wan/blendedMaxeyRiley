clear all
close all

global kx ky

load('C:\Users\15693\Desktop\NS_384x384_Re_200_k4\traj_3957')
% flag1 = xxp_history;
% flag2 = AA;
% 
kx1 = [0:n1/2-1 -n1/2:-1];
ky1 = [0:n2/2-1 -n2/2:-1];

[kx,ky] = meshgrid(kx1,ky1);
% 
% load('C:\Users\15693\Desktop\NS_384x384_Re_200_k4\traj_2359')
% flag3 = xxp_history;
% flag4 = AA;
% load('C:\Users\15693\Desktop\NS_384x384_Re_200_k4\traj_2198')
% xxp_history_new(:,:,1:200) = xxp_history;
% xxp_history_new(:,:,201:400) = flag1;
% xxp_history_new(:,:,401:600) = flag3;
% AA_new(1:200,:) = AA;
% AA_new(201:400,:) = flag2;
% AA_new(401:600,:) = flag4;
% np = 600;
% clear xxp_history AA
% xxp_history = xxp_history_new;
% AA = AA_new;

nx = n1;
nz = n2;
xx = linspace(x_left,x_right,n1);
yy = linspace(y_left,y_right,n2);
dx = xx(2)-xx(1);
aa = AA(1,1);

ff = linspace(0,2*pi);

% eta2D = eta;
% ff = linspace(0,2*pi);
% xx_flag = zeros(nx,nz);
% for ii=1:nz
%     xx_flag(:,ii) = linspace(x_left,x_right,nx);
% end

eta_video = VideoWriter('C:\Users\15693\Desktop\prop_history.mp4','MPEG-4');
eta_video.FrameRate = 27;
open(eta_video);
% zz(1:nx) = -1.2;
% XX = [xx,fliplr(xx)];
% YY = [-yy, zz];
% x_lim = [0., 95]; y_lim = [-0.1, 0.5];
fontlat = 12;
figeta = figure(2); figeta.Name = 'Particle trajectories using Maxey-Riley t = '; figeta.Color = 'w';
scrsz = get(groot,'ScreenSize');
figeta.Position = [1 scrsz(4)/15 scrsz(3)/1.2 scrsz(4)/1.2];
%figeta.Position = [ 1719          85         1054         591];
timestart = 1;
for ii=timestart:1:timestart+1108-1
   title_string = ['History of Free-Surface elevation at time t = ' num2str(dt*double(ii),3+floor(log(dt*double(ii)))) , ' [sec]'];
   figeta.Name = title_string;
   string_flag = ['C:\Users\15693\Desktop\512x512_Re200_k12\turb_u_' sprintf('%4.4i',ii)];
   load(string_flag);
   vel = u1.^2+u2.^2;
   w=vort(u1,u2);
%     if (ii == timestart)
%        hf1 = pcolor(xx,yy,vel);
%        set(hf1,'EdgeColor','none')
%     end
%     set(hf1,'CData',vel(:,:)');
    hf1 = pcolor(xx,yy,w);
    set(hf1,'EdgeColor','none')
    colorbar
    %caxis([0.,100.])
    c1 = colorbar;
    %colormap(flipud(bone));
    ylabel(c1,'$\omega_z$','interpreter','latex','fontsize',14,'rot',0)
    %ylabh = get(c1,'YLabel');
    %set(ylabh,'Position',[2*pi-3.5 10. 0])
    hold on
   %set(subplot(2,1,1),'Position', [0.05, 0.46, 0.92, 0.50])
   %ii1 = mod(ii,99);
   %ii1 = max(ii1,1);
%    ii1 = ii;
%    amp_init = max(eta2D(:,1));
   
%    for jj=1:nx
%        depth(jj,:) = linspace(-hh(jj),eta2D(jj,ii1),nz); 
%    end
   
   subplot(1,1,1)
   
   %contourf(xx_flag,depth,veloc(:,:,ii1),'LineStyle','none','LevelList',[0.:0.02:0.8])
   %hold on
   %plot(xx,eta2D(:,ii1),'Color',[0.1,0.1,0.1],'linewidth',1.5);
   %hold on
   
   
   for pp = 1:np
       line_flag = max(1,ii-1);
       if (abs(xxp_history(1,line_flag,pp)-xxp_history(1,ii,pp)) > xx(nx)-xx(100))
          line_flag = max(1,ii); 
       end
       if (abs(xxp_history(2,line_flag,pp)-xxp_history(2,ii,pp)) > xx(nx)-xx(100))
          line_flag = max(1,ii); 
       end
       line_flag = max(1,ii); 
   if (pp < 201)
   plot(xxp_history(1,line_flag:ii,pp),xxp_history(2,line_flag:ii,pp),'Color',[1.0,0,0],'linewidth',1.5)
   hold on
   aa_size = 1.*dx/aa;
   plot(xxp_history(1,ii,pp)+aa_size*AA(pp,ii)*cos(ff),xxp_history(2,ii,pp)+aa_size*AA(pp,ii)*sin(ff),'Color',[1.0,0,0],'linewidth',1.5)
   elseif (pp < 401)
   plot(xxp_history(1,line_flag:ii,pp),xxp_history(2,line_flag:ii,pp),'Color',[1.0,0,0],'linewidth',1.5)
   hold on
   aa_size = 1.*dx/aa;
   plot(xxp_history(1,ii,pp)+aa_size*AA(pp,ii)*cos(ff),xxp_history(2,ii,pp)+aa_size*AA(pp,ii)*sin(ff),'Color',[1.0,0,0],'linewidth',1.5)  
   elseif (pp < 601)
   plot(xxp_history(1,line_flag:ii,pp),xxp_history(2,line_flag:ii,pp),'Color',[1.0,0,0],'linewidth',1.5)
   hold on
   aa_size = 1.*dx/aa;
   plot(xxp_history(1,ii,pp)+aa_size*AA(pp,ii)*cos(ff),xxp_history(2,ii,pp)+aa_size*AA(pp,ii)*sin(ff),'Color',[1.0,0,0],'linewidth',1.5) 
   end
   
   hold on
   end
   
   %plot(xx,-hh,'Color',[0.1,0.1,0.1],'linewidth',1.5)
   %hold on
   %plot([x_left,x_right],[amp_init,amp_init],'k--')
   set(gca,'Color',[0,0,0])
   hold off
   grid on
   %hold on
   %plot([0., 100.],[-0.1, 1.5],'k--')
   %hold off
   
   max_aa = (max(AA(:,ii)))/aa;
   min_aa = (min(AA(:,ii)))/aa;
   
   %str = ['Particle trajectories using Maxey-Riley, at time: ',num2str(dt*(double(ii)-1.),'%4.2f'),' $(sec)$ with $\rho_p = 1025 [kg/m^3]$'];
   str = ['Maxey-Riley simulation with $\rho_p = 1 [kg/m^3]$ with radius  $R_0 = 1 \cdot 10^{-3}$'];
   %str = ['MR & RP simulation, at time: '];
   title(str,'interpreter','latex','color','k','fontsize',14)
   %ylabel('$\frac{\eta}{h_{0}}$','interpreter','latex','fontsize',20,'rot',0)
   ylabel('$y$','interpreter','latex','fontsize',16,'rot',0)
   %xlabel('$\frac{x}{h_{0}}$','interpreter','latex','fontsize',20,'rot',0)
   xlabel('$x$','interpreter','latex','fontsize',16,'rot',0)
   %set(gca,'YTick',[0,1,2,3,4,5,6] );
   xlim([x_left,x_right])
   ylim([y_left,y_right])
   
   frame = getframe(figeta);
   writeVideo(eta_video,frame);
end

close(eta_video);