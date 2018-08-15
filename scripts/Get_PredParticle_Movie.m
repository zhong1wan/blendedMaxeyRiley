clear all
close all
clc
global kx ky

% load file that contains parameters (grid, wavenumbers, etc.)
load(['../Re250_force1_k4/turb_u_' sprintf('%4.4i',100)]);

kx1 = [0:n1/2-1 -n1/2:-1];
ky1 = [0:n2/2-1 -n2/2:-1];
[kx,ky] = meshgrid(kx1,ky1);
[x,y,kx,ky]=gvars(n1,n2,x_right,y_right);

% load flow snapshots
load('../Re250_force1_k4/snap2500.mat');
dt_flow = 2e-3;

% load particle positions
% load('../particle_results_organized/test2/particlex_1000by1200steps_ep005.mat');
% load('../particle_results_organized/test2/multistep_pred_1000pts.mat');
load('../particle_results_organized/test2/multistep_pred_2kpts_rk4.mat');
dt_traj = 4e-3;
xt = xt(201:end,:,:);
xp = xp(1:2:end,:,:);

% video writer parameters
eta_video = VideoWriter('../trajectory_data/particle1000_truth.mp4','MPEG-4');

dt_plot = .04;
p = 1000;
ndt_flow = dt_plot/dt_flow;    % this needs to be an integer
ndt_traj = dt_plot/dt_traj;
idx1_flow = dt_traj/dt_flow + 401;    % index of first flow snapshots to plot
eta_video.FrameRate = 25;

open(eta_video);
figeta = figure(1);
timestart = 1;
times = size(xp,1)/ndt_traj - 1;
times = 50;

for tt=timestart:1:times
    
    % plot flow field and particles
    idx_flow = idx1_flow + (tt-1)*ndt_flow;
    idx_part = 1 + (tt-1)*ndt_traj;
    w = vort(squeeze(U1(idx_flow,:,:)), squeeze(U2(idx_flow,:,:)));
    if (tt == timestart)
        subplot(1,2,1);
        hf1 = pcolor(x,y,w(:,:));
        hold on;
        set(hf1,'EdgeColor','none')
        sc1 = scatter(xt(idx_part,1:p,1), xt(idx_part,1:p,2), 20, 'r', 'filled');
        title('truth','interpreter','latex','fontsize',16);
        
        subplot(1,2,2);
        hf2 = pcolor(x,y,w(:,:));
        hold on;
        set(hf2,'EdgeColor','none')
        sc2 = scatter(xp(idx_part,1:p,1), xp(idx_part,1:p,2), 20, 'r', 'filled');
        title('prediction','interpreter','latex','fontsize',16);
        
        set(gcf, 'Position', [100, 100, 1000, 400]);
    end
    
    set(gcf,'color','w');
    
    set(hf1,'CData',w(:,:));
    set(sc1,'XData',xt(idx_part,1:p,1),'YData',xt(idx_part,1:p,2));
    
    set(hf2,'CData',w(:,:));
    set(sc2,'XData',xp(idx_part,1:p,1),'YData',xp(idx_part,1:p,2));
    
    subplot(1,2,1);
    xlim([x_left,x_right]);
    ylim([x_left,x_right]);
    xlabel('$x$','interpreter','latex','fontsize',16,'rot',0)
    ylabel('$y$','interpreter','latex','fontsize',16,'rot',0)
    colorbar
    caxis([-60, 60])
    c1 = colorbar;
    ylabel(c1,'$\omega_z$','interpreter','latex','fontsize',16,'rot',0)
    
    subplot(1,2,2);
    xlim([x_left,x_right]);
    ylim([x_left,x_right]);
    xlabel('$x$','interpreter','latex','fontsize',16,'rot',0)
    ylabel('$y$','interpreter','latex','fontsize',16,'rot',0)
    colorbar
    caxis([-60, 60])
    c2 = colorbar;
    ylabel(c2,'$\omega_z$','interpreter','latex','fontsize',16,'rot',0)
    
    frame = getframe(figeta);
    writeVideo(eta_video,frame);
    disp(tt);
end

close(eta_video);

