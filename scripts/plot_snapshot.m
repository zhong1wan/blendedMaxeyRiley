clear all
close all
global kx ky

load(['../velocity_results/turb_u_' sprintf('%4.4i',100)]);
load('../Re250_force1_k4/snap2500.mat');

kx1 = [0:n1/2-1 -n1/2:-1];
ky1 = [0:n2/2-1 -n2/2:-1];

[kx,ky] = meshgrid(kx1,ky1);

[x,y,kx,ky]=gvars(n1,n2,x_right,y_right);
plot_id = [1, 500, 2500];
plot_time = [0, 1, 5];

figure(1);
set(gcf, 'Position', [300, 300, 1200, 280]);

for i = 1:length(plot_id)
    
    w=vort(squeeze(U1(plot_id(i),:,:)),squeeze(U2(plot_id(i),:,:)));
    subplot(1, length(plot_id), i);
    hf1 = pcolor(x, y, w);
    set(hf1,'EdgeColor','none')

    xlim([x_left,x_right]);
    ylim([x_left,x_right]);
    xlabel('$x$','interpreter','latex','fontsize',16,'rot',0)
    ylabel('$y$','interpreter','latex','fontsize',16,'rot',0)
    title(['$t=',num2str(plot_time(i)),'$'], 'interpreter', 'latex', 'fontsize', 16);
    colorbar
    caxis([-60, 60])
    c1 = colorbar;
    ylabel(c1,'$\omega_z$','interpreter','latex','fontsize',16,'rot',0)
end


