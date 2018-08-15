
%%% 
% read velocity from different files and plot pdf
% (contain initial conditions starting from different times)
%
%%%

% clc;
% close all;
clear all;

label = 7*(0:6);
Vdt = []; Vdp = []; Vdm = [];

for i = 1:length(label)
    rfile = ['../particle_results_organized/test2/velocity_results/velocities1000_t' sprintf('%2.2i',label(i)) '.mat'];
    load(rfile);
    Vdt = [Vdt; vd_t(:)]; 
    Vdp = [Vdp; vd_p(:)]; 
    Vdm = [Vdm; vd_m(:)];
end

[ft, xit] = ksdensity(Vdt,'NumPoints',200);
[fp, xip] = ksdensity(Vdp,'NumPoints',200);
[fm, xim] = ksdensity(Vdm,'NumPoints',200);

%%% smoothing %%%

fts = ft;
fts(1:90) = smoothdata(fts(1:90),'gaussian',15);
fts(end-90:end) = smoothdata(fts(end-90:end),'gaussian',15);
fts(80:end-80) = smoothdata(ft(80:end-80),'gaussian',3);

fps = fp;
fps(1:90) = smoothdata(fps(1:90),'gaussian',15);
fps(end-90:end) = smoothdata(fps(end-90:end),'gaussian',15);
fps(80:end-80) = smoothdata(fp(80:end-80),'gaussian',3);

fms = fm;
fms(1:90) = smoothdata(fms(1:90),'gaussian',20);
fms(end-90:end) = smoothdata(fms(end-90:end),'gaussian',20);
fms(80:end-80) = smoothdata(fm(80:end-80),'gaussian',5);

%%% plotting %%%

figure; box off; hold on;
plot(xit, fts, 'k-', 'linewidth', 1.25);
plot(xip, fps, 'linewidth', 1.25);
plot(xim, fms, 'r-', 'linewidth', 1.25);
set(gca, 'YScale', 'log', 'LineWidth', 1.5, 'FontSize', 14, ...
    'TickLabelInterpreter','latex');
xlim([-15, 15]);
ylim([1e-4, 1e1]);
xlabel('$v-u$','interpreter','latex');
ylabel('probability density','interpreter','latex');
legend({'true','RNN manifold','order 1 slow manifold'},'interpreter','latex');
legend boxoff;
