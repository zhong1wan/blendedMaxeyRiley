% clc;
% close all;
clear all;

load('../particle_results_organized/test2/multistep_pred_10000pts.mat');
load('../particle_results_organized/test2/acf_m.mat');

dt = 4e-3; fs = 1/dt;
numLags = 500;
lag_max = dt*numLags;
numPoints = size(vt,1);
acfx = zeros(numPoints,numLags+1); acfy = acfx;

%%% acf for truth %%%
for i = 1:numPoints
    datax = vt(:,i,1);
    datay = vt(:,i,2); 
    acfx(i,:) = autocorr(datax, numLags);
    acfy(i,:) = autocorr(datay, numLags);
end

acfx_m = mean(acfx,1); acfy_m = mean(acfy,1);
acf_t = .5*(acfx_m + acfy_m);

%%% acf for prediction %%%
for i = 1:numPoints
    datax = vp(:,i,1);
    datay = vp(:,i,2); 
    acfx(i,:) = autocorr(datax, numLags);
    acfy(i,:) = autocorr(datay, numLags);
end

acfx_m = mean(acfx,1); acfy_m = mean(acfy,1);
acf_p = .5*(acfx_m + acfy_m);

%%% plotting %%%
figure(1);
set(gcf,'Position',[300,300,1000,400]);
subplot(1,2,1);
box off; hold on;
plot((0:numLags)*dt, acf_t, 'k-', 'linewidth', 1);
plot((0:numLags)*dt, acf_p, 'linewidth', 1);
plot((0:numLags)*dt, acf_m, 'r', 'linewidth', 1);

xlim([0, lag_max]);
xlabel('$\tau$','interpreter','latex');
ylabel('$C_{vv}(\tau)$','interpreter','latex');
set(gca, 'LineWidth', 1.5, 'FontSize', 14, 'TickLabelInterpreter','latex');

%%%%%% spectrum %%%%%%
Sw_t = real(fft(acf_t, 501));
Sw_p = real(fft(acf_p, 501));
Sw_m = real(fft(acf_m, 501));
f = fs/501*(0:500);

figure(1);
subplot(1,2,2);
hold on; box off;
plot(f(1:251), Sw_t(1:251), 'k-', 'linewidth', 1);
plot(f(1:251), Sw_p(1:251), 'linewidth', 1);
plot(f(1:251), Sw_m(1:251), 'linewidth', 1);
xlim([0, 20]);
xlabel('$\omega$ (Hz)','interpreter','latex');
ylabel('$S(\omega)$','interpreter','latex');
legend({'true','RNN manifold','order 1 slow manifold'},'interpreter','latex');
legend boxoff;
set(gca, 'LineWidth', 1.5, 'FontSize', 14, 'TickLabelInterpreter','latex');
