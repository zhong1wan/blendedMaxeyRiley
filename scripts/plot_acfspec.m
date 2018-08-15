%%% 
% read velocity from different files and plot pdf
% (contain initial conditions starting from different times)
%
%%%

clear all;

% label = 10*(0:1);
% label = [0, 10, 20, 30, 40, 50, 60, 69];
label = [0];
numLags = 500; dt = 4e-3; fs = 1/dt;
Vt = []; Vp = []; Vm = [];
acf_t = zeros(numLags+1,1); acf_p = acf_t; acf_m = acf_t;
vdp = 0;

for i = 1:length(label)

    rfile = ['../particle_results_organized/test2/velocity_results_rk4/v2000_2ksteps_t' sprintf('%2.2i',label(i)) '.mat'];
    load(rfile);
    Mt = [squeeze(vt(:,:,1)), squeeze(vt(:,:,2))];
    Mp = [squeeze(vp(1:2:end,:,1)), squeeze(vp(1:2:end,:,2))];
    Mm = [squeeze(vm(:,:,1)), squeeze(vm(:,:,2))];
% %     mu = mean(mean(Mt)); sigma2 = mean(mean((Mt-mu).^2));
% %     acf_t = ac(Mt-mu, 250);
    acf_t = (i-1)/i*acf_t + ac(Mt,numLags)/i;
    acf_p = (i-1)/i*acf_p + ac(Mp,numLags)/i;
    acf_m = (i-1)/i*acf_m + ac(Mm,numLags)/i;
    
    disp(i)
end

% [tau,CCt] = reflectx((0:numLags)*dt, acf_t/acf_t(1));
% [tau,CCp] = reflectx((0:numLags)*dt, acf_p/acf_p(1));
% [tau,CCm] = reflectx((0:numLags)*dt, acf_m/acf_m(1));
figure; hold on;
% plot(tau, CCt);
% plot(tau, CCp);
% plot(tau, CCm);
plot((0:numLags)*dt, acf_t);
plot((0:numLags)*dt, acf_p);
plot((0:numLags)*dt, acf_m);
legend({'true','RNN manifold','order 1 slow manifold'},'interpreter','latex');
legend boxoff;
set(gca, 'LineWidth', 1.5, 'FontSize', 14, 'TickLabelInterpreter','latex');

%%%%%% spectrum %%%%%%
Sw_t = real(fft(acf_t));
Sw_p = real(fft(acf_p));
Sw_m = real(fft(acf_m));
f = fs/501*(0:500);
% [ff,Swt] = reflectx(f(1:251), Sw_t(1:251));
% [ff,Swp] = reflectx(f(1:251), Sw_p(1:251));
% [ff,Swm] = reflectx(f(1:251), Sw_m(1:251));

figure; 
hold on; box off;
plot(f(1:251), Sw_t(1:251), 'k-', 'linewidth', 1);
plot(f(1:251), Sw_p(1:251), 'linewidth', 1);
plot(f(1:251), Sw_m(1:251), 'linewidth', 1);
% plot(ff, Swt, 'k-', 'linewidth', 1);
% plot(ff, Swp, 'linewidth', 1);
% plot(ff, Swm, 'linewidth', 1);
xlim([0, 20]);
xlabel('$\omega$ (Hz)','interpreter','latex');
ylabel('$S(\omega)$','interpreter','latex');
legend({'true','RNN manifold','order 1 slow manifold'},'interpreter','latex');
legend boxoff;
set(gca, 'LineWidth', 1.5, 'FontSize', 14, 'TickLabelInterpreter','latex');


function acf = ac(V, lags)
% assume V has dimension [time, realizations]
    acf = zeros(lags+1,1);
    m = mean(mean(V));
    for tau = 0:lags
        A = (V(1:end-tau,:)-m).*(V(1+tau:end,:)-m);
        acf(tau+1) = mean(mean(A));
    end
end

function [x2, y2] = reflectx(x, y)
    x2 = [-fliplr(x(2:end)), x];
    y2 = [flipud(y(2:end)); y];
end
