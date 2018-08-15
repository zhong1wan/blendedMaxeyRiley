clc;
clear;
global n1 n2 L1 L2 kx ky x y nu dt phase1 phase2 kappa_min kappa_max
L1=2*pi; n1=128; % Length of domain in the x1-direction and number of points.
L2=2*pi; n2=128; % Length of domain in the x2-direction and number of points. 

[x,y,kx,ky]=gvars(n1,n2,L1,L2);

phase1 = rand(n1,n2);
phase2 = rand(n1,n2);

tinit =0;  % Initial time instant
tend  =10.0; % Final time instant
dwrite =2e-3; % time-step
dt = dwrite;
ic=2; % Initial condition indicator (=2 for random with gaussian decay)
Re=250; nu=1/Re; % Reynolds number

% Forcing will be applied at wavenumbers kappa_min < sqrt{kx^2+ky^2} <
% kappa_max
kappa_max = 5;
kappa_min = 4;

times = (tend-tinit)/dwrite;

% [u1, u2]=u0(ic,0);
[u1,u2] = u0(1,'../turb_u_0100.mat');
[u1,u2] = dns2d(x,y,u1,u2,tinit,tend,dwrite,times);

%plot_turb(tinit,tend,dwrite)