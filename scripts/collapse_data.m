clc; clear all; close all;

snapshots = 4999;
U1 = zeros(snapshots,128,128); U2 = zeros(snapshots,128,128);
DU1Dt = zeros(snapshots,128,128); DU2Dt = zeros(snapshots,128,128);

for i = 1:snapshots
    load(['../Re250_force1_k4/turb_u_' sprintf('%4.4i',i)]);
    U1(i,:,:) = u1; U2(i,:,:) = u2;
    DU1Dt(i,:,:) = Du1Dt; DU2Dt(i,:,:) = Du2Dt;
end

save('../Re250_force1_k4/snap4999.mat','U1','U2','DU1Dt','DU2Dt');
