clear all
close all
clc

%load('C:\Users\15693\Desktop\velocities.mat');
load(['C:\Users\15693\Desktop\New_kn_forcing\turb_u_' sprintf('%4.4i',1)]);
%u1 = zeros(n1,n2);
%u2 = zeros(n1,n2);

eta_video = VideoWriter('C:\Users\15693\Desktop\eta_history.mp4','MPEG-4');
eta_video.FrameRate = 24;
open(eta_video);
figeta = figure(1);
timestart = 100;
for tt=timestart:150
    load(['C:\Users\15693\Desktop\New_kn_forcing\turb_u_' sprintf('%4.4i',tt)]);
    %u1(:,:) = uu1(tt,:,:);
    %u2(:,:) = uu2(tt,:,:);
    vel = u1.^2+u2.^2;
    if (tt == timestart)
       hf1 = pcolor(x,y,vel);
       set(hf1,'EdgeColor','none')
    end
    set(hf1,'CData',vel(:,:)');
    frame = getframe(figeta);
    writeVideo(eta_video,frame);
end

close(eta_video);



