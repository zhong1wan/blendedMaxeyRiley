function [x,y,kx,ky]=gvars(n1,n2,L1,L2)

X=linspace(0,L1-L1/n1,n1);
Y=linspace(0,L2-L2/n2,n2);

[x,y]=meshgrid(X,Y);

%%% old way
% dx=L1/n1; dy=L2/n2;
% kx0 = (mod(1/2 + (0:(n1-1))/n1 , 1) - 1/2); kx1 = kx0 * (2*pi/dx);
% ky0 = (mod(1/2 + (0:(n2-1))/n2 , 1) - 1/2); ky1 = ky0 * (2*pi/dy);

kx1 = [0:n1/2-1 -n1/2:-1]*(2*pi/L1);
ky1 = [0:n2/2-1 -n2/2:-1]*(2*pi/L2);

[kx,ky] = meshgrid(kx1,ky1);