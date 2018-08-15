 function [f1,f2] = dealiase(f1,f2)
% 2/3 dealiasing
global kx ky n1 n2

RorC = isreal(f1);
if RorC==1
    f1=fft2(f1);
    f2=fft2(f2);
end

% k=sqrt(kx.^2+ky.^2);
% k0=max(max(k))/3;
% k_cutoff=(k<k0);

% Chandler & Kerswell use this
ksq=kx.^2+ky.^2;
k0=(n1/3)^2;
k_cutoff=(ksq<=k0);

% k_cutoff=(abs(kx)<n1/3 & abs(ky)<n2/3);

f1=k_cutoff.*f1; f2=k_cutoff.*f2;

if RorC==1
    f1=ifft2(f1,'symmetric');
    f2=ifft2(f2,'symmetric');
end