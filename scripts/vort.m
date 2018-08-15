function w=vort(u1,u2)
global kx ky

fu1=fft2(u1); % (me) Fast fourier transform on u1-velocity component
fu2=fft2(u2); % (me) Fast fourier transform on u2-velocity component
% (me) Formula for Fourier transform of vorticity?
fw=complex(0,kx).*fu2-complex(0,ky).*fu1; 
% (me) Inverse of fast fourier transform to derive vorticity field.
% (me) Question: What don't we immediately use w = u2-u1 ???
w=ifft2(fw,'symmetric'); 