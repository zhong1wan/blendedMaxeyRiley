function [u1, u2]=u0(init_cond,inputfile)
global n1 n2 kx ky x y nu

if init_cond==1 %read from file
%     wfile=['turb_u_' sprintf('%4.4i',inputfile)];
%     load(wfile);
    load(inputfile);
elseif init_cond==2 %random with gaussian decay
    k=sqrt(kx.^2+ky.^2); index=find(k==0); k(index)=1;
    k0=2;
    phi=rand(n2,n1);
    E=k.^4.*exp(-(k./k0).^2);
    E0=2e1*sqrt(2)/sum(sum(E));
    fu1=E0*(-ky./k).*sqrt(E./(2*pi*k)).*exp(complex(0,2*pi*phi));
    fu2=E0*( kx./k).*sqrt(E./(2*pi*k)).*exp(complex(0,2*pi*phi));
    fu1=fu1.*(n1*n2); fu2=fu2.*(n1*n2);
    [fu1,fu2]=dealiase(fu1,fu2);
    fu1(index)=0; fu2(index)=0;
    u1 = ifft2(fu1,'symmetric');
    u2 = ifft2(fu2,'symmetric');
    wfile=['C:\Users\15693\Desktop\turb_u_' sprintf('%4.4i',0)]; 
%    save(wfile, 'u1','u2'); 
elseif init_cond==3 % random phase - exp decay
    k=sqrt(kx.^2+ky.^2); index=find(k==0); k(index)=1;
    k0=1;
    phi=rand(n2,n1);
    E=exp(-k./k0);
    E0=-sqrt(2);
    fu1=E0*(-ky./k).*sqrt(E./(2*pi*k)).*exp(complex(0,2*pi*phi));
    fu2=E0*( kx./k).*sqrt(E./(2*pi*k)).*exp(complex(0,2*pi*phi));
    
    % set I=I_0
    I0=.1;
    fu1(5,1)=-complex(0,I0);
    fu1(125,1) = conj(fu1(5,1));
    
    fu1=fu1*(n1*n2); fu2=fu2*(n1*n2);
    [fu1,fu2]=dealiase(fu1,fu2);
    fu1(index)=0; fu2(index)=0;
    % Fourier to Phys.
    u1 = ifft2(fu1,'symmetric');
    u2 = ifft2(fu2,'symmetric');
    wfile=['turb_u_' sprintf('%4.4i',0)];
    save(wfile, 'u1','u2');
elseif init_cond==4 % exponential decay - uniformly random phase - gaussian distributed amp
    k=sqrt(kx.^2+ky.^2); index=find(k==0); k(index)=1;
    k0=2;
    rng('shuffle');
    phi=rand(n2,n1);
    Phi=normrnd(0*k,exp(-k./k0));
    E0=sqrt(2);
    fu1=E0*(-ky./k).*sqrt(Phi.^2./(2*pi*k)).*exp(complex(0,2*pi*phi));
    fu2=E0*( kx./k).*sqrt(Phi.^2./(2*pi*k)).*exp(complex(0,2*pi*phi));
    
    fu1(index)=0; fu2(index)=0;
    
    % set I=I_0
    I0=.1;
    fu1(5,1)=-complex(0,I0);
    fu1(125,1) = conj(fu1(5,1));
    
    fu1=fu1*(n1*n2); fu2=fu2*(n1*n2);
    [fu1,fu2]=dealiase(fu1,fu2);
    
    % Fourier to Phys.
    u1 = ifft2(fu1,'symmetric');
    u2 = ifft2(fu2,'symmetric');
    wfile=['turb_u_' sprintf('%4.4i',0)];
    save(wfile, 'u1','u2');
elseif init_cond==5
    u1 = (1/(16*nu))*sin(4*y);
    u2 = 0*u1;
    wfile=['C:\Users\15693\Desktop\turb_u_' sprintf('%4.4i',0)];
    save(wfile, 'u1','u2');
end

function [u1,u2]=gvortex(xx,yy,x0,y0,delta,gamma0)
r2=(xx-x0).^2+(yy-y0).^2;
index=find(r2==0);
u1 = gamma0*(1-exp(-r2/delta^2))./(2*pi*r2).*(y0-yy);
u1(index)=gamma0/(2*pi);
u2 = gamma0*(1-exp(-r2/delta^2))./(2*pi*r2).*(xx-x0);
u2(index)=gamma0/(2*pi);
