function []=plot_turb(t0,tf,dwrite)
global x y
nt = round((tf-t0)/dwrite)+1;

figure;
for j=0:nt-1
    j
    rfile = ['C:\Users\15693\Desktop\turb_u_' sprintf('%4.4i',j) '.mat'];
    load(rfile);
    w=vort(u1,u2);
    
    contourf(x,y,w,linspace(min(w(:)),max(w(:))-.01,15));
    axis equal tight; colorbar
    pause(1e-5);
end
