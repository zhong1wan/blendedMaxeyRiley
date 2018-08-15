function [u1, u2]=dns2d(x,y,u1,u2,tinit,tend,dwrite,times)
global n1 n2

u1_vec = reshape(u1,n1*n2,1);
u2_vec = reshape(u2,n1*n2,1);

x_left = 0.;
x_right = 2.*pi;
y_left = 0.;
y_right = 2.*pi;
dt = dwrite;

options = odeset('RelTol',1e-5,'AbsTol',1e-5);


for t=tinit:dwrite:(tend-dwrite)
    % Use ode45 to solve the problem for t, t+dt/2 t+dt
    [~, F]=ode45(@oderhs,[t t+dwrite/2 t+dwrite],[u1_vec;u2_vec],options);
    % Take the solution for t+dt
    u1 = reshape(F(3,1:n1*n2),         n2, n1);
    u2 = reshape(F(3,n1*n2+1:2*n1*n2), n2, n1);
    % calculate total derivative
    [Du1Dt,Du2Dt] = ddt(u1,u2,t+dwrite);

    % save the vel field
    j=round(t/dwrite)+1;
    wfile=['../velocity_results/turb_u_' sprintf('%4.4i',j)];

%     save(wfile, 'x','y','u1','u2','times','n1','n2','x_left','x_right','y_left','y_right','dt');
    save(wfile, 'u1','u2','Du1Dt', 'Du2Dt','dt');

    % Use the values at t+dt as initial conditions for the next time step
    u1_vec = F(3,1:n1*n2)'; 
    u2_vec = F(3,n1*n2+1:2*n1*n2)';
    disp(t);
end

save(wfile, 'x','y','u1','u2','times','n1','n2','x_left','x_right','y_left','y_right','dt');
u1 = reshape(F(3,1:n1*n2),         n2, n1);
u2 = reshape(F(3,n1*n2+1:2*n1*n2), n2, n1);