function rho = bt_rhoFcn(u,v)

nx = length(u(:,1));

rho = [zeros(nx,1) ones(nx,1).*v];



