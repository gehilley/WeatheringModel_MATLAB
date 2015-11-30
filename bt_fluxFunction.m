function F = bt_fluxFunction(u,v)

FeO = u(:,1);
O2 = u(:,2);

nx = length(u);

F(:,1) = zeros(nx,1);
F(:,2) = O2.*v;

