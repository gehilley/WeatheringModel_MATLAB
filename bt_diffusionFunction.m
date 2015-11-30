function dudt = bt_diffusionFunction(u,D,dx)

FeO = u(:,1);
O2 = u(:,2);

nx = length(FeO)-2;

dFeOdt = zeros(nx,1);

dO2dt = D.*(O2(3:end) - 2.*O2(2:end-1) + O2(1:end-2))./dx.^2;

dudt = [dFeOdt dO2dt];

