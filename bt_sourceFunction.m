function dudt = bt_sourceFunction(u,k,S,r,Oxide)

FeO = u(:,1);
O2 = u(:,2);

nx = length(FeO);

dFeOdt = (O2 >= 0).*-0.5.*k.*S.*FeO.^2.*O2.^(r) ./ Oxide;

dO2dt = (O2 >= 0).*(-0.5.*r.*k.*S.*(FeO).^2.*O2.^(r) ./ Oxide);

dudt = [dFeOdt dO2dt];

