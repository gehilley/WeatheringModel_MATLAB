function dydt = muscl(t,y,boundaryConditionFcn,fluxFcn,sourceFcn,rhoFcn,limiterFcn,diffusionFcn,n,dx)

t / (60*60*24*365)

nx = length(y) ./ n;

u = reshape(y,nx,n);
us = u;

% Add BCs:

[LBC,RBC] = boundaryConditionFcn(u);

u = [LBC;u;RBC];

% Calculate edge values:

ui = u(3:end-2,:);
uip1 = u(4:end-1,:);
uip2 = u(5:end,:);
uim1 = u(2:end-3,:);

ri = (u(2:end-1,:) - u(1:end-2,:))./(u(3:end,:) - u(2:end-1,:));

limiter = limiterFcn(ri);

uip12L = ui + 0.5.*limiter(2:end-1,:).*(uip1-ui);
uip12R = uip1 - 0.5.*limiter(3:end,:).*(uip2-uip1);
uim12L = uim1 + 0.5.*limiter(1:end-2,:).*(ui-uim1);
uim12R = ui - 0.5.*limiter(2:end-1,:).*(uip1-ui);

rhoi = rhoFcn(ui);
rhoim1 = rhoFcn(uim1);
rhoip1 = rhoFcn(uip1);

aim12 = (rhoi > rhoim1).*rhoi + (rhoi <= rhoim1).*rhoim1;
aip12 = (rhoi > rhoip1).*rhoi + (rhoi <= rhoip1).*rhoip1;

Fim12R = fluxFcn(uim12R);
Fim12L = fluxFcn(uim12L);
Fip12R = fluxFcn(uip12R);
Fip12L = fluxFcn(uip12L);

Fim12 = 0.5.*(Fim12R + Fim12L - aim12.*(uim12R - uim12L));
Fip12 = 0.5.*(Fip12R + Fip12L - aip12.*(uip12R - uip12L));

S = sourceFcn(us);

D = diffusionFcn(u(2:end-1,:));

dudt = -(Fip12 - Fim12)./dx + S + D;

dydt = reshape(dudt,nx*n,1);


