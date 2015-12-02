% Biotite Oxidation Simulation
% 
% Weathering Model example to demonstrate how to encode reactions and
% kinetics into weathering model.  Either a rule-based or a formal reactive
% transport approach may be employed due to model flexibility.

% Create the weathering model:

n = 100;
dx = 0.01;
waterBalance = 0.1;
erosionRate = 0;

waterBalance = waterBalance ./ (60*60*24*365);

WM = WM.weatheringModel(n,dx,waterBalance,erosionRate)

T = 3; % tortuosity
D = 3e-9;
D = D ./ (T.^2);

% Add components:

WM.addComponent(WM,'FeO','Solid', 181.5, 181.5);
WM.addComponent(WM,'O2','Fluid',2.6e-1,2.6e-1,D);
WM.addComponent(WM,'Fe(OH)3','Solid',0,0);
WM.addComponent(WM,'AllOxides','Solid',39171,39171);

% Create kinetic description:

k = 10^(-10.7);

KD = KD.makeEmptyKineticDescription(WM);
KD.addConstant(WM,k);

S = 4.7.*433.5;
KD.addLinearFunction(WM, 'AllOxides', S, 0);


r = 1/4;
% Fastest case:




S = S_m2_g.*433.5; %mol/m2
O2o = 2.6e-1; % mol/m3
FeOo = 181.5; %mol/m3
Oxide = 39171; % mol/m3
v = 0.1; % velocity in m/yr
D = 3e-9; % diffusivity in m2/s


dx = 0.01; % in m.
nx = 1000;

O2o_vec = ones(1,nx).*O2o;
FeOo_vec = ones(1,nx).*FeOo;

n = 2;




yinit = [FeOo_vec O2o_vec]';

uo = [FeOo O2o];

limiterFcn = @(u) superbee(u);

boundaryConditionFcn = @(u) bt_boundaryConditions(u,uo);
fluxFcn = @(u) bt_fluxFunction(u,v);
rhoFcn = @(u) bt_rhoFcn(u,v);
sourceFcn = @(u) bt_sourceFunction(u,k,S,r,Oxide);
diffusionFcn = @(u) bt_diffusionFunction(u,D,dx);

f = @(t,y) muscl(t,y,boundaryConditionFcn,fluxFcn,sourceFcn,rhoFcn,limiterFcn,diffusionFcn,n,dx);

tspan_years = 0:1:100;

tspan_s = tspan_years .* 60 .* 60 .* 24 .* 365;

[t,y] = ode45(f,tspan_s,yinit);

x = 0:dx:(nx-1).*dx;

FeO = y(:,1:nx);
O2 = y(:,nx+1:end);

subplot(3,1,1)

plot(x,FeO./FeOo)

subplot(3,1,2)

plot(x,O2./O2o)

subplot(3,1,3)
plot(tspan_years,O2(:,10)./O2o,'k-');
hold on
plot(tspan_years,FeO(:,10)./FeOo,'r--');


