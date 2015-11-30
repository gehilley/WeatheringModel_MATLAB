function dYdt = advectQuantitiesWithVelocity(t,Y,modelInfo)

% Y is ordered solid phase first, fluid phase second.

numberOfNodes = modelInfo.numberOfNodes;
npars_solid = modelInfo.numberOfSolidComponents;
npars_fluid = modelInfo.numberOfFluidComponents;
dx = modelInfo.dx;

totalNumberOfSolidNodes = numberOfNodes .* npars_solid;
totalNumberOfFluidNodes = numberOfNodes .* npars_fluid;

Ysolid = reshape(Y(1:totalNumberOfSolidNodes),numberOfNodes,npars_solid);
Yfluid = reshape(Y(totalNumberOfSolidNodes+1:end),numberOfNodes,npars_fluid);

Ysolid = [Ysolid;reshape(modelInfo.solidBoundaryConditions,1,npars_solid)];
Yfluid = [reshape(modelInfo.fluidBoundaryConditions,1,npars_fluid);Yfluid];

dYsoliddx = diff(Ysolid)./dx;
dYfluiddx = diff(Yfluid)./dx;

dYsoliddt = modelInfo.erosionRate.*dYsoliddx;
dYfluiddt = -modelInfo.waterBalance.*dYfluiddx;

dYdt = [reshape(dYsoliddt,prod(size(dYsoliddt)),1);reshape(dYfluiddt,prod(size(dYfluiddt)),1)];






