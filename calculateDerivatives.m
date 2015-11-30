function dYdt = calculateDerivatives(t,Y,modelInfo)

dYdt = transferComponents(t,Y,modelInfo) + advectQuantitiesWithVelocity(t,Y,modelInfo);
