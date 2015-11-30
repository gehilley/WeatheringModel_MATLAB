function dYdt = transferComponents(t,Y,modelInfo)

% Y is ordered solid phase first, fluid phase second.

numberOfNodes = modelInfo.numberOfNodes;
npars_solid = modelInfo.numberOfSolidComponents;
npars_fluid = modelInfo.numberOfFluidComponents;
dx = modelInfo.dx;

totalNumberOfSolidNodes = numberOfNodes .* npars_solid;
totalNumberOfFluidNodes = numberOfNodes .* npars_fluid;

Ysolid = reshape(Y(1:totalNumberOfSolidNodes),numberOfNodes,npars_solid);
Yfluid = reshape(Y(totalNumberOfSolidNodes+1:end),numberOfNodes,npars_fluid);

dYsoliddt = zeros(numberOfNodes,npars_solid);
dYfluiddt = zeros(numberOfNodes,npars_fluid);

% Loop through transfer functions to calculate dYdt:

for(i=1:modelInfo.numberOfTransferFunctions)
    
    thisTransferFunction = modelInfo.transferFunctions{i};
    
    if(strcmp(thisTransferFunction.outComponentType,'S'))
        
        outValues = Ysolid(:,thisTransferFunction.outComponentNumber).^(thisTransferFunction.stoichiometryIn./thisTransferFunction.stoichiometryOut);
        
        outValues = (Ysolid(:,thisTransferFunction.outComponentNumber) > 0).*outValues;
        
    else
        
        outValues = Yfluid(:,thisTransferFunction.outComponentNumber).^(thisTransferFunction.stoichiometryIn./thisTransferFunction.stoichiometryOut);
        
        outValues = (Yfluid(:,thisTransferFunction.outComponentNumber) > 0).*outValues;
        
    end
    
    k = thisTransferFunction.ks .* ones(numberOfNodes,1);
    
    % Loop through scaling functions:
    
    for(j=1:thisTransferFunction.numberOfScalingFunctions)
        
        thisScalingFunction = thisTransferFunction.scalingFunctions{j};
        
        if(strcmp(thisScalingFunction.componentType,'S'))
            
            scaleValues = Ysolid(:,thisScalingFunction.componentNumber);
            
        else
            
            scaleValues = Yfluid(:,thisScalingFunction.componentNumber);
            
        end
        
        kscale = thisScalingFunction.a.*scaleValues + thisScalingFunction.b.*exp(thisScalingFunction.c.*scaleValues) + thisScalingFunction.d;
        
        kscale = (kscale >= thisScalingFunction.f).*thisScalingFunction.f + ...
            ((kscale < thisScalingFunction.f) & ...
            (kscale > thisScalingFunction.e)).*kscale + ...
            (kscale <= thisScalingFunction.e).*thisScalingFunction.e;
        
        k = k.*kscale;
        
    end
    
    dYdtOutNoStoichiometry = -k.*outValues;
        
    dYdtOut = dYdtOutNoStoichiometry .* thisTransferFunction.stoichiometryOut;
    
    dYdtIn = - dYdtOutNoStoichiometry .* thisTransferFunction.stoichiometryIn;
    
    if(strcmp(thisTransferFunction.outComponentType,'S'))
        
        dYsoliddt(:,thisTransferFunction.outComponentNumber) = dYsoliddt(:,thisTransferFunction.outComponentNumber) + dYdtOut;
        
    else
        
        dYfluiddt(:,thisTransferFunction.outComponentNumber) = dYfluiddt(:,thisTransferFunction.outComponentNumber) + dYdtOut;
        
    end
    
    if(strcmp(thisTransferFunction.inComponentType,'S'))
        
        dYsoliddt(:,thisTransferFunction.inComponentNumber) = dYsoliddt(:,thisTransferFunction.inComponentNumber) + dYdtIn;
        
    elseif(strcmp(thisTransferFunction.inComponentType,'F'))
        
        dYfluiddt(:,thisTransferFunction.inComponentNumber) = dYfluiddt(:,thisTransferFunction.inComponentNumber) + dYdtIn;
        
    end
    
end

% Restack dYdt:

dYdt = [reshape(dYsoliddt,numberOfNodes.*npars_solid,1);reshape(dYfluiddt,numberOfNodes.*npars_fluid,1)];


        
        