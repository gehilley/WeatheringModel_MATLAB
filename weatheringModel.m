classdef weatheringModel < handle
    
    properties(SetAccess = private)
        
        dx = 0;
        numberOfNodes = 0;
        numberOfSolidComponents = 0;
        numberOfFluidComponents = 0;
        
        solidComponentNames = {};
        fluidComponentNames = {};
        
        solidComponents = {};
        fluidComponents = {};
        
        waterBalance = 0;
        erosionRate = 0;
        
        numberOfTransferFunctions = 0;
        transferFunctions = {};
        transferFunctionNames = {};
        
    end
    
    
    methods
        
        function WM = weatheringModel(n,dx,waterBalance,erosionRate)
            
            if(n <= 0 || n ~= floor(n))
                error('bioModel: Invalid number of nodes');
            end
            
            if(dx <= 0)
                error('bioModel: Invald node spacing.');
            end
                
            WM.numberOfNodes = n;
            WM.dx = dx;
            WM.waterBalance = waterBalance;
            WM.erosionRate = erosionRate;
            
        end
                
        function addSolidComponent(WM,solidComponentName,initialCondition,boundaryConditions)
            
            if(length(initialCondition) ~= 1 && length(initialCondition) ...
                    ~= WM.numberOfNodes)
                error(['addSolidComponent: Number of Values Specified in initialCondition is incompatible with ' num2str(WM.numberOfNodes) ' nodes in model.']);
            end
            
            WM.numberOfSolidComponents = WM.numberOfSolidComponents + 1;
            WM.solidComponentNames{WM.numberOfSolidComponents} = solidComponentName;
            
            thisComponent.name = solidComponentName;
            if(length(initialCondition) == 1)
                thisComponent.values = ones(WM.numberOfNodes,1).*initialCondition;
            else
                thisComponent.values = reshape(initialCondition,WM.numberOfNodes,1);
            end
            
            thisComponent.boundaryConditions = boundaryConditions;
            
            WM.solidComponents{WM.numberOfSolidComponents} = thisComponent;
            
        end
        
        function addFluidComponent(WM,fluidComponentName,initialCondition,boundaryConditions)
            
            if(length(initialCondition) ~= 1 && length(initialCondition) ...
                    ~= WM.numberOfNodes)
                error(['addFluidComponent: Number of Values Specified in initialCondition is incompatible with ' num2str(WM.numberOfNodes) ' nodes in model.']);
            end
            
            WM.numberOfFluidComponents = WM.numberOfFluidComponents + 1;
            WM.fluidComponentNames{WM.numberOfFluidComponents} = fluidComponentName;
            
            thisComponent.name = fluidComponentName;
            if(length(initialCondition) == 1)
                thisComponent.values = ones(WM.numberOfNodes,1).*initialCondition;
            else
                thisComponent.values = reshape(initialCondition,WM.numberOfNodes,1);
            end
            
            thisComponent.boundaryConditions = boundaryConditions;
            
            WM.fluidComponents{WM.numberOfFluidComponents} = thisComponent;
           
        end
        
        function addTransferFunction(WM,transferFunctionName,componentsIn,componentsOut,stoichiometriesIn,stoichiometriesOut,ks,useCollisionTheory)
            
            % Look for componentsOut and componentsIn in solidComponents:
            
            for(i=1:length(componentsIn))
                component = componentsIn(i);
                if(isempty(find(strcmp(WM.solidComponentNames,component))) && isempty(find(strcmp(WM.fluidComponentNames,component))))
                  error(['addTransferFunction: no component named ' componentOut]);
                end
            end
            
            for(i=1:length(componentsOut))
                component = componentsOut(i);
                if(isempty(find(strcmp(WM.solidComponentNames,component))) && isempty(find(strcmp(WM.fluidComponentNames,component))) && ~strcmp(component,'NULL'))
                   error(['addTransferFunction: no component named ' componentIn]);
                end
            end
            
            if((length(componentsIn) ~= length(stoichiometriesIn)) || (length(componentsOut) ~= length(stoichiometriesOut)))
                error('addTransferFunction: number of stoichiometric coefficients does not match number of components.');
            end
            
            thisTransferFunction = [];
            indexOfComponents = 1;
            for(i=1:length(componentsIn)) 
                component = componentsIn(i);
                thisComponent = [];
                thisComponent.stoichiometry = stoichiometriesIn(i);
                if(~isempty(find(strcmp(WM.solidComponentNames,componentIn))))
                    thisComponent.type = 'S';
                    thisComponent.componentNumber = find(strcmp(WM.solidComponentNames,component));
                elseif(~isempty(find(strcmp(WM.fluidComponentNames,componentIn))))
                    thisComponent.type = 'F';
                    thisComponent.componentNumber = find(strcmp(WM.fluidComponentNames,component));
                else
                    thisComponent.type = 'N';
                    thisComponent.componentNumber = 0;
                end
                
                thisTransferFunction.inComponents(indexOfComponents) = thisComponent;
                indexOfComponents = indexOfComponents + 1;
                
            end
            
            indexOfComponents = 1;
            
            for(i=1:length(componentsOut))
                component = componentsOut(i);
                thisComponent = [];
                thisComponent.stoichiometry = stoichiometriesOut(i);
                if(~isempty(find(strcmp(WM.solidComponentNames,componentOut))))
                    thisComponent.type = 'S';
                    thisComponent.componentNumber = find(strcmp(WM.solidComponentNames,component));
                elseif(~isempty(find(strcmp(WM.fluidComponentNames,componentOut))))
                    thisComponent.type = 'F';
                    thisComponent.componentNumber = find(strcmp(WM.fluidComponentNames,component));
                else
                    thisComponent.type = 'N';
                    thisComponent.componentNumber = 0;
                end
                
                thisTransferFunction.outComponents(indexOfComponents) = thisComponent;
                indexOfComponents = indexOfComponents + 1;
                
            end
            
            thisTransferFunction.ks = ks;
            thisTransferFunction.name = transferFunctionName;
            thisTransferFunction.numberOfScalingFunctions = 0;
            thisTransferFunction.scalingFunctions = {};
            thisTransferFunction.useCollisionTheory = useCollisionTheory;
            
            WM.numberOfTransferFunctions = WM.numberOfTransferFunctions + 1;
            WM.transferFunctions{WM.numberOfTransferFunctions} = thisTransferFunction;
            WM.transferFunctionNames{WM.numberOfTransferFunctions} = transferFunctionName;
            
        end
        
        function addScalingFunctionToTransferFunction(WM,transferFunctionName,scalingComponent,a,b,c,d,e,f)
            
            indexOfTransferFunction = find(strcmp(transferFunctionName,WM.transferFunctionNames));
            if(isempty(indexOfTransferFunction))
                error(['addScalingFunctionToTransferFunction: no transfer function called ' transferFunctionName ' exists.']);
            end
            
            if(isempty(find(strcmp(scalingComponent,WM.solidComponentNames))) && isempty(find(strcmp(scalingComponent,WM.fluidComponentNames))))
                error(['addScalingFunctionToTransferFunction: no component named ' scalingComponent ' exists.']);
            end
            
            thisTransferFunction = WM.transferFunctions{indexOfTransferFunction};
            thisTransferFunction.numberOfScalingFunctions = thisTransferFunction.numberOfScalingFunctions + 1;
            
            scalingFunction = [];
            
            if(~isempty(find(strcmp(scalingComponent,WM.solidComponentNames))))
                scalingFunction.componentType = 'S';
            else
                scalingFunction.componentType = 'F';
            end
            
            if(strcmp(scalingFunction.componentType,'S'))
                scalingFunction.componentNumber = find(strcmp(scalingComponent,WM.solidComponentNames));
            else
                scalingFunction.componentNumber = find(strcmp(scalingComponent,WM.fluidComponentNames));
            end
            
            scalingFunction.a = a;
            scalingFunction.b = b;
            scalingFunction.c = c;
            scalingFunction.d = d;
            scalingFunction.e = e;
            scalingFunction.f = f;
            
            thisTransferFunction.scalingFunctions{thisTransferFunction.numberOfScalingFunctions} = scalingFunction;
            
            WM.transferFunctions{indexOfTransferFunction} = thisTransferFunction;
            
            
        end
            
            
        function plotComponents(WM)
            
            xvec = -([0:WM.dx:((WM.numberOfNodes-1).*WM.dx)]);
            subplot(1,2,1)
            sc = [];
            for(i=1:WM.numberOfSolidComponents) 
                thiscomponent = WM.solidComponents{i};
                sc = [sc thiscomponent.values];
            end
            fc = [];
            for(i=1:WM.numberOfFluidComponents) 
                thiscomponent = WM.fluidComponents{i};
                fc = [fc thiscomponent.values];
            end
            
            plot(xvec,sc);
            view(90,-90);
            
            subplot(1,2,2)
            
            plot(xvec,fc);
            view(90,-90);
            
        end
        
        function WMout = integrateToTime(WM,time)
            
            modelInfo.numberOfNodes = WM.numberOfNodes;
            modelInfo.numberOfSolidComponents = WM.numberOfSolidComponents;
            modelInfo.numberOfFluidComponents = WM.numberOfFluidComponents;
            modelInfo.erosionRate = WM.erosionRate;
            modelInfo.waterBalance = WM.waterBalance;
            modelInfo.dx = WM.dx;
            modelInfo.numberOfTransferFunctions = WM.numberOfTransferFunctions;
            modelInfo.transferFunctions = WM.transferFunctions;
            modelInfo.solidComponentNames = WM.solidComponentNames;
            modelInfo.fluidComponentNames = WM.fluidComponentNames;
            
            Yo = [];
            
            for(i = 1:WM.numberOfSolidComponents)
                
              thisComponent = WM.solidComponents{i};
              Yo = [Yo;reshape(thisComponent.values,WM.numberOfNodes,1)];
              modelInfo.solidBoundaryConditions(i) = thisComponent.boundaryCondition;
              
            end
            
            for(i = 1:WM.numberOfFluidComponents)
                
                thisComponent = WM.fluidComponents{i};
                Yo = [Yo;reshape(thisComponent.values,WM.numberOfNodes,1)];
                modelInfo.fluidBoundaryConditions(i) = thisComponent.boundaryCondition;
                
            end
            
            fun = @(t,Y) calculateDerivatives(t,Y,modelInfo);

            [T,Y] = ode45(fun,[0 time],Yo);
            
            Yout = Y(end,:)';
            
            WMout = WM;
            
            numberOfSolidNodes = WM.numberOfNodes.*WM.numberOfSolidComponents;
            numberOfFluidNodes = WM.numberOfNodes.*WM.numberOfFluidComponents;
            
            YoutSolids = reshape(Yout(1:numberOfSolidNodes),WM.numberOfNodes,WM.numberOfSolidComponents);
            YoutFluids = reshape(Yout(numberOfSolidNodes+1:end),WM.numberOfNodes,WM.numberOfFluidComponents);
            
            for(i = 1:WM.numberOfSolidComponents)
                
                thisComponent = WM.solidComponents{i};
                thisComponent.values = YoutSolids(:,i);
                WMout.solidComponents{i} = thisComponent;
                
            end
            
            for(i=1:WM.numberOfFluidComponents)
                
                thisComponent = WM.fluidComponents{i};
                thisComponent.values = YoutFluids(:,i);
                WMout.fluidComponents{i} = thisComponent;
            
            end
            
        end
        
            
    end
    
end
