classdef weatheringModel < handle
    
    properties(Access = private)
        
        dx = 0;
        numNodes = 0;
           
        components = {};
                
        waterBalance = 0;
        erosionRate = 0;
        
        reactions = {};
        
        numSolidComponents = 0;
        numFluidComponents = 0;
        numConstantComponents = 0;
        
        numReactions = 0;
        
        indexesOfSolidComponents = [];
        indexesOfFluidComponents = [];
        indexesOfConstantComponents = [];
        
    end
    
    methods(Access=private)
    
        function u = getUMatrixFromCurrentValues(WM)
            
            numComponents = length(WM.components);
            u = zeros(WM.numberOfNodes(),numComponents);
            
            for(i=1:numComponents)
                u(:,i) = WM.components{i}.values;
            end
            
        end
        
        function [RBC,LBC] = getBCs(WM)
            
            numComponents = length(WM.components);
            RBC = zeros(2,numComponents);
            LBC = zeros(2,numComponents);
            
            for(i=1:numComponents)
                RBC(:,i) = repmat(WM.components{i}.boundaryConditions(1),2,1);
                LBC(:,i) = repmat(WM.components{i}.boundaryConditions(2),2,1);
            end
        end
        
        function v = getVelocityMatrix(WM)
            
            numComponents = length(WM.components);
            v = zeros(WM.numberOfNodes(),numComponents);
            for(i=1:numComponents)
                if(strcmp(WM.components{i}.type,'Fluid'))
                    v(:,i) = ones(WM.numberOfNodes(),1).*WM.waterBalance;
                elseif(strcmp(WM.components{i}.type,'Solid'))
                    v(:,i) = -ones(WM.numberOfNodes(),1).*WM.erosionRate;
                end
                
            end
        end
           
        function mat = getMatrixFromVector(WM,y)
            
            nx = WM.numberOfNodes();
            numComponents = length(WM.components);
            
            mat = reshape(y,nx,numComponents);
            
        end
        
        function y = getVectorFromMatrix(WM,mat)
            
            nx = WM.numberOfNodes();
            numComponents = length(WM.components);
            
            y = reshape(mat,nx.*numComponents,1);
            
        end
        
        function limiter = superbee(WM, u)
            
            n = length(WM.numComponents);
            nx = length(WM.numberOfNodes());
            
            zer = zeros(nx,n);
            one = ones(nx,n);
            twos = 2.*ones(nx,n);
            
            % Precondition u for nans and infs:
            
            i = find(isnan(u));
            u(i) = -2;
            
            i = find(isinf(u) & (u > 0));
            u(i) = 2;
            i = find(isinf(u) & (u < 0));
            u(i) = -2;
            
            term1 = zer;
            term2 = (2.*u > one).*one + (2.*u <= one).*u;
            term3 = (u > twos).*twos + (u <= twos).*u;
            
            limiter1 = (term1 > term2).*term1 + (term1 <= term2).*term2;
            limiter = (limiter1 > term3).*limiter1 + (limiter1 <= term3).*term3;


        end
        
        function S = calculateSourceFunction(WM, y)
            
            mat = WM.getMatrixFromVector(y);
            
            nx = WM.numberOfNodes();
            numComponents = length(WM.components);
            
            S = zeros(nx,numComponents);
            
            % Loop through reactions:
            
            for(i=1:WM.numReactions)
                
                thisReaction = WM.reaction{i};
                
                % calculate kinetics:
                
                k = zeros(nx,1);
                
                kinetics = thisReaction.kinetics;
                
                if(length(thisReaction.kinetics) == 1)
                    kineticDescriptions{1} = kinetics;
                else
                    kineticDescriptions = kinetics;
                end
                
                for(j=1:length(kineticDescriptions))
                    
                    thisKineticDescription = kineticDescriptions{j};
                    
                    k = k + thisKineticDescription.evaluateKineticsWithWeatheringModel(WM, mat);
                    
                end
                
                for(j = 1:length(thisReaction.reactants))
                    
                    thisReactant = thisReaction.reactants{j};
                    thisComponent = thisReactant.component;
                    thisStoichiometry = thisReactant.stoichiometry;
                    
                    [index, values] = WM.getIndexAndValuesForComponent(thisComponent);
                    
                    S(:,index) = -k.*thisStoichiometry;
                end
                
                for(j = 1:length(thisReaction.products))
                    
                    thisProduct = thisReaction.products{j};
                    thisComponent = thisProduct.component;
                    thisStoichiometry = thisProduct.stoichiometry;
                    
                    [index, values] = W>getIndexAndValuesForComponent(thisComponent);
                    
                    S(:,index) = k.*thisStoichiometry;
                    
                end
                
                
            end
            
            S = WM.getVectorFromMatrix(S);
            
        end
        
        function D = calculateDiffusionFunction(WM, y)
            
            mat = WM.getMatrixFromVector(y);
            [UBC,BBC] = WM.getBCs();
            
            nx = WM.numberOfNodes();
            numComponents = length(WM.components);
            
            D = zeros(nx,numComponents);
            
            [indexes, names] = WM.getIndexesAndNamesForComponentType('Fluid');
            
            for(i=1:length(names))
                
                componentName = names{i};
                [indexes, values] = WM.getIndexAndValuesForComponent(componentName);
                values = [UBC(1);values;LBC(1)];
                Diffusivity = WM.components{i}.diffusivity;
                D(:,indexes) = Diffusivity.*(values(3:end) - 2.*values(2:end-1) + values(1:end-2)) ./ WM.dx;
            end
            
            D = WM.getVectorFromMatrix(D);
            
        end
        
        function F = calculateFluxes(WM, y)
            
            mat = WM.getMatrixFromVector(y);
            
            v = WM.getVelocityMatrix();
            
            % Add BCs:

            [UBC,BBC] = WM.getBCs;

            mat = [UBC;mat;BBC];

            % Calculate edge values:

            ui = mat(3:end-2,:);
            uip1 = mat(4:end-1,:);
            uip2 = mat(5:end,:);
            uim1 = mat(2:end-3,:);

            ri = (u(2:end-1,:) - u(1:end-2,:))./(u(3:end,:) - u(2:end-1,:));

            limiter = WM.superbee(ri);

            uip12L = ui + 0.5.*limiter(2:end-1,:).*(uip1-ui);
            uip12R = uip1 - 0.5.*limiter(3:end,:).*(uip2-uip1);
            uim12L = uim1 + 0.5.*limiter(1:end-2,:).*(ui-uim1);
            uim12R = ui - 0.5.*limiter(2:end-1,:).*(uip1-ui);

            rhoi = v;
            rhoim1 = v;
            rhoip1 = v;

            aim12 = (rhoi > rhoim1).*rhoi + (rhoi <= rhoim1).*rhoim1;
            aip12 = (rhoi > rhoip1).*rhoi + (rhoi <= rhoip1).*rhoip1;

            Fim12R = uim12R.*v;
            Fim12L = uim12L.*v;
            Fip12R = uip12R.*v;
            Fip12L = uip12L.*v;

            Fim12 = 0.5.*(Fim12R + Fim12L - aim12.*(uim12R - uim12L));
            Fip12 = 0.5.*(Fip12R + Fip12L - aip12.*(uip12R - uip12L));
            
            F = -(Fip12 - Fim12)./WM.dx;
            
            F = WM.getVectorFromMatrix(F);
            
        end
        
        function dydt = calculateDyDt(WM, t, y)
            
            S = WM.calculateSourceFunction(WM, y);
            D = WM.calculateDiffusionFunction(WM, y);
            F = WM. calculateFluxes(WM,y);
            
            dydt = S + D + F;
            
        end
            
        function u = getICMatrix(WM)
            
            nx = WM.numberOfNodes();
            numComponents = length(WM.components);
            
            u = zeros(nx, numComponents);
            
            for(i=1:length(WM.components))
                u(:,i) = WM.components{i}.values;
            end
            
        end
        
        function maxDt = CFLTimeStep(WM)
            
            maxDt = WM.dx ./ max([abs(WM.waterBalance) abs(WM.erosionRate)]);
        
        end

    end
    
    
    methods
        
        function WM = weatheringModel(n,dx,waterBalance,erosionRate)
            
            if(n <= 0 || n ~= floor(n))
                error('bioModel: Invalid number of nodes');
            end
            
            if(dx <= 0)
                error('bioModel: Invald node spacing.');
            end
                
            WM.numNodes = n;
            WM.dx = dx;
            WM.waterBalance = waterBalance;
            WM.erosionRate = erosionRate;
            
        end
                
        function number = numberOfSolidComponents(WM)
            number = WM.numSolidComponents;
        end
        
        function number = numberOfFluidComponents(WM)
            number = WM.numFluidComponents;
        end
        
        function number = numberOfConstantComponents(WM)
            number = WM.numConstantComponents;
        end
        
        function number = numberOfReactions(WM)
            number = WM.numReactions;
        end
        
        function number = numberOfNodes(WM)
            number = WM.numberOfNodes;
        end
        
        function test = componentExists(WM,componentName)
            test = 0;
            for(i=1:length(WM.components))
                if(strcmp(WM.components{i}.name, componentName))
                    test = 1;
                    return;
                end
            end
            
        end
            
        function addComponent(WM,componentName,componentType, initialCondition,boundaryConditions, diffusivity)
            
            if(length(initialCondition) ~= 1 && length(initialCondition) ...
                    ~= WM.numberOfNodes)
                error(['addComponent: Number of Values Specified in initialCondition is incompatible with ' num2str(WM.numberOfNodes) ' nodes in model.']);
            end
            
            if(WM.componentExists(componentName)) 
                error('addComponent: Component already exists.');
            end
            
            if((componentType ~= 'Solid') && (componentType ~= 'Fluid') && (componentType ~= 'Constant'))
                error('addComponent: Unrecognized component type.');
            end
            
            if((componentType ~= 'Fluid') && exists(diffusivity))
                error('addComponent: Diffusivity of components only applies to Fluid phase.');
            end
            
            if(~exists(diffusivity))
                diffusivity = 0;
            end
            
            thisComponent.name = componentName;
            thisComponent.type = componentType;
            
            if(length(initialCondition) == 1)
                thisComponent.values = ones(WM.numberOfNodes,1).*initialCondition;
            else
                thisComponent.values = reshape(initialCondition,WM.numberOfNodes,1);
            end
            
            thisComponent.boundaryConditions = boundaryConditions;
            indexOfThisComponent = length(WM.components)+1;
            
            if(strcmp(componentType,'Solid'))
                WM.numSolidComponents = WM.numSolidComponents + 1;
                WM.indexesOfSolidComponents = [WM.indexesOfSolidComponents, indexOfThisComponent];
            elseif(strcmp(componentType,'Fluid'))
                WM.numFluidComponents = WM.numFluidComponents + 1;
                WM.indexesOfFluidComponents = [WM.indexesOfFluidComponents, indexOfThisComponent];
                thisComponent.diffusivity = diffusivity;
            elseif(strcmp(componentType, 'Constant'))
                WM.numConstantComponents = WM.numConstantComponents + 1;
                WM.indexesOfConstantComponents = [WM.indexesOfConstantComponents, indexOfThisComponent];
            end
            
            WM.components{indexOfThisComponent} = thisComponent;
            
        end
        
        function addReaction(WM, reactants, reactantStoichiometries, products, productStoichiometries, kinetics)
            
            if( (length(reactants) ~= length(reactantStoichiometries)) || (length(reactantStoichiometries) ~= length(productStoichiometries)) )
                error('addReaction: number of components must match number of stoichiometric coefficients.');
            end
            
            for(i=1:length(reactants))
                if(~WM.componentExists(reactants{i}))
                    error('addReaction: reactant component does not exist.');
                end
            end
            
            for(i=1:length(products))
                if(~WM.componentExists(products{i}))
                    error('addReaction: product component does not exist.');
                end
            end
            
            thisReaction = {};
            
            for(i=1:length(reactants))
                thisReactant.stoichiometry = reactantStoichiometries(i);
                thisReactant.component = reactants{i};
                thisReaction.reactants{i} = thisReactant;
            end
            
            for(i=1:length(products))
                thisProduct.stoichiometry = productStoichiometries(i);
                thisProduct.component = productStoichiometries{i};
                thisReaction.products{i} = thisProduct;
            end
            
            thisReaction.kinetics = kinetics;
            
            WM.numReactions = WM.numReactions + 1;
            WM.reactions{WM.numReactions} = thisReaction;
        
        end
        
        function [index, values] = getIndexAndValuesForComponent(WM, componentName)
            if(~WM.componentExists(componentName))
                error('getValuesForComponent: Component does not exist.');
            end
            
            for(i=1:length(WM.components))
                if(strcmp(WM.components{i}.name,componentName))
                    index = i;
                    values = WM.components{i}.values;
                    return;
                end
            end
        end
        
        function [indexes, names] = getIndexesAndNamesForComponentType(WM,componentType)
            
            if((componentType ~= 'Solid') && (componentType ~= 'Fluid') && (componentType ~= 'Constant'))
                error('getIndexesAndNamesForComponentType: Unrecognized component type.');
            end
            
            indexes = [];
            names = {};
            counter = 1;
            for(i=1:length(WM.components))
                if(strcmp(WM.components{i}.type, componentType))
                    indexes(counter) = i;
                    names{counter} = WM.components{i}.name;
                    counter = counter + 1;
                end
            end
            
        end
                
        function yout = solveForTime(WM, tvec)
            
            yinit = WM.getICs();
            
            % Setup ODE solver (implicit):
            
            options = odeset('MaxStep',WM.CFLTimeStep());
            
            [tout, yout] = ode23s(WM.calculateDyDt,tvec,yinit,options);
            
        end             
            
    end
    
end
