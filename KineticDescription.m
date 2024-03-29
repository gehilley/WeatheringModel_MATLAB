classdef KineticDescription
    %KineticDescription Class that defines the kinetics of a particular
    %reaction.
    %   This class is used to create a kinetic description that is
    %   teathered to a reaction.  It is defined to be flexible enough to
    %   allow a rule-based empirical description of component transfer or a
    %   rigorous reactive-transport-based framework.
    
    properties(Access = private)
        
        components = {};
        
    end
    
    methods
        
        function KD = KineticDescription(WM)
            KD.components = {};
        end
        
        function number = numberOfComponents(KD)
            number = length(KD.components);
        end
        
        function KD = addConstant(KD,WM,a)
            
            thisComponent.type = 'Constant';
            thisComponent.value = a;
            KD.components{KD.numberOfComponents()+1} = thisComponent;
            
        end
        
        function KD = addPowerFunction(KD, WM, componentName, exponent)
            
            if(~WM.componentExists(componentName))
                error('addPowerFunction: no component with componentName exists in weathering model.');
            end
            
            thisComponent.type = 'PowerFunction';
            thisComponent.exponent = exponent;
            thisComponent.name = componentName;
            
            KD.components{KD.numberOfComponents()+1} = thisComponent;
            
        end
        
        function KD = addLinearFunction(KD, WM, componentName, a, b)
            
            if(~WM.componentExists(componentName))
                error('addLinearFunction: no component with componentName exists in weathering model.');
            end
            
            thisComponent.type = 'LinearFunction';
            thisComponent.name = componentName;
            thisComponent.a = a;
            thisComponent.b = b;
            
            KD.components{KD.numberOfComponents()+1} = thisComponent;
            
        end
        
        function KD = addExponentialFunction(KD, WM, componentName, foldingValue, scale, offset)
            
            if(~WM.componentExists(componentName))
                error('addExponentialFunction: no component with componentName exists in weathering model.');
            end
            
            thisComponent.type = 'ExponentialFunction';
            thisComponent.name = componentName;
            thisComponent.lambda = foldingValue;
            thisComponent.b = offset;
            thisComponent.a = scale;
            
            KD.components{KD.numberOfComponents()+1} = thisComponent;
            
        end
        
        function k = evaluateKineticsWithWeatheringModel(KD, WM, mat)
            
            k = ones(WM.numberOfNodes(),1);
            
            for(i=1:length(KD.components))
                
                thisComponent = KD.components{i};
                
                switch thisComponent.type
                    case 'Constant'
                        k = k.*thisComponent.value;
                    case 'PowerFunction'
                        scaledComponentName = thisComponent.name;
                        [index, values] = WM.getIndexAndValuesForComponent(scaledComponentName);
                        values = (mat(:,index) > 0).*(mat(:,index));
                        k = k.*values.^thisComponent.exponent;
                    case 'LinearFunction'
                        scaledComponentName = thisComponent.name;
                        [index, values] = WM.getIndexAndValuesForComponent(scaledComponentName);
                        values = (mat(:,index) > 0).*(mat(:,index));
                        values = values.*thisComponent.a + thisComponent.b;
                        values = (values >= 0).*values;
                        k = k.*values;
                    case 'ExponentialFunction'
                        scaledComponentName = thisComponent.name;
                        [index, values] = WM.getIndexAndValuesForComponent(scaledComponentName);
                        values = (mat(:,index) > 0).*(mat(:,index));
                        values = thisComponent.a.*(exp(thisComponent.lambda.*values) - thisComponent.b);
                        k = k.*values;
                end
                
            end
            
            
        end
        
        
    end
    
end

