classdef Orbit < matlab.mixin.Copyable
    % Important note on handle classes: handle objects behave as handles/pointers
    % To create a different object with the same values use COPY method
    
    properties
        odefun
        
        x0
        
        x
        t
        
        xf
        tf
        
        odeSol
        odeOptions
        odeIntegrator
    end
    
    methods
        function obj = Orbit(x0, odefun)
            if(nargin ==2)
                obj.odefun = odefun;
                obj.x0 = x0;
                obj.x = obj.x0;
                obj.t = 0;
            elseif(nargin ==0)
                %empty constructor
            else
                warning('number of inputs not recognized');
            end
        end
        
        function addDeltaX0(obj, deltaX0)
            if(size(deltaX0,1)>numel(obj.x0))   %update vector may contain final time
                obj.x0 = obj.x0 + deltaX0(1:numel(obj.x0));
                obj.tf = obj.tf + deltaX0(numel(obj.x0)+1);
            else
                obj.x0 = obj.x0 + deltaX0;
            end
        end
        
        function obj = integrateX0(obj, tEnd)
            if(nargin < 2)
                tEnd = obj.tf;
            end
            
            obj.odeSol = ode45(obj.odefun, [0 tEnd], obj.x0, obj.odeOptions);
            obj.t = obj.odeSol.x;
            obj.x = obj.odeSol.y;
            
            obj.tf = obj.t(end);
            obj.xf = obj.x(:,end);
        end
        
        function phi = STM(obj, tEnd, dimensions)
            % STM found by variation of the initial conditions
            n = numel(obj.x0);
            
            if(nargin < 2)
                tEnd = obj.tf;
            end
            if(nargin < 3)
                activeDimensions = ones(1,n);
                withTime = 0;
            else
                % transform from indexes to an indicator vector
                withTime = find(dimensions==n+1);
                activeDimensions(dimensions) = ones(1,numel(dimensions));
                activeDimensions = activeDimensions(1:n);
            end
            
            phi = eye(n);
            delta = 1e-3;
            
            % Reference solution
            obj.odeSol = ode45(obj.odefun, [0 tEnd], obj.x0, obj.odeOptions);
            obj.t = obj.odeSol.x;
            obj.x = obj.odeSol.y;
            
            obj.tf = obj.t(end);
            obj.xf = obj.x(:,end);
            
            % Compute STM
            for i=find(activeDimensions==1)
                perturbation = delta * circshift([1; zeros(n-1,1)], i-1);
                xDelta0 = obj.x0 + perturbation;
                [~,xDeltaf] = ode45(obj.odefun, [0 tEnd], xDelta0, obj.odeOptions);
                phi(:,i) = (xDeltaf(end,:)' - obj.xf)/delta;
            end
            for i=find(activeDimensions==0)
                phi(:,i) = zeros(n,1);
            end
            
            if(withTime)
                phi(:,end+1) = obj.odefun([], obj.xf);
            end
        end
    end
    
end