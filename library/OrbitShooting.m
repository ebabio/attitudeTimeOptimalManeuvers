classdef OrbitShooting < matlab.mixin.Copyable
    
    properties
        maxIter
        epsilon
        
        nIntervals
        
        gEqualityConstraint
        variationMap
        
        orbit
        count
        k
    end
    
    methods
        function obj = OrbitShooting()
            obj.maxIter = 30;
            obj.epsilon = 1e-12;
            obj.k = .5;
        end
        
        function [orbit, count] = target(obj, gEqualityConstraint, initialOrbit, variationMap)
            obj.gEqualityConstraint = gEqualityConstraint;
            obj.variationMap = variationMap;
            
            % Setup orbits by Intervals
            iterOrbit(obj.nIntervals) = Orbit();
            tfVector = initialOrbit.tf * (1:obj.nIntervals)./obj.nIntervals;
            for i=1:obj.nIntervals
                iterOrbit(i) = copy(initialOrbit);
                if(i==1)
                    iterOrbit(i).tf = tfVector(i);
                else
                    iterOrbit(i).x0 = deval(initialOrbit.odeSol, tfVector(i-1));
                    iterOrbit(i).tf = tfVector(i) - tfVector(i-1);
                end
                iterOrbit(i).integrateX0();
            end
            
            % First iteration
            error = 1;
            obj.count = 0;
            
            % Iterate
            while(error > obj.epsilon)
                obj.count = obj.count+1;
                
                [hessian, psi] = obj.computeHessian(iterOrbit);
                
                tolerance = error*(10^-(2+obj.count/10));
                dx0 = -obj.k*pinv(hessian, tolerance)*psi; %lstsqr
                
                for i=1:obj.nIntervals
                    dim = numel(iterOrbit(i).x0);
                    dx0i = dx0((i-1)*dim+1: i*dim);
                    if(i==obj.nIntervals)
                        dx0i(dim+1) = dx0(end);
                    end
                    iterOrbit(i).addDeltaX0(dx0i);
                end
                
                error = norm(psi);
                disp(['iteration ' num2str(obj.count) ' error:' num2str(error)]);
                
                if(iterOrbit(end).tf < 0)
                    error('time for the last interval is negative');
                end
                
                if(obj.count > obj.maxIter)
                    obj.count = -1;
                    break
                end
            end
            
            if(obj.count == -1)
                warning('Max number of iterations reached')
            else
                disp('Boundary conditions satisfied')
            end
            
            % Prepare final orbit
            if(obj.nIntervals == 1)
                tPrevious =0;
            else
                tPrevious = tfVector(end-1);
            end
            obj.orbit = copy(iterOrbit(1));
            obj.orbit.tf = tPrevious+iterOrbit(end).tf;
            obj.orbit.integrateX0();
            orbit = copy(obj.orbit);
            count = obj.count;
        end
        
        
        function [hessian, psi0] = computeHessian(obj, orbit)
            % STM found by variation of the initial conditions
            dim = numel(orbit(1).x0);
            auxConst = obj.gEqualityConstraint(orbit(1).x0, orbit(end).xf);
            nBCs = numel(auxConst);
            
            % transform from indexes to an indicator vector
            withTime = find(obj.variationMap==dim+1);
            activeDimensions(obj.variationMap) = ones(1,numel(obj.variationMap));
            activeDimensions = activeDimensions(1:dim);
            
            hessian = zeros(nBCs+dim*(obj.nIntervals-1),dim*obj.nIntervals+1);
            delta = 1e-3;
            
            % Reference solution
            for i=1:obj.nIntervals
                orbit(i).integrateX0(orbit(i).tf);
            end
            
            constraintErr0 = obj.gEqualityConstraint(orbit(1).x0, orbit(end).xf);
            continuityErr0 = obj.continuityError(orbit);
            psi0 = [constraintErr0; continuityErr0];
            
            % Compute STM
            for j=1:obj.nIntervals
                if(j==1)
                    dimensions = find(activeDimensions==1);
                else
                    dimensions = 1:dim;
                end
                
                index0 = dim*(j-1);
                for i=dimensions;
                    variedOrbit = copy(orbit);
                    perturbation = delta * circshift([1; zeros(dim-1,1)], i-1);
                    variedOrbit(j).x0 = variedOrbit(j).x0 + perturbation;
                    variedOrbit(j).integrateX0(variedOrbit(j).tf);
                    
                    constraintErr = obj.gEqualityConstraint(variedOrbit(1).x0, variedOrbit(end).xf);
                    continuityErr = obj.continuityError(variedOrbit);
                    psi = [constraintErr; continuityErr];
                    hessian(:,index0+i) = (psi - psi0)./delta;
                end
            end
            
            
            if(withTime)
                variedOrbit = copy(orbit);
                variedOrbit(end).integrateX0(variedOrbit(end).tf+delta);
                
                constraintErr = obj.gEqualityConstraint(variedOrbit(1).x0, variedOrbit(end).xf);
                hessian(1:nBCs,end) = (constraintErr - constraintErr0)./delta;
            end
        end
        
        function fVec = continuityError(obj, orbit)
            dim = numel(orbit(1).x0);
            
            fVec = zeros(dim * (obj.nIntervals-1),1);
            
            interval = 2;
            while(interval <= obj.nIntervals)
                index0 = dim*(interval-2);
                fVec(index0+1:index0+dim) = orbit(interval).x0 - orbit(interval-1).xf;
                interval = interval+1;
            end
        end
    end
end