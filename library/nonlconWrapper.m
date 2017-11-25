function [g, h] = nonlconWrapper(x, shootingBCs, parameters)

%% Assign parameters
nIntervals = parameters.nIntervals;
setupOrbit = parameters.setupOrbit;

dim = (numel(x)-1)/nIntervals;
%% Evaluate orbits

tf = x(end);

iterOrbit(nIntervals) = Orbit();
tfVector = tf * (1:nIntervals)./nIntervals;
for i=1:nIntervals
    iterOrbit(i) = copy(setupOrbit);
    if(i==1)
        iterOrbit(i).x0 = x(1:dim);
        iterOrbit(i).tf = tfVector(i);
    else
        index0 = (i-1)*dim;
        iterOrbit(i).x0 = x(index0+1:index0+dim);
        iterOrbit(i).tf = tfVector(i) - tfVector(i-1);
    end
    iterOrbit(i).integrateX0();
end

%% Evaluate constraints

% Trajectory constraints
constraintErr = shootingBCs(iterOrbit(1).x0, iterOrbit(end).xf);

% Defects
defectErr = zeros(dim * (nIntervals-1),1);
for i  = 2:nIntervals
    index0 = dim*(i-2);
    defectErr(index0+1:index0+dim) = iterOrbit(i).x0 - iterOrbit(i-1).xf;
end

%% Wrap up for returning
g = 0;
h = [constraintErr; defectErr];