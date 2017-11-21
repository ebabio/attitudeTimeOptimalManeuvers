function [BCs, pathCs] = eigenaxisBCs(x0BCs, x0, xfBCs, xf, xDotFun)
BCs = zeros(5,1);
pathCs = zeros(2,1);

%% Initial state constraints
theta0 = x0(1);
omega0 = x0(2);
lambda0 = x0(3:4);

% Initial Hamiltonian
xDot = xDotFun([],x0);
H0 = lambda0' * xDot(1:2);

BCs(1:2)  = [theta0-x0BCs(1); omega0-x0BCs(2)];
pathCs(1) = H0+1;
%% Final state constraints
% Reassign variables
thetaf = xf(1);
omegaf = xf(2);
lambdaf = xf(3:4);

% Final Hamiltonian
xDot = xDotFun([],xf);
Hf = lambdaf' * xDot(1:2);  

BCs(3:4) = [thetaf-xfBCs(1); omegaf-xfBCs(2)];
pathCs(2) = Hf+1;
BCs(5) = pathCs(2);
return