function [BCs, pathCs] = attitudeBCs(x0BCs, x0, xfBCs, xf, xDotFun, tf)
BCs = zeros(15,1);
pathCs = zeros(2,1);

if(nargin == 6)
    collocation = 1;
else
    collocation = 0;
end

%% Allow for xDotFun collocation parameters
if(collocation)
    xDotWrapper = @(tau, x) (xDotFun(tau, x, tf));
else
    xDotWrapper = @(t,x) (xDotFun(t, x));
end

%% Initial state constraints
% Reassign variables
q0 = x0(1:4);
omega0 = x0(5:7);
lambda_q0 = x0(8:12);
lambda_omega0 = x0(13:14);

% Hamiltonian
xDot = xDotWrapper([],x0);
H0 = [lambda_q0;lambda_omega0]' * xDot(1:7);

BCs(1:7)  = [q0-x0BCs(1:4); omega0-x0BCs(5:7)];
pathCs(1) = H0+1;
BCs(15) = norm(lambda_q0)-1;

%% Final State Constraints
qf = xf(1:4);
omegaf = xf(5:7);
lambda_qf = xf(8:12);
lambda_omegaf = xf(13:14);

% Hamiltonian
xDot = xDotWrapper([],xf);
Hf = [lambda_qf;lambda_omegaf]' * xDot(1:7);

BCs(8:14)  = [qf-xfBCs(1:4); omegaf-xfBCs(5:7)];
pathCs(2) = Hf+1;
% BCs(15) = pathCs(2);
return