function dxdt = eigenaxisODEs(~, x, parameters)
k = parameters;

% Reassign variables
theta = x(1);
omega = x(2);
lambda_theta = x(3);
lambda_omega = x(4);

% Synthethize control
u = -sign(k*lambda_omega);

% Compute derivatives
thetaDot = omega;
omegaDot = u;
lambdaDot_theta = 0;
lambdaDot_omega = -lambda_theta;

% Actual update
dxdt = [thetaDot; omegaDot; lambdaDot_theta; lambdaDot_omega];
return