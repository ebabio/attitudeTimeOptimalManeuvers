function dxdt = attitudeODEs(t, x, control)

% Reassign variables
x = x';         %since quaternion and row vectors in matlab
q = x(1:4);
omega = x(5:7);
lambda_q = x(8:11);
% lambda_omega = x(12:14);

% Synthethize control
u = control(t,x);

% Compute derivatives
qDot = 1/2 * quatmultiply([0 omega], q);
omegaDot = u;
lambdaDot_q = 1/2 * quatmultiply([0 omega], lambda_q);
lambdaDot_omega = -1/2 * quatmultiply(quatconj(q), lambda_q);

dxdt = [qDot omegaDot lambdaDot_q lambdaDot_omega(2:4)]';
return