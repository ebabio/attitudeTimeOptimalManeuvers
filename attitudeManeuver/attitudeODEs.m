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
% a. Quaternion multiplication way
qDot = 1/2 * quatmultiply(q, [0 omega]); % check derivation
omegaDot = u;
lambdaDot_q = 1/2 * quatmultiply(lambda_q, [0 omega]); % check derivation
lambdaDot_omega = -1/2 * quatmultiply(quatconj(q), lambda_q); % check derivation

dxdt = [qDot omegaDot lambdaDot_q lambdaDot_omega(2:4)]';
return