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
lambdaDot_omega = 1/2 * quatmultiply(lambda_q, quatconj(q)); % check derivation

% b. Element wise
qDotEW(2) = 1/2 * ( omega(1)*q(1)   - omega(2)*q(4)     + omega(3)*q(3));
qDotEW(3) = 1/2 * ( omega(1)*q(4)   + omega(2)*q(1)     - omega(3)*q(2));
qDotEW(4) = 1/2 * (-omega(1)*q(3)   + omega(2)*q(2)     + omega(3)*q(1));
qDotEW(1) = 1/2 * (-omega(1)*q(2)   - omega(2)*q(3)     - omega(3)*q(4));

lambdaDot_omegaEW(1) = - 1/2 * (lambda_q(2)*q(1)    + lambda_q(3)*q(4)  - lambda_q(4)*q(3)  - lambda_q(1)*q(2));
lambdaDot_omegaEW(2) = - 1/2 * (-lambda_q(2)*q(4)   + lambda_q(3)*q(1)  + lambda_q(4)*q(2)  - lambda_q(1)*q(3));
lambdaDot_omegaEW(3) = - 1/2 * (lambda_q(2)*q(3)    - lambda_q(3)*q(2)  + lambda_q(4)*q(1)  - lambda_q(1)*q(4));

lambdaDot_qEW(2) = -1/2 * (lambda_q(2)*0        -lambda_q(3)*omega(3)   + lambda_q(4)*omega(2)  - lambda_q(1)*omega(1));
lambdaDot_qEW(3) = -1/2 * (lambda_q(2)*omega(3) -lambda_q(3)*0          - lambda_q(4)*omega(1)  - lambda_q(1)*omega(2));
lambdaDot_qEW(4) = -1/2 * (-lambda_q(2)*omega(2)+lambda_q(3)*omega(1)   + lambda_q(4)*0         - lambda_q(1)*omega(3));
lambdaDot_qEW(1) = -1/2 * (lambda_q(2)*omega(1) +lambda_q(3)*omega(2)   + lambda_q(4)*omega(3)  - lambda_q(1)*0);

dxdt = [qDot omegaDot lambdaDot_q lambdaDot_omega(2:4)]';
return