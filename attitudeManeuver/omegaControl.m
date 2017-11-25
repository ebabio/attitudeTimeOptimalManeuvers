function u = omegaControl(t, x, parameters)

lNorm = parameters.lNorm;
k = parameters.k;

% Reassign variables
q = x(1:4);
omega = x(5:7);
lambda_q = x(8:11);
lambda_omega = x(12:14);

% Modified control
lambda_omega_exp = @(i) abs(lambda_omega(i))^(1/(lNorm-1));
den = (lambda_omega_exp(1)^lNorm + lambda_omega_exp(2)^lNorm + lambda_omega_exp(3)^lNorm)^(1/lNorm);

fun = @(x) sign(x);

u(1) = -fun(k*lambda_omega(1)) * lambda_omega_exp(1) / den;
u(2) = -fun(k*lambda_omega(2)) * lambda_omega_exp(2) / den;
u(3) = -fun(k*lambda_omega(3)) * lambda_omega_exp(3) / den;
