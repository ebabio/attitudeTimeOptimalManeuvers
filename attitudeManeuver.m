%% Eigenaxis Test

%% Setup workspace
addpath(genpath(pwd))
clear all
clc

%% Define boundary conditions

% Initial / final conditions
roll = 0;
pitch  = deg2rad(1);
yaw = deg2rad(180);

q0 = [1 0 0 0]';
qf = eul2quat([yaw pitch roll])';

omegaBCs = zeros(3,1);

x0BCs = [q0; omegaBCs];
xfBCs = [qf; omegaBCs];

parameters.lNorm = 2;
parameters.k = 1e2;

% First guess
thetaf = 2*acos(qf(1));
n = qf(2:4)/norm(qf(2:4));
tf = 2*sqrt(thetaf);

lambda_q = -[0; n];
lambda_omega = -tf/4*n;
lambda0 = [lambda_q; lambda_omega];

%% Setup

% 1. Trajectory ODEs
control = @(t,x) omegaControl(t, x, parameters);
realOde = @(t, x) (attitudeODEs(t, x, control));
collocationOde = @(tau, x, tf) (tf*realOde(tau,x));

% 2. Boundary Conditions
collocationBCs = @(x0, xf, tf) (attitudeBCs(x0BCs, x0, xfBCs, xf, collocationOde, tf));

% 3. Compute initial trajectory
x0 = [ x0BCs; lambda0];
orbit = Orbit(x0, realOde);
orbit.odeOptions = odeset('RelTol',1e-7 ,'AbsTol',1e-10);
orbit.integrateX0(tf);

% 4. Collocation algorithm
Nt = 20;
tau = linspace(0,1,Nt)';

solinit = bvpinit(tau, x0, tf);
solinit.x = orbit.odeSol.x;
solinit.y = orbit.odeSol.y;
solinit.parameters = tf;

options = bvpset('NMax',1e3, 'RelTol',1e-2 ,'AbsTol',1e-3); % default: ('NMax',1e3, 'RelTol',1e-3 ,'AbsTol',1e-6);

%% Solution

sol3dattitude = bvp4c(collocationOde, collocationBCs, solinit, options);

%% Evaluate the solution at all times

% Initial guess
maneuver0.t = solinit.x;
maneuver0.x = solinit.y;

maneuver0.q = maneuver0.x(1:4,:);
maneuver0.eul = quat2eul(maneuver0.q')';
maneuver0.omega = maneuver0.x(5:7,:);
maneuver0.lambda_q =  maneuver0.x(8:11,:);
maneuver0.lambda_omega =  maneuver0.x(12:14,:);

maneuver0.u = zeros(3, numel(maneuver0.t));
for i=1:numel(maneuver0.t)
    maneuver0.u(:,i) = control(maneuver0.t, maneuver0.x(:,i));
end

% Final maneuver
maneuver.tf = sol3dattitude.parameters(1);
maneuver.t = sol3dattitude.x * maneuver.tf;
maneuver.x = sol3dattitude.y;

maneuver.q = maneuver.x(1:4,:);
maneuver.eul = quat2eul(maneuver.q')';
maneuver.omega = maneuver.x(5:7,:);
maneuver.lambda_q =  maneuver.x(8:11,:);
maneuver.lambda_omega =  maneuver.x(12:14,:);

maneuver.u = zeros(3, numel(maneuver.t));
for i=1:numel(maneuver.t)
    maneuver.u(:,i) = control(maneuver.t, maneuver.x(:,i));
end

%% Plot data

% Initial maneuver
figure(1)
clf reset

subplot(2,2,1)
plot(maneuver0.t, maneuver0.eul)
title('Euler Angles')
legend('Yaw', 'Roll', 'Pitch')
subplot(2,2,2)
plot(maneuver0.t, maneuver0.omega)
title('Angular Rates')
legend('\omega_x', '\omega_y', '\omega_z')
subplot(2,2,3)
plot(maneuver0.t, maneuver0.u)
title('Torque input')
legend('u_x', 'u_y', 'u_z')
subplot(2,2,4)
hold on
plot(maneuver0.t, maneuver0.lambda_omega)
plot(maneuver0.t, maneuver0.lambda_q)
hold off
title('costates')
legend('\lambda_\omega_x', '\lambda_\omega_y', '\lambda_\omega_z', ...
    '\lambda_q_0', '\lambda_q_x', '\lambda_q_y', '\lambda_q_z')

% Final maneuver
figure(2)
clf reset

subplot(2,2,1)
plot(maneuver.t, maneuver.eul)
title('Euler Angles')
legend('Yaw', 'Roll', 'Pitch')
subplot(2,2,2)
plot(maneuver.t, maneuver.omega)
title('Angular Rates')
legend('\omega_x', '\omega_y', '\omega_z')
subplot(2,2,3)
plot(maneuver.t, maneuver.u)
title('Torque input')
legend('u_x', 'u_y', 'u_z')
subplot(2,2,4)
hold on
plot(maneuver.t, maneuver.lambda_omega)
plot(maneuver.t, maneuver.lambda_q)
hold off
title('costates')
legend('\lambda_\omega_x', '\lambda_\omega_y', '\lambda_\omega_z', ...
    '\lambda_q_0', '\lambda_q_x', '\lambda_q_y', '\lambda_q_z')