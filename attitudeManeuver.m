%% 3d Attitude Test

%% Script mode

reuse = 1
%% Setup workspace
addpath(genpath(pwd))
if(reuse == 1 && exist('orbitFinal', 'var'))
    clearvars -except orbitFinal reuse continuation parameters fileName
else
    clearvars -except continuation fileName
end
clc

%% Define boundary conditions

% Initial / final conditions
roll = 0;
pitch  = deg2rad(0.1);
yaw = deg2rad(90);

q0 = [1 0 0 0]';
qf = eul2quat([yaw pitch roll])';

omegaBCs = zeros(3,1);

x0BCs = [q0; omegaBCs];
xfBCs = [qf; omegaBCs];

% Control parameters
if(exist('continuation', 'var') && exist('parameters', 'var') && continuation)
    continuationCoefficient = 1.03;
    parameters.lNorm = continuationCoefficient * parameters.lNorm;
    parameters.k = continuationCoefficient * parameters.k;
else
    parameters.lNorm = 3;
    parameters.k = 5e1;
end
display(['Looking for optimal control with norm ' num2str(parameters.lNorm) ':'])

%% Setup

% 1. Trajectory Integration routine
control = @(t,x) omegaControl(t, x, parameters);
shootingOde = @(t, x) (attitudeODEs(t, x, control));

orbit = Orbit(zeros(14,1), shootingOde);
orbit.odeOptions = odeset('RelTol',1e-11,'AbsTol',1e-14);

% 2. Boundary Conditions
shootingBCs = @(x0, xf) (attitudeBCs(x0BCs, x0, xfBCs, xf, shootingOde));

% 3. Shooting Algorithm
bvpShooting = OrbitShooting();
bvpShooting.epsilon = 1e-5;
bvpShooting.nIntervals = 1;
bvpShooting.k = .5;

%% Initial solution

if(~exist('orbitFinal', 'var'))
    % First guess
    thetaf = 2*acos(qf(1));
    n = qf(2:4)/norm(qf(2:4));
    tf = 2*sqrt(thetaf);
    
    lambda_q = [0; n];
    lambda_omega = -tf/4*n;
    lambda0 = [lambda_q; lambda_omega];
    
    % First actual integration
    orbit.x0 = [ x0BCs; lambda0];
    orbit.integrateX0(.98*tf);
else
    % Reuse previous solution as final guess
    orbit.x0 = orbitFinal.x0 + 0.01*randn(size(orbit.x0));
    orbit.integrateX0(.95*orbitFinal.tf);
end


%% Solve
% 1. Test for initial guess (initial point BCs should be satisfied)
f0 = shootingBCs(orbit.x0, orbit.xf);
orbit.x0(8:14)

% 2. Iterate shooting
try
    tic
    [orbitFinal, count] = bvpShooting.target(shootingBCs, orbit, 8:15);
catch
end
toc

% 3. Test for final guess (all BCs should be satisfied)
f = shootingBCs(orbitFinal.x0, orbitFinal.xf);
orbitFinal.x0(8:14)

%% Evaluate the initial solution at all times

maneuver0.t = orbit.t;
maneuver0.x = orbit.x;

maneuver0.q = maneuver0.x(1:4,:);
maneuver0.eul = quat2eul(maneuver0.q')';
maneuver0.omega = maneuver0.x(5:7,:);
maneuver0.lambda_q =  maneuver0.x(8:11,:);
maneuver0.lambda_omega =  maneuver0.x(12:14,:);

maneuver0.u = zeros(3, numel(maneuver0.t));
for i=1:numel(maneuver0.t)
    maneuver0.u(:,i) = control(maneuver0.t, maneuver0.x(:,i));
end

% plot Data
figure(1)
clf reset

subplot(2,2,1)
plot(maneuver0.t, maneuver0.eul)
title('Euler Angles')
legend('Yaw', 'Roll', 'Pitch', 'Location', 'best')
subplot(2,2,2)
plot(maneuver0.t, maneuver0.omega)
title('Angular Rates')
legend('\omega_x', '\omega_y', '\omega_z', 'Location', 'best')
subplot(2,2,3)
plot(maneuver0.t, maneuver0.u)
title('Torque input')
legend('u_x', 'u_y', 'u_z')
subplot(2,2,4)
plot(maneuver0.t, maneuver0.lambda_omega)
title('omegacostates')
legend('\lambda_\omega_x', '\lambda_\omega_y', '\lambda_\omega_z', 'Location', 'best')


%% Evaluate the final solution at all times
maneuver.t = orbitFinal.t;
maneuver.x = orbitFinal.x;

maneuver.q = maneuver.x(1:4,:);
maneuver.eul = quat2eul(maneuver.q')';
maneuver.omega = maneuver.x(5:7,:);
maneuver.lambda_q =  maneuver.x(8:11,:);
maneuver.lambda_omega =  maneuver.x(12:14,:);

maneuver.u = zeros(3, numel(maneuver.t));
for i=1:numel(maneuver.t)
    maneuver.u(:,i) = control(maneuver.t, maneuver.x(:,i));
end

%plot data
figure(2)
clf reset

subplot(2,2,1)
plot(maneuver.t, maneuver.eul)
title('Euler Angles')
legend('Yaw', 'Roll', 'Pitch', 'Location', 'best')
subplot(2,2,2)
plot(maneuver.t, maneuver.omega)
title('Angular Rates')
legend('\omega_x', '\omega_y', '\omega_z', 'Location', 'best')
subplot(2,2,3)
plot(maneuver.t, maneuver.u)
title('Torque input')
legend('u_x', 'u_y', 'u_z', 'Location', 'best')
subplot(2,2,4)
plot(maneuver.t, maneuver.lambda_omega)
title('omegacostates')
legend('\lambda_\omega_x', '\lambda_\omega_y', '\lambda_\omega_z', 'Location', 'best')