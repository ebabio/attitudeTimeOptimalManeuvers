%% Eigenaxis Test

%% Setup workspace
addpath(genpath(pwd))
clearvars all
clc

%% Setup initial conditions and guess

% Initial Conditions
theta0 = 0;
omega0 = 0;
x0BCs = [theta0; omega0];

% Final conditions
thetaf = pi;
omegaf = 0;
xfBCs = [thetaf; omegaf];

% Control parameters
parameters = 1e2;

% First guess
tf = 2*sqrt(thetaf) + 0.1*randn(1);
lambda0 = -[2/tf; 1] + 0.1*randn(2,1);

%% Setup 
% 1. Trajectory ODEs
x0 = [ x0BCs; lambda0];

realOde = @(t, x) (eigenaxisODEs(t, x, parameters));
collocationOde = @(tau, x, tf) (tf*realOde(tau,x));

% 2. Boundary Conditions
collocationBCs = @(x0, xf, tf) (eigenaxisBCs(x0BCs, x0, xfBCs, xf, collocationOde, tf));

% 3. Compute initial trajectory
orbit = Orbit(x0, realOde);
orbit.integrateX0(tf);

% 4. Collocation algorithm
Nt = 20;
tau = linspace(0,1,Nt)';

solinit = bvpinit(tau, x0, tf);
solinit.x = orbit.odeSol.x;
solinit.y = orbit.odeSol.y;
solinit.parameters = tf;

options = bvpset('NMax',1e3, 'RelTol',1e-3 ,'AbsTol',1e-6); % default

%% Solution

solEigenaxis = bvp4c(collocationOde, collocationBCs, solinit, options);

%% Evaluate solution

maneuver.tf = solEigenaxis.parameters(1);

maneuver.t = solEigenaxis.x * maneuver.tf;
maneuver.x = solEigenaxis.y;

%% Display results

figure(1)
clf reset
subplot(2,1,1)
plot(orbit.t, orbit.x(1:2,:))
subplot(2,1,2)
plot(orbit.t, orbit.x(3:4,:))

figure(2)
clf reset
subplot(2,1,1)
plot(maneuver.t, maneuver.x(1:2,:))
subplot(2,1,2)
plot(maneuver.t, maneuver.x(3:4,:))

