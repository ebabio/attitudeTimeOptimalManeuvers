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
thetaf = .1;
omegaf = 0;
xfBCs = [thetaf; omegaf];

% Control parameters
parameters = 1e2;

% First guess
tf = 2*sqrt(thetaf) + 0.1*randn(1);
lambda0 = -[2/tf; 1] + 0.1*randn(2,1);

%% Setup 

% 1. Trajectory Integration routine

x0 = [ x0BCs; lambda0];
shootingOde = @(t, x) (eigenaxisODEs(t, x, parameters));

orbit = Orbit(x0, shootingOde);
orbit.odeOptions = odeset('RelTol',1e-12,'AbsTol',1e-15);

% 2. Boundary Conditions
shootingBCs = @(x0, xf) (eigenaxisBCs(x0BCs, x0, xfBCs, xf, shootingOde));

% 3. Shooting Algorithm
bvpShooting = OrbitShooting();
bvpShooting.epsilon = 1e-8;
bvpShooting.nIntervals = 1;

%% Test1: integration

orbit.integrateX0(tf);

figure(1)
clf reset
subplot(2,1,1)
plot(orbit.t, orbit.x(1:2,:))
subplot(2,1,2)
plot(orbit.t, orbit.x(3:4,:))

%% Test2: achieve targeting

% 1. Test for initial guess (initial point BCs should be satisfied)
f = shootingBCs(orbit.x0, orbit.xf);

% 2. Iterate shooting
tic
[orbit, count] = bvpShooting.target(shootingBCs, orbit, [3 4 5]);
toc

% 3. Test for final guess (all BCs should be satisfied)
f = shootingBCs(orbit.x0, orbit.xf)

figure(2)
clf reset
subplot(2,1,1)
plot(orbit.t, orbit.x(1:2,:))
subplot(2,1,2)
plot(orbit.t, orbit.x(3:4,:))

