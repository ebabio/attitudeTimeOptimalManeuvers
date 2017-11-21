function [x0, tf] = firstAttitudeGuess(qf)

%% Get axis info
thetaf = 2*acos(qf(1));
n = qf(2:4)/norm(qf(2:4));

%% Setup initial conditions and guess

% Initial Conditions
theta0 = 0;
omega0 = 0;
x0BCs = [theta0; omega0];

% Final conditions
% theta f has just been defined above
omegaf = 0;
xfBCs = [thetaf; omegaf];

% Control parameters
parameters = 1e2;

% First guess
tf = 2*sqrt(thetaf);
lambda0 = -[2/tf; 1];

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
bvpShooting.nIntervals = 2;

%% Solution

% First integrate initial guess
orbit.integrateX0(tf);

% Solve 
[orbit, count] = bvpShooting.target(shootingBCs, orbit, [3 4 5]);

%% Transform to 3d attitude problem
tf = orbit.tf;

q = [cos(theta/2); n*sin(theta/2)];
omega = zeros(3,1);
lambda_omega = -n;
lambda_q = [0; 2*lambda_omega]/tf;

x0 = [q; omega; lambda_q; lambda_omega;];

