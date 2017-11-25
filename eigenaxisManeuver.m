%% Eigenaxis Test

%% Setup workspace
addpath(genpath(pwd))
clear
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
tf = 2*sqrt(thetaf) + 0.3*randn(1);
lambda0 = -[2/tf; 1] + 0.3*randn(2,1);

%% Setup

% 1. Trajectory Integration routine

x0 = [ x0BCs; lambda0];
shootingOde = @(t, x) (eigenaxisODEs(t, x, parameters));

orbit = Orbit(x0, shootingOde);
orbit.odeOptions = odeset('RelTol',1e-12,'AbsTol',1e-15);

% 2. Boundary Conditions
shootingBCs = @(x0, xf) (eigenaxisBCs(x0BCs, x0, xfBCs, xf, shootingOde));

% 3. Optimizer setup
fObjective = @(x) x(end);

struct.nIntervals = 3;
struct.setupOrbit = orbit;
gConstraints = @(x) ( nonlconWrapper(x, shootingBCs, struct) );

lb = [x0BCs - eps*ones(size(x0BCs)); -100*ones(size(lambda0)); xfBCs - eps*ones(size(x0BCs)); -100*ones(size(lambda0)); 0];
ub = [x0BCs + eps*ones(size(x0BCs)); +100*ones(size(lambda0)); xfBCs + eps*ones(size(x0BCs)); +100*ones(size(lambda0)); 1.5*tf];
%% Get initial guess
% Compute one whole trajectory
orbit.integrateX0(tf);

% Setup independent shootings
dim = numel(x0);
tfVector = orbit.tf * (1:struct.nIntervals)./struct.nIntervals;
for i=1:struct.nIntervals
    if(i==1)
        x0Opt = orbit.x0(1:dim);
    else
        index0 = (i-1)*dim;
        x0Opt(index0+1:index0+dim) = deval(orbit.odeSol, tfVector(i-1));
    end
end
x0Opt(end+1) = tf;

% Display data
figure(1)
clf reset
subplot(2,1,1)
plot(orbit.t, orbit.x(1:2,:))
subplot(2,1,2)
plot(orbit.t, orbit.x(3:4,:))

%% Test2: achieve targeting

% 1. Test for initial guess (initial point BCs should be satisfied)
f = shootingBCs(orbit.x0, orbit.xf)

% 2. Iterate shooting
tic
[xf, tf, exitflaf, output] = fmincon(fObjective, x0Opt, [], [], [], [], [], [], gConstraints)
toc

% 3. Test for final guess (all BCs should be satisfied)
orbit.x0 = xf(1:4);
orbit.integrateX0(tf);
f = shootingBCs(orbit.x0, orbit.xf)

figure(2)
clf reset
subplot(2,1,1)
plot(orbit.t, orbit.x(1:2,:))
subplot(2,1,2)
plot(orbit.t, orbit.x(3:4,:))

