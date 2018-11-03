clearvars;
clc;    

addpath(genpath(pwd));

%% Todo
% Vehicle create F and M ode fun inside its class
% create read vehicle function
% plot inside vehicle class
% 

%% User Settings

conf.solver         = 'ode1';                       % 'ode1', 'ode2', 'ode45'
conf.time_step      = 0.001;                        % simulation time step
conf.vehicle_path   = fullfile(pwd,'vehicles');     % path of the vehicle mat folder
conf.vehicle_name   = 'veh.mat';                    % vehicle mat file name

conf.maneuver = 'SWD';



%% Simulate

% Load vehicle
fprintf('Loading Kinematics...')
veh = open(fullfile(conf.vehicle_path, conf.vehicle_name));
veh = veh.veh;
fprintf('OK\n')

% Plot kinematics analysis 
if false
    veh.plotKinAnalysis;
    veh.sus_fl.compare_K_analytic_numeric('cWH_3');
    veh.sus_fl.compare_K_analytic_numeric('bDS');
end

% Plot kinematics animation
if false
    veh.plotKinAnimation();
end

% Calculate Vehicle Dynamics
veh.calculateDynamics;

% veh.sim()






