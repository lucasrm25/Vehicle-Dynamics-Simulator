%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE A FORMULA SAE 2017 VEHICLE containing the kinematics and dynamics
% 
% Inputs: kinematics hardpoints
% Output: create a vehicle of class 'veh' and save into vehicles/veh.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc;

save_veh = true;
save_path = 'vehicles/veh.mat';

% if isempty(gcp('nocreate')), parpool(4); end


%% Kinematics data

env.grav = [0 0 -9.8]';         % [m/s2]

% Vehicle Parameters
DNA.track               = 1.200;                % [m]
DNA.wheelbase           = 1.55;                 % [m]
DNA.DIN_CO              = 1e-3*[1528 0 224.7]'; % [m]
DNA.CG_s                = 1e-3*[1000 0 300]';   % [m] CG position
DNA.susMass             = 180;                  % [kg]
DNA.CG_I                = diag([100 100 100]);  % [kg.m2]

tyr.unloadedRadius    = 8*0.0254;     % [m]
tyr.aspectRatio       = 45;           % -
tyr.width             = 0.220;        % [m]
tyr.compressionLength = 0.08;         % [m]
tyr.stiffness         = 500000;       % [N/m] 

% Mass
susFL.unsMass = 2;
% setup
susFL.setup.toe                 = 1*pi/180;     % [rad]
susFL.setup.camber              = -2*pi/180;    % [rad]
susFL.setup.kspring             = 200000;       % [N/m]
% Hard points
susFL.init.cWH_0 = [DNA.DIN_CO(1) DNA.track/2+susFL.setup.camber*tyr.unloadedRadius DNA.DIN_CO(3)]';
susFL.init.UACH1_0 = 1e-3*[1659.20 260.000 250.000]';
susFL.init.UACH2_0 = 1e-3*[1409.50 260.000 230.000]';
susFL.init.LACH1_0 = 1e-3*[1648.00 170.000 100.000]';
susFL.init.LACH2_0 = 1e-3*[1406.50 170.000 110.000]';
susFL.init.UAUR_0  = 1e-3*[1513.00 521.000 307.000]';
susFL.init.LAUR_0  = 1e-3*[1525.00 540.000 134.000]';
susFL.init.TRUR_0 = 1e-3*[1576.00 564.000 165.000]';
susFL.init.SRTR_0 = 1e-3*[1618.00 179.000 123.120]';
susFL.init.DSCH_0 = 1e-3*[1473.00 260.000 320.000]';
susFL.init.RKDS_0 = 1e-3*[1484.68 324.508 157.743]';
susFL.init.RKPR_0 = 1e-3*[1478.43 290.000 85.0000]';
susFL.init.RKCH_0 = 1e-3*[1460.33 190.000 100.000]';
susFL.init.PRUA_0 = 1e-3*[1515.00 492.000 275.000]';
susFL.init.RKAL_0 = 1e-3*[1472.09 255.000 130.000]';
susFL.init.ALAB_0 = 1e-3*[1472.09 255.000 70.0000]';
susFL.init.ALAR_0 = 1e-3*[1392.09 255.000 70.0000]';


% Mass
susFR.unsMass = 2;
% setup
susFR.setup.toe                 = -1*pi/180;     % [rad]
susFR.setup.camber              = 2*pi/180;    % [rad]
susFR.setup.kspring             = 200000;       % [N/m]
% Hard points
susFR.init = susFL.init;
names = fieldnames(susFR.init);
for i=1:numel( names )
    susFR.init.(names{i}) = susFR.init.(names{i}) .* [1 -1 1]';
end



% Mass
susRL.unsMass = 2;
% setup
susRL.setup.toe                 = -1*pi/180;      % [rad]
susRL.setup.camber              = -2*pi/180;      % [rad]
susRL.setup.kspring             = 200000;         % [N/m]
% Hard points
susRL.init.cWH_0 = 1e-3*[8.000 550.000 220.000]';
susRL.init.UACH1_0 = 1e-3*[350.000 325.000 269.000]';
susRL.init.UACH2_0 = 1e-3*[160.200 265.000 240.000]';
susRL.init.LACH1_0 = 1e-3*[350.000 315.000 150.000]';
susRL.init.LACH2_0 = 1e-3*[159.600 230.000 131.500]';
susRL.init.UAUR_0  = 1e-3*[15.0000 520.000 305.000]';
susRL.init.LAUR_0  = 1e-3*[3.00000 535.000 142.000]';
susRL.init.TRUR_0 = 1e-3*[-36.0000 520.000 301.879]';
susRL.init.SRTR_0 = 1e-3*[126.000 285.000 246.000]';
susRL.init.DSCH_0 = 1e-3*[215.000 270.000 355.000]';
susRL.init.RKDS_0 = 1e-3*[178.889 317.733 190.554]';
susRL.init.RKPR_0 = 1e-3*[183.226 312.000 130.000]';
susRL.init.RKCH_0 = 1e-3*[208.275 257.341 141.224]';
susRL.init.PRUA_0 = 1e-3*[41.0000 500.000 275.800]';
susRL.init.RKAL_0 = 1e-3*[202.896 286.000 165.000]';
susRL.init.ALAB_0 = 1e-3*[202.896 286.000 105.000]';
susRL.init.ALAR_0 = 1e-3*[96.8960 286.000 105.000]';

% Mass
susRR.unsMass = 2;
% setup
susRR.setup.toe                 = 1*pi/180;      % [rad]
susRR.setup.camber              = 2*pi/180;      % [rad]
susRR.setup.kspring             = 200000;        % [N/m]
% Hard points
susRR.init = susRL.init;
names = fieldnames(susRR.init);
for i=1:numel( names )
    susRR.init.(names{i}) = susRR.init.(names{i}) .* [1 -1 1]';
end


%% Run Kinematics and create vehicle

WH_stroke = [-2 2]*0.0254;
SR_stroke = [-0.02 0.02];
WH_steps = 15;
SR_steps = 10;

veh_Tyr = veh_tyr(tyr.unloadedRadius, tyr.aspectRatio, tyr.width, tyr.compressionLength, tyr.stiffness);


fprintf('Computing Kinematics - Parallel batch mode...')
j = {batch('veh_sus', 1, {susFL.init, susFL.setup, susFL.unsMass, veh_Tyr, WH_stroke, WH_steps, SR_stroke, SR_steps});
     batch('veh_sus', 1, {susFR.init, susFR.setup, susFR.unsMass, veh_Tyr, WH_stroke, WH_steps, SR_stroke, SR_steps});
     batch('veh_sus', 1, {susRL.init, susRL.setup, susRL.unsMass, veh_Tyr, WH_stroke, WH_steps, 0.1*SR_stroke, 2});
     batch('veh_sus', 1, {susRR.init, susRR.setup, susRR.unsMass, veh_Tyr, WH_stroke, WH_steps, 0.1*SR_stroke, 2})   };
for ji =1:numel(j)
    wait(j{ji});
end
veh_sus_fl = fetchOutputs(j{1}); veh_sus_fl = veh_sus_fl{1};
veh_sus_fr = fetchOutputs(j{2}); veh_sus_fr = veh_sus_fr{1};
veh_sus_rl = fetchOutputs(j{3}); veh_sus_rl = veh_sus_rl{1};
veh_sus_rr = fetchOutputs(j{4}); veh_sus_rr = veh_sus_rr{1};
fprintf('OK\n')

% veh_sus_fl = veh_sus(susFL.init, susFL.setup, susFL.unsMass, veh_Tyr, WH_stroke, WH_steps, SR_stroke, SR_steps);
% veh_sus_fr = veh_sus(susFR.init, susFR.setup, susFR.unsMass, veh_Tyr, WH_stroke, WH_steps, SR_stroke, SR_steps);
% veh_sus_rl = veh_sus(susRL.init, susRL.setup, susRL.unsMass, veh_Tyr, WH_stroke, WH_steps, 0.1*SR_stroke, 2);
% veh_sus_rr = veh_sus(susRR.init, susRR.setup, susRR.unsMass, veh_Tyr, WH_stroke, WH_steps, 0.1*SR_stroke, 2);


veh = veh(DNA, veh_sus_fl, veh_sus_fr, veh_sus_rl, veh_sus_rr);


%% Save to file

if save_veh
    fprintf('Saving Kinematics...')
    save(save_path, 'veh');
    fprintf('OK\n')
end