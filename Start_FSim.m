clearvars;
clc;

addpath(genpath(pwd));

%% User Settings

load_kinematics = true;
save_kinematics = ~load_kinematics;

conf.solver         = 'ode1';                       % 'ode1', 'ode2', 'ode45'
conf.time_step      = 0.001;                        % simulation time step
conf.vehicle_path   = fullfile(pwd,'vehicles');     % path of the vehicle mat file
conf.vehicle_name   = 'FSAE2017';                   % vehicle mat file name

conf.maneuver = 'SWD';


%% Read Vehicle


if ~load_kinematics
    % % % % % % % % % % % % % % read_vehicle(conf)
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
end

%% Create Vehicle Kinematics

if ~load_kinematics
    WH_stroke = [-2 2]*0.0254;
    SR_stroke = [-0.02 0.02];
    WH_steps = 15;
    SR_steps = 10;

    % Create Tyre
    veh_tyr = veh_tyr(tyr.unloadedRadius, tyr.aspectRatio, tyr.width, tyr.compressionLength, tyr.stiffness);

    % Create Suspension
    if isempty(gcp('nocreate'))
        parpool(4)
    end

    veh_sus_fl = veh_sus(susFL.init, susFL.setup, susFL.unsMass, veh_tyr, WH_stroke, WH_steps, SR_stroke, SR_steps);
    veh_sus_fr = veh_sus_fl.mirror();
    veh_sus_rl = veh_sus(susRL.init, susRL.setup, susRL.unsMass, veh_tyr, WH_stroke, WH_steps, [0 0], 1);
    veh_sus_rr = veh_sus_rl.mirror();

    % veh_sus_fl.generate_Kmatrix_numeric;

    % Create Vehicle
    veh = veh(DNA, veh_sus_fl, veh_sus_fr, veh_sus_rl, veh_sus_rr);
    
    if save_kinematics
        fprintf('Saving Kinematics...')
        save('vehicles/veh.mat', 'veh');
        fprintf('OK\n')
    end
else
    fprintf('Loading Kinematics...')
    veh = open('vehicles/veh.mat');
    veh = veh.veh;
    fprintf('OK\n')
end

% veh.plotKinAnalysis;
if false
    veh.sus_fl.compare_K_analytic_numeric('cWH_3');
    veh.sus_fl.compare_K_analytic_numeric('bDS');
end

%% Calculate Vehicle Dynamics

veh.calculateDynamics;
% veh.sim()


% Plot K matrices to file !!!
% aa = func2str(veh_sus_fl.Kmatrix.TYRO_3);
% fileID = fopen(fullfile(pwd,'test.txt'),'w');
% fprintf(fileID,'%s\n',aa);
% fclose(fileID);



%% Plot

movie = false;

fig = figure('Color','w');
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')

axis equal
xlim([-0.2 2])
ylim([-0.9 0.9])
zlim([-0.1 0.6])

% view(3)
az = 210;
el = 45;
view(az, el);

bSR_slider = 0;
FLAUR_3_slider = 0;
RLAUR_3_slider = 0;
hscrollbar_bSR    = uicontrol('style','slider','units','normalized','position',[0.05 0 0.90 .05],'callback', @(src,evt) assignin('base','bSR_slider',src.Value) );
hscrollbar_FLAUR  = uicontrol('style','slider','units','normalized','position',[0 .05 .05 0.95],'callback',@(src,evt) assignin('base','FLAUR_3_slider',src.Value) );
hscrollbar_RLAUR  = uicontrol('style','slider','units','normalized','position',[0.95 .05 .05 0.95],'callback',@(src,evt) assignin('base','RLAUR_3_slider',src.Value) );


WH_stroke = [-2 2]*0.0254;
SR_stroke = [-0.02 0.02];

loops = 40;
F(loops) = struct('cdata',[],'colormap',[]);
init=true;
% while true
for i = 1:loops
    FXq = veh.sus_fl.init.LAUR_0(3)+WH_stroke(1) + (WH_stroke(2)-WH_stroke(1))*(sin(i/loops*2*pi)/2+0.5); %*FLAUR_3_slider;
    FYq = -(SR_stroke(1) + (SR_stroke(2)-SR_stroke(1))*bSR_slider);
    RXq = veh.sus_rl.init.LAUR_0(3)+WH_stroke(1) + (WH_stroke(2)-WH_stroke(1))*(sin(i/loops*2*pi+pi/2)/2+0.5); %*RLAUR_3_slider;
    RYq = 0;
    
    q_tm1 = [0 0 0   0 0 0  FXq FXq RXq RXq   FYq RYq ];
    
    veh.sus_fl.q_val = q_tm1([7 11])';
    veh.sus_fr.q_val = q_tm1([8 11])';
    veh.sus_rl.q_val = q_tm1([9 12])';
    veh.sus_rr.q_val = q_tm1([10 12])';

    veh.sus_fl.plot(fig,init)
    veh.sus_fr.plot(fig,init)
	veh.sus_rl.plot(fig,init)
    veh.sus_rr.plot(fig,init)
   
%     axis vis3d
    drawnow();
    if movie
        F(i) = getframe(gcf);
    end
    init = false;
end
% end

if movie
    fig = figure;
    movie(fig,F,2);
end

%% Start Maneuver

% man_SWD;
% man_Joystick(conf);
% man_Signalbuilder(conf);

