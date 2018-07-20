% function [ maneuver_output ] = man_SWD(conf)

    % Settings
    
    conf_SWD.stf_frequency  = 0.7;        % Front steering wheel angle frequency [Hz]
    conf_SWD.stf_hold_time  = 0.5;        % Hold time in the sinus valey [s] 
    conf_SWD.stf_amplitude  = 120;        % Front steering wheel angle amplitude [º]       

    
    
    %% Configure simulation
    % https://de.mathworks.com/help/simulink/slref/model-parameters.html

    set_param('top_model',...
                'SolverType','Fixed-step',...
                'Solver', conf.solver,...
                'FixedStep',num2str(conf.time_step),...
                'MaxOrder','3',...
                'StartTime','0',...
                'StopTime','20')
            
    driver_choice = Simulink.Parameter(1);          % 1=Joystick 2=Signalbuilder 3=Driver

    %% Set initial conditions

    veh_veh_vx_0 = 0;
    veh_veh_sx_0 = 0;
    
    %% Start simulation
    
    maneuver_output = sim('top_model','timeout',100)

% end

