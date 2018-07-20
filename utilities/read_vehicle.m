function read_vehicle( conf )
    vehicle = load(fullfile(conf.vehicle_path,conf.vehicle_name));
    assignin('base','vehicle',vehicle.(conf.vehicle_name))
end

